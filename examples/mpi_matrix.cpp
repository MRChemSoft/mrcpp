#include "MRCPP/MWFunctions"
#include "MRCPP/Parallel"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const int min_scale = -4;
const int max_scale = 20;

const int order = 7;
const double prec = 1.0e-5;

int main(int argc, char **argv) {
    Timer tot_t;

    // Initialize MPI
    MPI_Comm comm;
    MPI_Request req_null;
    int wrank, wsize;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;
    req_null = MPI_REQUEST_NULL;
    MPI_Comm_rank(comm, &wrank);
    MPI_Comm_size(comm, &wsize);
#else
    comm = 0;
    req_null = 0;
    wrank = 0;
    wsize = 1;
#endif

    // Initialize printing
    int printlevel = 0;
    Printer::init(printlevel, wrank, wsize);
    Printer::printEnvironment();
    Printer::printHeader(0, "Non-blocking communication");

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(world, basis);

    // Setting up projector
    MWProjector<3> project(prec, max_scale);

    // Projecting vector of functions
    int nFuncs = 5;
    FunctionTreeVector<3> f_vec;
    for (int i = 0; i < nFuncs; i++) {
        auto f = [i] (const double *r) -> double {
            const double beta = 1.0*i*i;
            const double r_0[3] = {0.0, 0.0, 0.0};
            double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
            return exp(-beta*R*R);
        };
        FunctionTree<3> *tree = new FunctionTree<3>(MRA);
        if (i%wsize == wrank) {
            project(*tree, f);
            tree->normalize();
        }
        f_vec.push_back(tree);
    }

    std::vector<MPI_Request> requests;
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(nFuncs, nFuncs);
    for (int j = 0; j < f_vec.size(); j++) {
        int dst = j%wsize;
        FunctionTree<3> &f_j = *f_vec[j];
        for (int i = 0; i < f_vec.size(); i++) {
            int src = i%wsize;
            FunctionTree<3> &f_i = *f_vec[i];
            if (src != dst) {
                int tag = 1000000*dst;
                MPI_Request req = req_null;
                if (wrank == src) isend_tree(f_i, dst, tag, comm, &req);
                if (wrank == dst) recv_tree(f_i, src, tag, comm);
                requests.push_back(req);
            }
            // Compute my column(s) of the overlap matrix
            if (wrank == dst) S(i,j) = f_i.dot(f_j);
        }
    }

#ifdef HAVE_MPI
    for (int i = 0; i < requests.size(); i++) {
        MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
    }
#endif

    // Delete all trees
    f_vec.clear(true);

    // Print overlap matrix
    Printer::setPrecision(5);
    println(0, S);

    // Finalize MPI
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    tot_t.stop();
    Printer::printFooter(0, tot_t, 2);

    return 0;
}

