#include <tuple>

#include "MRCPP/MWFunctions"
#include "MRCPP/Parallel"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const int min_scale = -4;
const int max_depth = 25;

const int order = 7;
const double prec = 1.0e-5;

int main(int argc, char **argv) {
    mrcpp::Timer tot_t;

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
    mrcpp::Printer::init(printlevel, wrank, wsize);
    mrcpp::Printer::printEnvironment();
    mrcpp::Printer::printHeader(0, "Non-blocking communication");

    // Constructing world box
    const auto corner = std::array<int, 3>{-1, -1, -1};
    const auto boxes = std::array<int, 3>{2, 2, 2};
    mrcpp::BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    mrcpp::InterpolatingBasis basis(order);
    mrcpp::MultiResolutionAnalysis<3> MRA(world, basis, max_depth);

    // Projecting vector of functions
    int nFuncs = 5;
    mrcpp::FunctionTreeVector<3> f_vec;
    for (int i = 0; i < nFuncs; i++) {
        auto f = [i](const mrcpp::Coord<3> &r) -> double {
            const double beta = 1.0 * i * i;
            double R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
            return std::exp(-beta * R * R);
        };
        mrcpp::FunctionTree<3> *tree = new mrcpp::FunctionTree<3>(MRA);
        if (i % wsize == wrank) {
            mrcpp::project<3>(prec, *tree, f);
            tree->normalize();
        }
        f_vec.push_back(std::make_tuple(1.0, tree));
    }

    std::vector<MPI_Request> requests;
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(nFuncs, nFuncs);
    for (int j = 0; j < f_vec.size(); j++) {
        int dst = j % wsize;
        mrcpp::FunctionTree<3> &f_j = get_func(f_vec, j);
        for (int i = 0; i < f_vec.size(); i++) {
            int src = i % wsize;
            mrcpp::FunctionTree<3> &f_i = get_func(f_vec, i);
            if (src != dst) {
                int tag = 1000000 * dst;
                MPI_Request req = req_null;
                if (wrank == src) mrcpp::isend_tree(f_i, dst, tag, comm, &req);
                if (wrank == dst) mrcpp::recv_tree(f_i, src, tag, comm);
                requests.push_back(req);
            }
            // Compute my column(s) of the overlap matrix
            if (wrank == dst) S(i, j) = mrcpp::dot(f_i, f_j);
        }
    }

#ifdef HAVE_MPI
    for (int i = 0; i < requests.size(); i++) { MPI_Wait(&requests[i], MPI_STATUS_IGNORE); }
#endif

    // Delete all trees
    clear(f_vec, true);

    // Print overlap matrix
    mrcpp::Printer::setPrecision(5);
    println(0, S);

    // Finalize MPI
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    tot_t.stop();
    mrcpp::Printer::printFooter(0, tot_t, 2);

    return 0;
}
