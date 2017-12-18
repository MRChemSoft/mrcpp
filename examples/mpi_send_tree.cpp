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
    int wrank, wsize;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &wrank);
    MPI_Comm_size(comm, &wsize);
#else
    comm = 0;
    wrank = 0;
    wsize = 1;
#endif

    // Initialize printing
    int printlevel = 0;
    mrcpp::Printer::init(printlevel, wrank, wsize);
    mrcpp::Printer::printEnvironment();
    mrcpp::Printer::printHeader(0, "Blocking communication");

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    mrcpp::BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    mrcpp::InterpolatingBasis basis(order);
    mrcpp::MultiResolutionAnalysis<3> MRA(world, basis, max_depth);

    // Defining analytic function
    auto f = [] (const double *r) -> double {
        const double beta = 100.0;
        const double alpha = pow(beta/mrcpp::pi, 3.0/2.0);
        double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return alpha*exp(-beta*R*R);
    };

    // All ranks define the function
    mrcpp::FunctionTree<3> f_tree(MRA);

    // Only rank 0 projects the function
    if (wrank == 0) mrcpp::project(prec, f_tree, f);

    {   // Print data before send
        double integral = f_tree.integrate();
        double sq_norm = f_tree.getSquareNorm();
        mrcpp::Printer::printDouble(0, "Integral", integral);
        mrcpp::Printer::printDouble(0, "Square norm", sq_norm);
    }

    mrcpp::Timer send_t;
    // Send function to all other ranks
    int src = 0;
    for (int dst = 0; dst < wsize; dst++) {
        if (dst == src) continue;
        int tag = 11111*dst; // Unique tag for each communication
        if (wrank == src) mrcpp::send_tree(f_tree, dst, tag, comm);
        if (wrank == dst) mrcpp::recv_tree(f_tree, src, tag, comm);
    }
    send_t.stop();
    mrcpp::Printer::printSeparator(0, ' ');
    mrcpp::Printer::printTime(0, "Time sending tree", send_t);
    mrcpp::Printer::printSeparator(0, ' ');

    {   // Print data after send
        double integral = f_tree.integrate();
        double sq_norm = f_tree.getSquareNorm();
        mrcpp::Printer::printDouble(0, "Integral", integral);
        mrcpp::Printer::printDouble(0, "Square norm", sq_norm);
    }

    // Finalize MPI
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    tot_t.stop();
    mrcpp::Printer::printFooter(0, tot_t, 2);

    return 0;
}

