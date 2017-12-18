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
    MPI_Comm wcomm, scomm;
    int wrank, wsize, srank, ssize;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);

    wcomm = MPI_COMM_WORLD;
    MPI_Comm_rank(wcomm, &wrank);
    MPI_Comm_size(wcomm, &wsize);

    MPI_Comm_split_type(wcomm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &scomm);
    MPI_Comm_rank(scomm, &srank);
    MPI_Comm_size(scomm, &ssize);
#else
    wcomm = 0;
    wrank = 0;
    wsize = 1;

    scomm = 0;
    srank = 0;
    ssize = 1;
#endif

    // Initialize printing
    int printlevel = 0;
    mrcpp::Printer::init(printlevel, wrank, wsize);
    mrcpp::Printer::printEnvironment();
    mrcpp::Printer::printHeader(0, "Shared memory MPI");

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

    // Initialize a shared memory tree, max 100MB
    mrcpp::SharedMemory *shared_mem = new mrcpp::SharedMemory(scomm, 100);
    mrcpp::FunctionTree<3> f_tree(MRA, shared_mem);

    // Only first rank projects
    int frank = 0;
    if (srank == frank) mrcpp::project(prec, f_tree, f);
    mrcpp::share_tree(f_tree, frank, 0, scomm);

    {   // Print data after share
        double integral = f_tree.integrate();
        double sq_norm = f_tree.getSquareNorm();
        mrcpp::Printer::printDouble(0, "Integral", integral);
        mrcpp::Printer::printDouble(0, "Square norm", sq_norm);
    }

    // Last rank rescales the tree
    int lrank = ssize - 1;
    if (srank == lrank) f_tree *= 2.0;
    mrcpp::share_tree(f_tree, lrank, 0, scomm);

    {   // Print data after rescale
        double integral = f_tree.integrate();
        double sq_norm = f_tree.getSquareNorm();
        mrcpp::Printer::printDouble(0, "Integral", integral);
        mrcpp::Printer::printDouble(0, "Square norm", sq_norm);
    }

    // Must be deleted before MPI_Finalize
    delete shared_mem;

    // Finalize MPI
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    tot_t.stop();
    mrcpp::Printer::printFooter(0, tot_t, 2);

    return 0;
}
