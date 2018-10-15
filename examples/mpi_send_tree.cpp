#include "MRCPP/MWFunctions"
#include "MRCPP/Parallel"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const auto min_scale = -4;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-5;

const auto D = 3; // Dimensions

int main(int argc, char **argv) {
    mrcpp::Timer tot_t;

    // Initialize MPI
    auto comm = MPI_Comm();
    auto wrank = int();
    auto wsize = int();

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
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel, wrank, wsize);
    mrcpp::Printer::printEnvironment();
    mrcpp::Printer::printHeader(0, "Blocking communication");

    // Constructing world box
    auto corner = std::array<int, D>{-1,-1,-1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto world = mrcpp::BoundingBox<D>(min_scale, corner, boxes);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);

    // Defining analytic function
    auto f = [] (const mrcpp::Coord<D> &r) -> double {
        const auto beta = 100.0;
        const auto alpha = std::pow(beta/mrcpp::pi, 3.0/2.0);
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return alpha*std::exp(-beta*R*R);
    };

    // All ranks define the function
    mrcpp::FunctionTree<D> f_tree(MRA);

    // Only rank 0 projects the function
    if (wrank == 0) mrcpp::project<D>(prec, f_tree, f);

    {   // Print data before send
        auto integral = f_tree.integrate();
        auto sq_norm = f_tree.getSquareNorm();
        mrcpp::Printer::printDouble(0, "Integral", integral);
        mrcpp::Printer::printDouble(0, "Square norm", sq_norm);
    }

    auto send_t = mrcpp::Timer();
    // Send function to all other ranks
    auto src = 0;
    for (auto dst = 0; dst < wsize; dst++) {
        if (dst == src) continue;
        auto tag = 11111*dst; // Unique tag for each communication
        if (wrank == src) mrcpp::send_tree(f_tree, dst, tag, comm);
        if (wrank == dst) mrcpp::recv_tree(f_tree, src, tag, comm);
    }
    send_t.stop();
    mrcpp::Printer::printSeparator(0, ' ');
    mrcpp::Printer::printTime(0, "Time sending tree", send_t);
    mrcpp::Printer::printSeparator(0, ' ');

    {   // Print data after send
        auto integral = f_tree.integrate();
        auto sq_norm = f_tree.getSquareNorm();
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
