#include "MRCPP/MWFunctions"
#include "MRCPP/Parallel"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include <memory>

const auto min_scale = -4;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-5;

const auto D = 3;

int main(int argc, char **argv) {
    auto tot_t = mrcpp::Timer();

    // Initialize MPI

    auto wcomm = MPI_Comm();
    auto scomm = MPI_Comm();

    auto wrank = int();
    auto wsize = int();
    auto srank = int();
    auto ssize = int();

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
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel, wrank, wsize);
    mrcpp::print::environment(0);
    mrcpp::print::header(0, "Shared memory MPI");

    // Constructing world box
    auto corner = std::array<int, D>{-1, -1, -1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto world = mrcpp::BoundingBox<D>(min_scale, corner, boxes);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);

    // Defining analytic function
    auto f = [](const mrcpp::Coord<D> &r) -> double {
        const auto beta = 100.0;
        const auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return alpha * std::exp(-beta * R * R);
    };

    // Initialize a shared memory tree, max 100MB
    auto shared_mem = new mrcpp::SharedMemory(scomm, 100);
    mrcpp::FunctionTree<D> f_tree(MRA, shared_mem);

    // Only first rank projects
    auto frank = 0;
    if (srank == frank) mrcpp::project<D>(prec, f_tree, f);
    mrcpp::share_tree(f_tree, frank, 0, scomm);

    { // Print data after share
        auto integral = f_tree.integrate();
        auto sq_norm = f_tree.getSquareNorm();
        mrcpp::print::value(0, "Integral", integral);
        mrcpp::print::value(0, "Square norm", sq_norm);
    }

    // Last rank rescales the tree
    auto lrank = ssize - 1;
    if (srank == lrank) f_tree.rescale(2.0);
    mrcpp::share_tree(f_tree, lrank, 0, scomm);

    { // Print data after rescale
        auto integral = f_tree.integrate();
        auto sq_norm = f_tree.getSquareNorm();
        mrcpp::print::value(0, "Integral", integral);
        mrcpp::print::value(0, "Square norm", sq_norm);
    }

    // Must be deleted before MPI_Finalize
    delete shared_mem;

    // Finalize MPI
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    mrcpp::print::footer(0, tot_t, 2);

    return 0;
}
