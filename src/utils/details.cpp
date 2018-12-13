#include <algorithm>
#include <array>

namespace mrcpp {
namespace details {


/** @brief checks if all elements of an array of doubles are equal */
template<int D>
bool are_all_equal(const std::array<double, D> &exponent) {
        return std::all_of(exponent.begin(), exponent.end(),
               [ex = std::begin(exponent)](double i) {return i == *ex; });
    }

template bool are_all_equal<1>(const std::array<double, 1> &exponent);
template bool are_all_equal<2>(const std::array<double, 2> &exponent);
template bool are_all_equal<3>(const std::array<double, 3> &exponent);
} // namespace details
} // namespace mrcpp
