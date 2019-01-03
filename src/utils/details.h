namespace mrcpp {
namespace details {

template<int D>
bool are_all_equal(const std::array<double, D> &exponent);

template<typename T, int D>
std::array<T, D> convert_to_std_array(T *arr);
} // namespace details
} // namespace utils
