#ifndef HYPERPLANEFINDER_MATH_HPP
#define HYPERPLANEFINDER_MATH_HPP


#include <type_traits>

namespace math {
	template <typename T, typename U>
	inline constexpr T pow(T base, U exponent) {
		static_assert(std::is_integral_v<U>);
		return exponent == 0 ? 1 : base * pow(base, exponent - 1);
	}
}


#endif //HYPERPLANEFINDER_MATH_HPP
