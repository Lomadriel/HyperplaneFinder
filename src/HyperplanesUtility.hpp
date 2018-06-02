#ifndef HYPERPLANEFINDER_HYPERPLANESUTILITY_HPP
#define HYPERPLANEFINDER_HYPERPLANESUTILITY_HPP


#include <bitset>
#include <vector>
#include <tuple>

#include "PermutationGenerator.hpp"
#include "index_repetition.hpp"
#include "math.hpp"
#include "impossible.hpp"

// Declarations
namespace segre {

	/*------------------------------------------------------------------------*//**
	 * @brief      Convert an hyperplane representation from a std::bitset to a
	 *             std::vector.
	 *
	 * @details    In the std::bitset representation, each bit represent a point
	 *             of the hyperplane, included if equal true, excluded if false.
	 *
	 *             In the std::vector representation, the vector contain the
	 *             index of all included points in increasing order.
	 *
	 * @param[in]  hyperplane  The hyperplane to convert
	 *
	 * @tparam     NbrPoints   Number of points of the geometry
	 *
	 * @return     The hyperplane converted to std::vector
	 */
	template<size_t NbrPoints>
	std::vector<unsigned int> bitsetToVector(const std::bitset<NbrPoints>& hyperplane);

	/*------------------------------------------------------------------------*//**
	 * @brief      Convert an hyperplane representation from a std::vector to a
	 *             std::bitset.
	 *
	 * @details    In the std::vector representation, the vector contain the
	 *             index of all included points in increasing order.
	 *
	 *             In the std::bitset representation, each bit represent a point
	 *             of the hyperplane, included if equal true, excluded if false.
	 *
	 * @param[in]  hyperplane  The hyperplane to convert
	 *
	 * @tparam     NbrPoints   Number of points of the geometry
	 *
	 * @return     The hyperplane converted to std::bitset
	 */
	template<size_t NbrPoints>
	std::bitset<NbrPoints> vectorToBitset(const std::vector<unsigned int>& hyperplane);

	/*------------------------------------------------------------------------*//**
	 * @brief      Utility struct for passing number from
	 *             makeMultiPermutationsGenerator() to
	 *             makeMultiPermutationsGenerator_impl().
	 *
	 * @tparam     N     Number to pass
	 */
	template<size_t N>
	struct nothing {
	};

	/*------------------------------------------------------------------------*//**
	 * @brief      Makes a multi permutations generator with @p Dimension times
	 *             a permutation of @p NbrPointsPerLine elements and one time a
	 *             permutation of @p Dimension elements.
	 *
	 * @details    The multi permutations generator returned is adapted for
	 *             generating all permutations of an hyperplane of a geometry of
	 *             dimension @p Dimension and with @p NbrPointsPerLine points
	 *             per lines.
	 *
	 * @tparam     Dimension         Dimension of the geometry
	 * @tparam     NbrPointsPerLine  Number of points per lines of the geometry
	 *
	 * @return     A multi permutations generator
	 */
	template<size_t Dimension, size_t NbrPointsPerLine>
	auto makeMultiPermutationsGenerator();

	/*------------------------------------------------------------------------*//**
	 * @brief      Implementation of makeMultiPermutationsGenerator(), needed
	 *             to generate an index_repetition
	 */
	template<size_t Dimension, size_t... R>
	MultiPermutationGenerator<R..., Dimension> makeMultiPermutationsGenerator_impl(nothing<Dimension>, index_repetition<R...>);

	/*------------------------------------------------------------------------*//**
	 * @brief      Call @p func on all elements of @p tuple in the increasing
	 *             order.
	 *
	 * @param[in]  func   The function / lambda to call
	 * @param[in]  tuple  The tuple to iterate over
	 *
	 * @tparam     Func   Type of @p func
	 * @tparam     T      Types contained by the tuple
	 */
	template<typename Func, typename... T>
	void iterateOnTuple(Func func, const std::tuple<T...>& tuple);

	/*------------------------------------------------------------------------*//**
	 * @brief      Implementation of iterateOnTuple(), needed to generate an
	 *             std::index_sequence
	 */
	template<typename Func, typename... T, size_t... Is>
	void iterateOnTuple_impl(Func func, const std::tuple<T...>& tuple, std::index_sequence<Is...>);

	/*------------------------------------------------------------------------*//**
	 * @brief      Convert a permutation from a tuple of [@p Dimension times a
	 *             permutation of @p NbrPointsPerLine elements and one time a
	 *             permutation of @p Dimension elements] to a tuple of [an array
	 *             with @p Dimension permutations of @p NbrPointsPerLine
	 *             elements] and [a permutation of @p Dimension elements]
	 *
	 * @param[in]  permutation       The permutation to convert
	 *
	 * @tparam     Dimension         Dimension of the geometry
	 * @tparam     NbrPointsPerLine  Number of points per lines of the geometry
	 * @tparam     Permutation       Type of @p permutation
	 *
	 * @return     The permutation converted
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, typename Permutation>
	std::tuple<std::array<std::array<unsigned int, NbrPointsPerLine>, Dimension>, std::array<unsigned int, Dimension>> convertPermutation(const Permutation& permutation);

	/*------------------------------------------------------------------------*//**
	 * @brief      Apply a permutation to an hyperplane.
	 *
	 * @param[in]  hyperplane        The hyperplane to permute
	 * @param[in]  permutation       The permutation to apply
	 *
	 * @tparam     Dimension         Dimension of the geometry
	 * @tparam     NbrPointsPerLine  Number of points per lines of the geometry
	 * @tparam     Permutation       Type of @p permutation, std::tuple\< @p
	 *                               Dimension times std::array\<unsigned int,
	 *                               @p NbrPointsPerLine\>, std::array\<unsigned
	 *                               int, @p Dimension\>
	 *
	 * @return     The permuted hyperplane
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, typename Permutation>
	std::vector<unsigned int> applyPermutation(const std::vector<unsigned int>& hyperplane, const Permutation& permutation);

	/*------------------------------------------------------------------------*//**
	 * @brief      Apply a coordinates permutation to an hyperplane.
	 *
	 * @param[in]  hyperplane        The hyperplane to permute
	 * @param[in]  permutation       The permutation to apply
	 *
	 * @tparam     Dimension         Dimension of the geometry
	 * @tparam     NbrPointsPerLine  Number of points per lines of the geometry
	 * @tparam     Permutation       Type of @p permutation, std::tuple\< @p
	 *                               Dimension times std::array\<unsigned int,
	 *                               @p NbrPointsPerLine\>\>
	 *
	 * @return     The permuted hyperplane
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, typename Permutation>
	std::vector<unsigned int> applyCoordPermutation(const std::vector<unsigned int>& hyperplane, const Permutation& permutation);

	/*------------------------------------------------------------------------*//**
	 * @brief      Apply a dimensions permutation to an hyperplane.
	 *
	 * @param[in]  hyperplane        The hyperplane to permute
	 * @param[in]  permutation       The permutation to apply
	 *
	 * @tparam     Dimension         Dimension of the geometry
	 * @tparam     NbrPointsPerLine  Number of points per lines of the geometry
	 *
	 * @return     The permuted hyperplane
	 */
	template<size_t Dimension, size_t NbrPointsPerLine>
	std::vector<unsigned int> applyDimensionPermutation(const std::vector<unsigned int>& hyperplane, const std::array<unsigned int, Dimension>& permutation);

	/*------------------------------------------------------------------------*//**
	 * @brief      Calculates the hyperplane stabilisation permutations.
	 *
	 * @details    The hyperplane stabilisation permutations are the
	 *             permutations that doesn't modify the hyperplane, the
	 *             hyperplane stay the same on permutation application.
	 *
	 * @param[in]  hyperplane        The hyperplane
	 *
	 * @tparam     Dimension         Dimension of the geometry
	 * @tparam     NbrPointsPerLine  Number of points per lines of the geometry
	 * @tparam     NbrPoints         Number of points of the geometry
	 *
	 * @return     The hyperplane stabilisation permutations.
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<std::tuple<std::array<std::array<unsigned int, NbrPointsPerLine>, Dimension>, std::array<unsigned int, Dimension>>> computeHyperplaneStabilisationPermutations(std::bitset<NbrPoints> hyperplane);

	/*------------------------------------------------------------------------*//**
	 * @brief      Makes the permutations table, the table of the result of
	 *             applying all possible permutation on each hyperplane of @p
	 *             hyperplanes.
	 *
	 * @details    The element at position (@c n, @c m) is the number of the
	 *             hyperplane from the @p hyperplanes that is the result of
	 *             applying permutation number @c m on hyperplane number @c n
	 *             from the @p hyperplanes.
	 *
	 * @param[in]  hyperplanes       The hyperplanes to generate the table
	 *
	 * @tparam     Dimension         Dimension of the geometry
	 * @tparam     NbrPointsPerLine  Number of points per lines of the geometry
	 * @tparam     NbrPoints         Number of points of the geometry
	 *
	 * @return     The permutations table
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<std::vector<unsigned int>> makePermutationsTable(const std::vector<std::bitset<NbrPoints>>& hyperplanes);
}

// Implementations
namespace segre {

	template<size_t NbrPoints>
	std::vector<unsigned int> bitsetToVector(const std::bitset<NbrPoints>& hyperplane) {
		std::vector<unsigned int> coords;
		coords.reserve(hyperplane.count());
		for(unsigned int i = 0; i < NbrPoints; ++i) {
			if(hyperplane[i]) {
				coords.push_back(i);
			}
		}
		return coords;
	}

	template<size_t NbrPoints>
	std::bitset<NbrPoints> vectorToBitset(const std::vector<unsigned int>& coords) {
		std::bitset<NbrPoints> hyperplane;
		for(unsigned int coord : coords) {
			hyperplane[coord] = true;
		}
		return hyperplane;
	}

	template<size_t Dimension, size_t NbrPointsPerLine>
	auto makeMultiPermutationsGenerator() {
		return makeMultiPermutationsGenerator_impl(nothing<Dimension>{}, make_index_repetition<Dimension, NbrPointsPerLine>());
	}

	template<size_t Dimension, size_t... R>
	MultiPermutationGenerator<R..., Dimension> makeMultiPermutationsGenerator_impl(nothing<Dimension>, index_repetition<R...>) {
		return MultiPermutationGenerator<R..., Dimension>();
	}

	template<typename Func, typename... T>
	void iterateOnTuple(Func func, const std::tuple<T...>& tuple) {
		iterateOnTuple_impl(func, tuple, std::make_index_sequence<sizeof...(T)>());
	}

	template<typename Func, typename... T, size_t... Is>
	void iterateOnTuple_impl(Func func, const std::tuple<T...>& tuple, std::index_sequence<Is...>) {
		(func(std::get<Is>(tuple)), ...);
	}

	template<size_t Dimension, size_t NbrPointsPerLine, typename Permutation>
	std::tuple<std::array<std::array<unsigned int, NbrPointsPerLine>, Dimension>, std::array<unsigned int, Dimension>> convertPermutation(const Permutation& permutation) {
		std::tuple<std::array<std::array<unsigned int, NbrPointsPerLine>, Dimension>, std::array<unsigned int, Dimension>> permutation_value;
		size_t i = 0;
		iterateOnTuple([&](const auto& sub_permutation) {
			if(i < Dimension) {
				// Always true but make the compiler happy as previous if isn't constexpr
				if constexpr (std::is_same<const std::array<unsigned int, NbrPointsPerLine>&, decltype(sub_permutation)>::value) {
					std::get<0>(permutation_value)[i] = sub_permutation;
				}
			}
			else {
				// Always true but make the compiler happy as previous if isn't constexpr
				if constexpr (std::is_same<const std::array<unsigned int, Dimension>&, decltype(sub_permutation)>::value) {
					std::get<1>(permutation_value) = sub_permutation;
				}
			}
			++i;
		}, permutation);
		return permutation_value;
	}

	template<size_t Dimension, size_t NbrPointsPerLine, typename Permutation>
	std::vector<unsigned int> applyPermutation(const std::vector<unsigned int>& hyperplane, const Permutation& permutation) {
		std::vector<unsigned int> permuted_hyperplane;
		permuted_hyperplane.reserve(hyperplane.size());
		for(unsigned int point : hyperplane) {
			std::array<unsigned int, Dimension> permuted_point_coords;
			unsigned int i = 0;
			iterateOnTuple([&](const auto& sub_permutation) {
				if(i < Dimension) {
					// swap coord
					permuted_point_coords[i] = sub_permutation[point % NbrPointsPerLine];
					point /= NbrPointsPerLine;
				}
				else {
					// swap dimension
					for(unsigned int j = 0; j < Dimension; ++j) {
						std::swap(permuted_point_coords[j], permuted_point_coords[sub_permutation[j]]);
					}
				}
				++i;
			}, permutation);
			unsigned int permuted_point = 0;
			for(unsigned int j = 0; j < Dimension; ++j) {
				permuted_point += permuted_point_coords[j] * math::pow(static_cast<unsigned int>(NbrPointsPerLine), j);
			}
			permuted_hyperplane.push_back(permuted_point);
		}
		std::sort(permuted_hyperplane.begin(), permuted_hyperplane.end());
		return permuted_hyperplane;
	}

	template<size_t Dimension, size_t NbrPointsPerLine, typename Permutation>
	std::vector<unsigned int> applyCoordPermutation(const std::vector<unsigned int>& hyperplane, const Permutation& permutation) {
		std::vector<unsigned int> permuted_hyperplane;
		permuted_hyperplane.reserve(hyperplane.size());
		for(unsigned int point : hyperplane) {
			unsigned int permuted_point = 0;
			unsigned int i = 0;
			iterateOnTuple([&](const auto& sub_permutation) {
				permuted_point += sub_permutation[point % NbrPointsPerLine] * math::pow(static_cast<unsigned int>(NbrPointsPerLine), i);
				point /= NbrPointsPerLine;
				++i;
			}, permutation);
			permuted_hyperplane.push_back(permuted_point);
		}
		std::sort(permuted_hyperplane.begin(), permuted_hyperplane.end());
		return permuted_hyperplane;
	}

	template<size_t Dimension, size_t NbrPointsPerLine>
	std::vector<unsigned int> applyDimensionPermutation(const std::vector<unsigned int>& hyperplane, const std::array<unsigned int, Dimension>& permutation) {
		std::vector<unsigned int> permuted_hyperplane;
		permuted_hyperplane.reserve(hyperplane.size());
		for(unsigned int point : hyperplane) {
			unsigned int permuted_point = 0;
			for(unsigned int i = 0; i < Dimension; ++i){
				permuted_point += (point % NbrPointsPerLine) * math::pow(static_cast<unsigned int>(NbrPointsPerLine), permutation[i]);
				point /= NbrPointsPerLine;
			}
			permuted_hyperplane.push_back(permuted_point);
		}
		std::sort(permuted_hyperplane.begin(), permuted_hyperplane.end());
		return permuted_hyperplane;
	}

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints>
	std::vector<std::tuple<std::array<std::array<unsigned int, NbrPointsPerLine>, Dimension>, std::array<unsigned int, Dimension>>> computeHyperplaneStabilisationPermutations(std::bitset<NbrPoints> hyperplane) {
		std::vector<unsigned int> points = bitsetToVector(hyperplane);
		std::vector<std::tuple<std::array<std::array<unsigned int, NbrPointsPerLine>, Dimension>, std::array<unsigned int, Dimension>>> hyperplane_stabilisation_permutations;

		auto multi_permutations_generator = makeMultiPermutationsGenerator<Dimension, NbrPointsPerLine>();
		while(!multi_permutations_generator.isFinished()) {
			const auto current_permutation = multi_permutations_generator.nextPermutation();
			if(points == applyPermutation<Dimension, NbrPointsPerLine>(points, multi_permutations_generator.nextPermutation())) {
				hyperplane_stabilisation_permutations.push_back(convertPermutation<Dimension, NbrPointsPerLine>(current_permutation));
			}
		}

		return hyperplane_stabilisation_permutations;
	}

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints>
	std::vector<std::vector<unsigned int>> makePermutationsTable(const std::vector<std::bitset<NbrPoints>>& hyperplanes) {
		std::vector<std::vector<unsigned int>> permutations_table;
		permutations_table.reserve(hyperplanes.size());

		for(const auto& vPoint : hyperplanes) {
			auto multi_permutations_generator = segre::makeMultiPermutationsGenerator<Dimension, NbrPointsPerLine>();
			std::vector<unsigned int> hyperplane_permutations;
			hyperplane_permutations.reserve(decltype(multi_permutations_generator)::getPermutationsNumber());
			const std::vector<unsigned int> points = segre::bitsetToVector(vPoint);

			while(!multi_permutations_generator.isFinished()) {
				const std::bitset<NbrPoints> hyperplane_permutation = segre::vectorToBitset<NbrPoints>(segre::applyPermutation<Dimension, NbrPointsPerLine>(points, multi_permutations_generator.nextPermutation()));
				const ptrdiff_t pos = std::find(hyperplanes.cbegin(), hyperplanes.cend(), hyperplane_permutation) - hyperplanes.cbegin();
				if(pos > static_cast<ptrdiff_t>(hyperplanes.size())) {
					IMPOSSIBLE;
				}
				hyperplane_permutations.push_back(static_cast<unsigned int>(pos));
			}
			permutations_table.push_back(std::move(hyperplane_permutations));
		}

		return permutations_table;
	}

}


#endif //HYPERPLANEFINDER_HYPERPLANESUTILITY_HPP
