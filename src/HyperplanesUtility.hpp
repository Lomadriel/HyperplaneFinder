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

	template<size_t NbrPoints>
	std::vector<unsigned int> bitsetToVector(const std::bitset<NbrPoints>& hyperplane);

	template<size_t NbrPoints>
	std::bitset<NbrPoints> vectorToBitset(const std::vector<unsigned int>& hyperplane);

	template<size_t N>
	struct nothing {
	};

	template<size_t Dimension, size_t NbrPointsPerLine>
	auto makeMultiPermutationsGenerator();

	template<size_t Dimension, size_t... R>
	MultiPermutationGenerator<R..., Dimension> makeMultiPermutationsGenerator_impl(nothing<Dimension>, index_repetition<R...>);

	template<typename Func, typename... T>
	void iterateOnTuple(Func func, const std::tuple<T...>& tuple);

	template<typename Func, typename... T, size_t... Is>
	void iterateOnTuple_impl(Func func, const std::tuple<T...>& tuple, std::index_sequence<Is...>);

	template<size_t Dimension, size_t NbrPointsPerLine, typename Perm>
	std::tuple<std::array<std::array<unsigned int, NbrPointsPerLine>, Dimension>, std::array<unsigned int, Dimension>> convertPermutation(const Perm& permutation);

	template<size_t Dimension, size_t NbrPointsPerLine, typename Permutation>
	std::vector<unsigned int> applyPermutation(const std::vector<unsigned int>& points, const Permutation& permutation);

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<std::tuple<std::array<std::array<unsigned int, NbrPointsPerLine>, Dimension>, std::array<unsigned int, Dimension>>> computeHyperplaneStabilisationPermutations(std::bitset<NbrPoints> hyperplane);

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

	template<size_t Dimension, size_t NbrPointsPerLine, typename Perm>
	std::tuple<std::array<std::array<unsigned int, NbrPointsPerLine>, Dimension>, std::array<unsigned int, Dimension>> convertPermutation(const Perm& permutation) {
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
	std::vector<unsigned int> applyPermutation(const std::vector<unsigned int>& points, const Permutation& permutation) {
		std::vector<unsigned int> permuted_points;
		permuted_points.reserve(points.size());
		for(unsigned int point : points) {
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
			permuted_points.push_back(permuted_point);
		}
		std::sort(permuted_points.begin(), permuted_points.end());
		return permuted_points;
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
