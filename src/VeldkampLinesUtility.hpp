#ifndef HYPERPLANEFINDER_VELDKAMPLINESUTILITY_HPP
#define HYPERPLANEFINDER_VELDKAMPLINESUTILITY_HPP


#include <vector>

#include "PointGeometry.hpp"
#include "HyperplanesUtility.hpp"

// Declarations
namespace segre{

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::array<unsigned int, NbrPointsPerLine> getPermutation(
	  const std::array<unsigned int, NbrPointsPerLine>& line,
	  const std::vector<std::vector<unsigned int>>& hyp_permutations_table,
	  size_t permutation_number
	);
}

// Implementations
namespace segre {

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints>
	std::array<unsigned int, NbrPointsPerLine> getPermutation(const std::array<unsigned int, NbrPointsPerLine>& line, const std::vector<std::vector<unsigned int>>& hyp_permutations_table, size_t permutation_number) {
		std::array<unsigned int, NbrPointsPerLine> permuted_line;
		for(size_t i = 0; i < NbrPointsPerLine; ++i){
			permuted_line[i] = hyp_permutations_table[line[i]][permutation_number];
		}
		return permuted_line;
	}
}


#endif //HYPERPLANEFINDER_VELDKAMPLINESUTILITY_HPP
