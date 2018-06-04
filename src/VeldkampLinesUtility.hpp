#ifndef HYPERPLANEFINDER_VELDKAMPLINESUTILITY_HPP
#define HYPERPLANEFINDER_VELDKAMPLINESUTILITY_HPP


#include <vector>

#include "PointGeometry.hpp"
#include "HyperplanesUtility.hpp"

// Declarations
namespace segre {

	/*------------------------------------------------------------------------*//**
	 * @brief      Apply the permutation number @p permutation_number to the
	 *             Veldkamp line @p line using the hyperplanes permutation table
	 *             @p hyp_permutations_table.
	 *
	 * @param[in]  line                    The line to permute
	 * @param[in]  hyp_permutations_table  The hyperplanes permutations table
	 *                                     (see makePermutationsTable())
	 * @param[in]  permutation_number      The permutation number
	 *
	 * @tparam     Dimension               Dimension of the geometry
	 * @tparam     NbrPointsPerLine        Number of points per lines of the
	 *                                     geometry
	 * @tparam     NbrPoints               Number of points of the geometry
	 *
	 * @return     The permuted Veldkamp line.
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::array<unsigned int, NbrPointsPerLine> applyPermutation(
	  const std::array<unsigned int, NbrPointsPerLine>& line,
	  const std::vector<std::vector<unsigned int>>& hyp_permutations_table,
	  size_t permutation_number
	);

	/*------------------------------------------------------------------------*//**
	 * @brief      Separate a Veldkamp lines table entry by permutations.
	 *
	 * @details    Separation works as follow:
	 *
	 *             1. All lines are marked as unchecked
	 *
	 *             2. Create a subentry, take the first unchecked line from the
	 *             entry, apply all possible permutations to this line, the
	 *             permuted lines generated are equals to lines from the entry,
	 *             mark these lines as checked and add them to the subentry.
	 *
	 *             3. While there is unchecked lines in the entry, do 2
	 *
	 *             4. Return the subentries
	 *
	 * @param[in]  lines_table_entry       The Veldkamp lines table entry (with
	 *                                     lines)
	 * @param[in]  hyp_permutations_table  The hyperplanes permutations table
	 * @param[in]  permutations_number     The permutations number
	 *
	 * @tparam     Dimension               Dimension of the geometry
	 * @tparam     NbrPointsPerLine        Number of points per lines of the
	 *                                     geometry
	 * @tparam     NbrPoints               Number of points of the geometry
	 *
	 * @return     The Veldkamp lines table entries resulting of the separation
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<VeldkampLineTableEntry> separateByPermutations(
	  const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry,
	  const std::vector<std::vector<unsigned int>>& hyp_permutations_table,
	  unsigned int permutations_number
	);

	/*------------------------------------------------------------------------*//**
	 * @brief      Separate entries of a Veldkamp lines table by permutations
	 *
	 * @details    See separateByPermutations() on a Veldkamp line table entry
	 *             for separation method details
	 *
	 * @param[in]  lin_table_with_lines    The Veldkamp lines table (with lines)
	 * @param[in]  hyp_permutations_table  The hyperplanes permutations table
	 *                                     (see makePermutationsTable())
	 *
	 * @tparam     Dimension               Dimension of the geometry
	 * @tparam     NbrPointsPerLine        Number of points per lines of the
	 *                                     geometry
	 * @tparam     NbrPoints               Number of points of the geometry
	 *
	 * @return     The Veldkamp lines table resulting of the separation
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<VeldkampLineTableEntry> separateByPermutations(
	  const std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>>& lin_table_with_lines,
	  const std::vector<std::vector<unsigned int>>& hyp_permutations_table
	);
}

// Implementations
namespace segre {

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints>
	std::array<unsigned int, NbrPointsPerLine> applyPermutation(const std::array<unsigned int, NbrPointsPerLine>& line, const std::vector<std::vector<unsigned int>>& hyp_permutations_table, size_t permutation_number) {
		std::array<unsigned int, NbrPointsPerLine> permuted_line;
		for(size_t i = 0; i < NbrPointsPerLine; ++i){
			permuted_line[i] = hyp_permutations_table[line[i]][permutation_number];
		}
		return permuted_line;
	}

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints>
	std::vector<VeldkampLineTableEntry> separateByPermutations(const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry, const std::vector<std::vector<unsigned int>>& hyp_permutations_table, unsigned int permutations_number) {
		std::vector<VeldkampLineTableEntry> output_table;
		std::vector<std::array<unsigned int, NbrPointsPerLine>> sorted_lines(lines_table_entry.lines);
		std::for_each(sorted_lines.begin(), sorted_lines.end(), [](std::array<unsigned int, NbrPointsPerLine>& line){std::sort(line.begin(), line.end());});

		std::vector<bool> flags;
		flags.reserve(lines_table_entry.lines.size());
		for(size_t i = 0; i < lines_table_entry.lines.size(); ++i){
			flags.push_back(false);
		}
		while(true){
			const size_t pos = static_cast<const size_t>(std::find(flags.cbegin(), flags.cend(), false) - flags.cbegin());
			if(pos >= flags.size()) {
				break;
			}
			const std::array<unsigned int, NbrPointsPerLine>& line = lines_table_entry.lines[pos];
			flags[pos] = true;
			VeldkampLineTableEntry table_entry = lines_table_entry.entry;
			table_entry.count = 1;
			for(size_t i = 0; i < permutations_number; ++i){
				std::array<unsigned int, NbrPointsPerLine> permuted_line = applyPermutation<Dimension, NbrPointsPerLine>(line, hyp_permutations_table, i);
				std::sort(permuted_line.begin(), permuted_line.end());
				const size_t pos2 = static_cast<const size_t>(std::find(sorted_lines.cbegin(), sorted_lines.cend(), permuted_line) - sorted_lines.cbegin());
				if(pos2 < sorted_lines.size()) {
					if(!flags[pos2]){
						flags[pos2] = true;
						++table_entry.count;
					}
				}
				else{
					IMPOSSIBLE;
				}
			}
			output_table.push_back(std::move(table_entry));
		}

		return output_table;
	}

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<VeldkampLineTableEntry> separateByPermutations(
	  const std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>>& lin_table_with_lines,
	  const std::vector<std::vector<unsigned int>>& hyp_permutations_table
	){
		std::vector<VeldkampLineTableEntry> output_table;
		for(const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry : lin_table_with_lines){
			for(VeldkampLineTableEntry& entry : separateByPermutations<Dimension, NbrPointsPerLine>(lines_table_entry, hyp_permutations_table, decltype(makeMultiPermutationsGenerator<Dimension, NbrPointsPerLine>())::getPermutationsNumber())){
				output_table.push_back(std::move(entry));
			}
		}
		return output_table;
	}
}


#endif //HYPERPLANEFINDER_VELDKAMPLINESUTILITY_HPP
