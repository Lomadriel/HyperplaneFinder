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
	 *             entry, apply all permutations to this line, the permuted
	 *             lines generated are equals to lines from the entry, mark
	 *             these lines as checked and add them to the subentry.
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
	 * @brief      Separate a Veldkamp lines table entry by permutations.
	 *             Entries generated contains the Veldcamp lines.
	 *
	 * @details    For informations about separation method see
	 *             separateByPermutations().
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
	 * @return     The Veldkamp lines table entries (with lines) resulting of
	 *             the separation
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>> separateByPermutationsWithLines(
	  const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry,
	  const std::vector<std::vector<unsigned int>>& hyp_permutations_table,
	  unsigned int permutations_number
	);

	/*------------------------------------------------------------------------*//**
	 * @brief      Separate a Veldkamp lines table entry by 2 steps
	 *             permutations.
	 *
	 * @details    2 steps separation works as follow:
	 *
	 *             1. Step 1: Separates lines by coord permutation (see
	 *             separateByPermutationsWithLines() for details)
	 *
	 *             2. Step 2: For each line of each subentry, apply all
	 *             dimensions permutation to the line, if the permuted line is
	 *             in the same subentry, do nothing, if the permuted line is in
	 *             another subentry, join/merge the two subentries
	 *
	 * @param[in]  lines_table_entry                 The Veldkamp lines table
	 *                                               entry (with lines)
	 * @param[in]  hyp_coord_permutations_table      The hyperplanes coordinates
	 *                                               permutations table
	 * @param[in]  hyp_dimension_permutations_table  The hyperplanes dimensions
	 *                                               permutations table
	 *
	 * @tparam     Dimension                         Dimension of the geometry
	 * @tparam     NbrPointsPerLine                  Number of points per lines
	 *                                               of the geometry
	 * @tparam     NbrPoints                         Number of points of the
	 *                                               geometry
	 *
	 * @return     The Veldkamp lines table entries resulting of the separation
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<VeldkampLineTableEntry> separateBy2StepsPermutations(
	  const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry,
	  const std::vector<std::vector<unsigned int>>& hyp_coord_permutations_table,
	  const std::vector<std::vector<unsigned int>>& hyp_dimension_permutations_table
	);

	/*------------------------------------------------------------------------*//**
	 * @brief      Separate a Veldkamp lines table entry by 2 steps
	 *             permutations.
	 *
	 * @details    For informations about separation method see
	 *             separateBy2StepsPermutations().
	 *
	 * @param[in]  lines_table_entry                 The Veldkamp lines table
	 *                                               entry (with lines)
	 * @param[in]  hyp_coord_permutations_table      The hyperplanes coordinates
	 *                                               permutations table
	 * @param[in]  hyp_dimension_permutations_table  The hyperplanes dimensions
	 *                                               permutations table
	 *
	 * @tparam     Dimension                         Dimension of the geometry
	 * @tparam     NbrPointsPerLine                  Number of points per lines
	 *                                               of the geometry
	 * @tparam     NbrPoints                         Number of points of the
	 *                                               geometry
	 *
	 * @return     The Veldkamp lines table entries (with lines) resulting of
	 *             the separation
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>> separateBy2StepsPermutationsWithLines(
	  const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry,
	  const std::vector<std::vector<unsigned int>>& hyp_coord_permutations_table,
	  const std::vector<std::vector<unsigned int>>& hyp_dimension_permutations_table
	);

	/*------------------------------------------------------------------------*//**
	 * @brief      Separate entries of a Veldkamp lines table by permutations
	 *
	 * @details    See separateByPermutations() on a Veldkamp line table entry
	 *             for separation method details.
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

	/*------------------------------------------------------------------------*//**
	 * @brief      Separate entries of a Veldkamp lines table by 2 steps
	 *             permutations.
	 *
	 * @param[in]  lin_table_with_lines              The Veldkamp lines table
	 *                                               (with lines)
	 * @param[in]  hyp_coord_permutations_table      The hyperplanes coordinates
	 *                                               permutations table
	 * @param[in]  hyp_dimension_permutations_table  The hyperplanes dimensions
	 *                                               permutations table
	 *
	 * @tparam     Dimension                         Dimension of the geometry
	 * @tparam     NbrPointsPerLine                  Number of points per lines
	 *                                               of the geometry
	 * @tparam     NbrPoints                         Number of points of the
	 *                                               geometry
	 *
	 * @return     The Veldkamp lines table resulting of the separation
	 */
	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<VeldkampLineTableEntry> separateBy2StepsPermutations(
	  const std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>>& lin_table_with_lines,
	  const std::vector<std::vector<unsigned int>>& hyp_coord_permutations_table,
	  const std::vector<std::vector<unsigned int>>& hyp_dimension_permutations_table
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

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints>
	std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>> separateByPermutationsWithLines(const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry, const std::vector<std::vector<unsigned int>>& hyp_permutations_table, unsigned int permutations_number) {
		std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>> output_table;
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
			VeldkampLineTableEntryWithLines<NbrPointsPerLine> table_entry(lines_table_entry.entry);
			table_entry.entry.count = 1;
			table_entry.lines.push_back(line);
			for(size_t i = 0; i < permutations_number; ++i){
				std::array<unsigned int, NbrPointsPerLine> permuted_line = applyPermutation<Dimension, NbrPointsPerLine>(line, hyp_permutations_table, i);
				std::sort(permuted_line.begin(), permuted_line.end());
				const size_t pos2 = static_cast<const size_t>(std::find(sorted_lines.cbegin(), sorted_lines.cend(), permuted_line) - sorted_lines.cbegin());
				if(pos2 < sorted_lines.size()) {
					if(!flags[pos2]){
						flags[pos2] = true;
						++table_entry.entry.count;
						table_entry.lines.push_back(lines_table_entry.lines[pos2]);
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
	std::vector<VeldkampLineTableEntry> separateBy2StepsPermutations(
	  const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry,
	  const std::vector<std::vector<unsigned int>>& hyp_coord_permutations_table,
	  const std::vector<std::vector<unsigned int>>& hyp_dimension_permutations_table
	){
		// Separate by coord permutations
		std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>> coord_permuted_entries =
		  separateByPermutationsWithLines<Dimension, NbrPointsPerLine>(
		    lines_table_entry,
		      hyp_coord_permutations_table,
		      decltype(makeCoordPermutationsGenerator<Dimension, NbrPointsPerLine>())::getPermutationsNumber()
		  );

		// Prepare table of sorted entries lines
		std::vector<std::vector<std::array<unsigned int, NbrPointsPerLine>>> entries_sorted_lines;
		entries_sorted_lines.resize(coord_permuted_entries.size());
		for(unsigned int i_entry = 0; i_entry < coord_permuted_entries.size(); ++i_entry){
			entries_sorted_lines[i_entry].reserve(coord_permuted_entries[i_entry].lines.size());
			for(unsigned int i_line = 0; i_line < coord_permuted_entries[i_entry].lines.size(); ++i_line){
				std::array<unsigned int, NbrPointsPerLine> sorted_line = coord_permuted_entries[i_entry].lines[i_line];
				std::sort(sorted_line.begin(), sorted_line.end());
				entries_sorted_lines[i_entry].push_back(std::move(sorted_line));
			}
		}

		// Prepare processed entries joined flags and output table position
		std::vector<bool> entry_is_joined;
		std::vector<unsigned int> entry_output_table_pos;
		entry_is_joined.resize(coord_permuted_entries.size(), false);
		entry_output_table_pos.resize(coord_permuted_entries.size(), 0);

		// Process entries: join with dimensions permutations
		std::vector<VeldkampLineTableEntry> output_table;
		for(unsigned int i_entry = 0; i_entry < coord_permuted_entries.size(); ++i_entry){
			const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& entry = coord_permuted_entries[i_entry];
			if(entry_is_joined[i_entry]){
				output_table[entry_output_table_pos[i_entry]].count += coord_permuted_entries[i_entry].entry.count;
			}
			else{
				entry_is_joined[i_entry] = true;
				entry_output_table_pos[i_entry] = static_cast<unsigned int>(output_table.size());
				output_table.push_back(coord_permuted_entries[i_entry].entry);
			}

			for(const std::array<unsigned int, NbrPointsPerLine> line : entry.lines){
				for(size_t i_permutation = 0; i_permutation < PermutationGenerator<Dimension>::getPermutationsNumber(); ++i_permutation){
					std::array<unsigned int, NbrPointsPerLine> permuted_line = applyPermutation<Dimension, NbrPointsPerLine>(line, hyp_dimension_permutations_table, i_permutation);
					std::sort(permuted_line.begin(), permuted_line.end());

					bool found = false;
					for(unsigned int i_search_entry = 0; i_search_entry < coord_permuted_entries.size(); ++i_search_entry){
						const size_t found_pos = static_cast<const size_t>(std::find(entries_sorted_lines[i_search_entry].begin(), entries_sorted_lines[i_search_entry].end(), permuted_line) - entries_sorted_lines[i_search_entry].begin());
						if(found_pos < entries_sorted_lines[i_search_entry].size()){
							if(entry_is_joined[i_search_entry] && entry_output_table_pos[i_search_entry] != entry_output_table_pos[i_entry]){
								IMPOSSIBLE;
							}
							// Join i_search_entry with current entry (i_entry)
							entry_is_joined[i_search_entry] = true;
							entry_output_table_pos[i_search_entry] = entry_output_table_pos[i_entry];
							found = true;
							break;
						}
					}
					if(!found){
						IMPOSSIBLE;
					}
				}
			}
		}

		return output_table;
	}

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>> separateBy2StepsPermutationsWithLines(
	  const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry,
	  const std::vector<std::vector<unsigned int>>& hyp_coord_permutations_table,
	  const std::vector<std::vector<unsigned int>>& hyp_dimension_permutations_table
	){
		// Separate by coord permutations
		std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>> coord_permuted_entries =
		separateByPermutationsWithLines<Dimension, NbrPointsPerLine>(
		  lines_table_entry,
			hyp_coord_permutations_table,
			decltype(makeCoordPermutationsGenerator<Dimension, NbrPointsPerLine>())::getPermutationsNumber()
		);

		// Prepare table of sorted entries lines
		std::vector<std::vector<std::array<unsigned int, NbrPointsPerLine>>> entries_sorted_lines;
		entries_sorted_lines.resize(coord_permuted_entries.size());
		for(unsigned int i_entry = 0; i_entry < coord_permuted_entries.size(); ++i_entry){
			entries_sorted_lines[i_entry].reserve(coord_permuted_entries[i_entry].lines.size());
			for(unsigned int i_line = 0; i_line < coord_permuted_entries[i_entry].lines.size(); ++i_line){
				std::array<unsigned int, NbrPointsPerLine> sorted_line = coord_permuted_entries[i_entry].lines[i_line];
				std::sort(sorted_line.begin(), sorted_line.end());
				entries_sorted_lines[i_entry].push_back(std::move(sorted_line));
			}
		}

		// Prepare processed entries joined flags and output table position
		std::vector<bool> entry_is_joined;
		std::vector<unsigned int> entry_output_table_pos;
		entry_is_joined.resize(coord_permuted_entries.size(), false);
		entry_output_table_pos.resize(coord_permuted_entries.size(), 0);

		// Process entries: join with dimensions permutations
		std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>> output_table;
		for(unsigned int i_entry = 0; i_entry < coord_permuted_entries.size(); ++i_entry){
			const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& entry = coord_permuted_entries[i_entry];
			if(entry_is_joined[i_entry]){
				output_table[entry_output_table_pos[i_entry]].entry.count += coord_permuted_entries[i_entry].entry.count;
				std::move(coord_permuted_entries[i_entry].lines.begin(), coord_permuted_entries[i_entry].lines.end(), output_table[entry_output_table_pos[i_entry]].lines.end());
			}
			else{
				entry_is_joined[i_entry] = true;
				entry_output_table_pos[i_entry] = static_cast<unsigned int>(output_table.size());
				output_table.push_back(coord_permuted_entries[i_entry]);
			}

			for(const std::array<unsigned int, NbrPointsPerLine> line : entry.lines){
				for(size_t i_permutation = 0; i_permutation < PermutationGenerator<Dimension>::getPermutationsNumber(); ++i_permutation){
					std::array<unsigned int, NbrPointsPerLine> permuted_line = applyPermutation<Dimension, NbrPointsPerLine>(line, hyp_dimension_permutations_table, i_permutation);
					std::sort(permuted_line.begin(), permuted_line.end());

					bool found = false;
					for(unsigned int i_search_entry = 0; i_search_entry < coord_permuted_entries.size(); ++i_search_entry){
						const size_t found_pos = static_cast<const size_t>(std::find(entries_sorted_lines[i_search_entry].begin(), entries_sorted_lines[i_search_entry].end(), permuted_line) - entries_sorted_lines[i_search_entry].begin());
						if(found_pos < entries_sorted_lines[i_search_entry].size()){
							if(entry_is_joined[i_search_entry] && entry_output_table_pos[i_search_entry] != entry_output_table_pos[i_entry]){
								IMPOSSIBLE;
							}
							// Join i_search_entry with current entry (i_entry)
							entry_is_joined[i_search_entry] = true;
							entry_output_table_pos[i_search_entry] = entry_output_table_pos[i_entry];
							found = true;
							break;
						}
					}
					if(!found){
						IMPOSSIBLE;
					}
				}
			}
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

	template<size_t Dimension, size_t NbrPointsPerLine, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension)>
	std::vector<VeldkampLineTableEntry> separateBy2StepsPermutations(
	  const std::vector<VeldkampLineTableEntryWithLines<NbrPointsPerLine>>& lin_table_with_lines,
	  const std::vector<std::vector<unsigned int>>& hyp_coord_permutations_table,
	  const std::vector<std::vector<unsigned int>>& hyp_dimension_permutations_table
	){
		std::vector<VeldkampLineTableEntry> output_table;
		for(const VeldkampLineTableEntryWithLines<NbrPointsPerLine>& lines_table_entry : lin_table_with_lines){
			for(VeldkampLineTableEntry& entry : separateBy2StepsPermutations<Dimension, NbrPointsPerLine>(lines_table_entry, hyp_coord_permutations_table, hyp_dimension_permutations_table)){
				output_table.push_back(std::move(entry));
			}
		}
		return output_table;
	}
}


#endif //HYPERPLANEFINDER_VELDKAMPLINESUTILITY_HPP
