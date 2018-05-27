#ifndef HYPERPLANEFINDER_POINTGEOMETRY_HPP
#define HYPERPLANEFINDER_POINTGEOMETRY_HPP

#include <cstddef>
#include <array>
#include <bitset>
#include <vector>
#include <functional>
#include <algorithm>
#include <utility>
#include <map>
#include <iostream>
#include <set>

#include "CombinationGenerator.hpp"
#include "math.hpp"
#include "impossible.hpp"
#include "HyperplaneTableEntry.hpp"
#include "VeldkampLineTableEntry.hpp"

// Fixme : Replace std::bitset by a custom bitset
namespace std { // NOLINT
	template <size_t N>
	bool operator<(const std::bitset<N>& x, const std::bitset<N>& y) {
		for (size_t i = N - 1; i--;) {
			if (x[i] ^ y[i]) {
				return y[i];
			}
		}

		return false;
	}
}

namespace segre {
	template <std::size_t NbrPointsPerLine>
	struct VeldkampLines {
		explicit VeldkampLines(std::vector<std::array<unsigned int, NbrPointsPerLine>>&& exceptional_lines,
		                       std::vector<std::array<unsigned int, NbrPointsPerLine>>&& projectives_lines) noexcept;

		std::vector<std::array<unsigned int, NbrPointsPerLine>> exceptional;
		std::vector<std::array<unsigned int, NbrPointsPerLine>> projectives;
	};

	template <size_t N1, size_t N2>
	inline std::bitset<N1> copyBitset(const std::bitset<N2>& bs2) {
		std::bitset<N1> bs1;
		for (size_t i = 0; i < N2; i++) {
			bs1[i] = bs2[i];
		}

		return bs1;
	}

	const std::array<std::array<unsigned int, 2>, 4> TENSOR_2D = {{ {{1, 0}}, {{0, 1}}, {{1, 1}}, {{1, 2}} }};

	template <size_t Dimension, size_t NbrPointsPerLine,
			size_t NbrLines, size_t NbrPoints = math::pow(NbrPointsPerLine, Dimension),
			size_t TensorSize = math::pow(2UL, Dimension)>
	class PointGeometry {
	public:
		explicit PointGeometry(std::array<std::bitset<NbrPoints>, NbrLines>&& lines) noexcept
				: m_geometryLines(std::move(lines))
				, m_geometryPoints(TENSOR_2D)
				, m_subGeometriesMasks() {
			computeMasks();
		}

		explicit PointGeometry(std::array<std::bitset<NbrPoints>, NbrLines>&& lines,
		                       std::array<std::array<unsigned int, TensorSize>, NbrPoints>&& tensors) noexcept
				: m_geometryLines(std::move(lines))
				, m_geometryPoints(std::move(tensors))
		        , m_subGeometriesMasks() {
			computeMasks();
		}

		/**
		 * @details Computes the hyperplanes of the geometry by checking every combination of k points
		 * 	where k is in [0, NbrPoints - 1]
		 *
		 * @return a vector of hyperplanes
		 */
		std::vector<std::bitset<NbrPoints>> findHyperplanesByBruteforce() const noexcept {
			std::vector<std::bitset<NbrPoints>> hyperplanes;

			for (unsigned int j = 2; j < NbrPoints; ++j) {
				permutations<std::uint64_t>(NbrPoints, j, 0, hyperplanes);
			}

			return std::move(hyperplanes);
		}

		/**
		 * Checks if the given combination is an hyperplane.
		 *
		 * @param potentialHyperplane
		 * @return true if potentialHyperplane is an hyperplane, false otherwise.
		 */
		bool isHyperplane(const std::bitset<NbrPoints>& potentialHyperplane) const noexcept {
			for (const std::bitset<NbrPoints>& line : m_geometryLines) {
				std::bitset<NbrPoints> intersection = line & potentialHyperplane;
				std::size_t intersectionSize = intersection.count();

				if (intersectionSize == 0) {
					return false;
				}

				if (intersectionSize > 1) {
					if ((line & potentialHyperplane) != line) { // Test if line is included in potentialHyperplane.
						return false;
					}
				}
			}

			return true;
		}

		/**
		 * Computes the veldkamp lines of the geometry using the given hyperplanes.
		 *
		 * @param veldkampPoints the hyperplanes of the geometry.
		 * @return A struct containing the projective lines and the supposed exceptional lines.
		 */
		VeldkampLines<NbrPointsPerLine>
		computeVeldkampLines(const std::vector<std::bitset<NbrPoints>>& veldkampPoints) const noexcept {
			std::vector<std::array<unsigned int, NbrPointsPerLine>> supposedExceptional;
			std::vector<std::array<unsigned int, NbrPointsPerLine>> projectiveLines;

			CombinationGenerator gen;
			gen.initialize(static_cast<unsigned int>(veldkampPoints.size()), 2);

			while (!gen.isFinished()) {
				const std::vector<unsigned int>& currentCombination = gen.nextCombination();

				std::vector<unsigned int> sameCore;
				const std::bitset<NbrPoints>& h1 = veldkampPoints[currentCombination[0]];
				const std::bitset<NbrPoints>& h2 = veldkampPoints[currentCombination[1]];

				const std::bitset<NbrPoints> intersection12 = h1 & h2;

				for (size_t i = 0, n = veldkampPoints.size(); i < n; ++i) {
					const std::bitset<NbrPoints> intersection1i = h1 & veldkampPoints[i];
					const std::bitset<NbrPoints> intersection2i = h2 & veldkampPoints[i];

					if (intersection12 == intersection1i
					    && intersection12 == intersection2i
					    && intersection1i == intersection2i) {
						sameCore.emplace_back(i);
					}
				}

				CombinationGenerator gen2;
				gen2.initialize(static_cast<unsigned int>(sameCore.size()), 2);

				while (!gen2.isFinished()) {
					const std::vector<unsigned int>& currentCombination2 = gen2.nextCombination();

					if (sameCore[currentCombination2[0]] > currentCombination[1]) {
						const std::bitset<NbrPoints>& ha = veldkampPoints[sameCore[currentCombination2[0]]];
						const std::bitset<NbrPoints>& hb = veldkampPoints[sameCore[currentCombination2[1]]];

						const std::bitset<NbrPoints> intersectionAB = ha & hb;

						if (intersection12 == intersectionAB) {
							if (sameCore.size() == 2) {
								projectiveLines.emplace_back(
										std::array<unsigned int, NbrPointsPerLine>({currentCombination[0],
										                                            currentCombination[1],
										                                            sameCore[currentCombination2[0]],
										                                            sameCore[currentCombination2[1]]}));
							} else {
								supposedExceptional.emplace_back(
										std::array<unsigned int, NbrPointsPerLine>({currentCombination[0],
										                                            currentCombination[1],
										                                            sameCore[currentCombination2[0]],
										                                            sameCore[currentCombination2[1]]}));
							}
						}
					}
				}
			}

			return VeldkampLines<NbrPointsPerLine>{std::move(supposedExceptional), std::move(projectiveLines)};
		}

		/**
		 * Excludes the projectives lines from the list of exceptional lines.
		 * @param vLines a struct containing the projective and exceptional lines.
		 * @param vPoints the hyperplanes of the current geometry
		 * @param nextGeometry the next geometry used to build the matrix associated to the hyperplanes.
		 */
		void distinguishVeldkampLines(VeldkampLines<NbrPointsPerLine>& vLines,
		                              const std::vector<std::bitset<NbrPoints>>& vPoints,
		                              const PointGeometry<Dimension + 1,
		                                                  NbrPointsPerLine,
		                                                  math::pow(NbrPointsPerLine, Dimension) * (1 + Dimension)>& nextGeometry) const {
			std::vector<size_t> toRemove;
			constexpr size_t NewNbrPoints = math::pow(NbrPointsPerLine, Dimension + 1);

			// Checks the rank of the matrix associated to each hyperplane.
			// If the rank of the matrix is lesser than pow(2, Dimension + 1) then the line isn't exceptional.
			for (size_t index = 0; index < vLines.exceptional.size(); ++index) {
				std::bitset<NewNbrPoints> hyperplane;
				for (size_t i = 0; i < vLines.exceptional[index].size(); ++i) {
					hyperplane |= copyBitset<NewNbrPoints>(vPoints[vLines.exceptional[index][i]]) <<= (i * NbrPoints);
				}

				// Checks if the matrix associated to the hyperplane live in the projective space.
				if (getRank(nextGeometry.buildMatrix(hyperplane)) < math::pow(2UL, Dimension + 1)) {
					toRemove.push_back(index);
				}
			}

			vLines.projectives.reserve(vLines.projectives.size() + toRemove.size());

			// Remove the projective lines from the list of exceptional lines.
			while (!toRemove.empty()) {
				vLines.projectives.push_back(vLines.exceptional[toRemove.back()]);
				vLines.exceptional.erase(std::next(vLines.exceptional.begin(), static_cast<unsigned int>(toRemove.back())));

				toRemove.pop_back();
			}
		}

		size_t getRank(std::vector<std::array<unsigned int, math::pow(2UL, Dimension + 1)>>&& matrix) const {
			size_t rank = math::pow(2UL, Dimension + 1);

			for (size_t i = 0; i < math::pow(2UL, Dimension + 1); ++i) {
				if (matrix[i][i] == 0) {
					size_t k = i;
					while (k < matrix.size() && matrix[k][i] == 0) {
						++k;
					}

					if (k != matrix.size()) {
						std::swap(matrix[i], matrix[k]);
					} else {
						rank = i;
						break;
					}
				}

				if (matrix[i][i] == 2) {
					for (size_t j = 0; j < math::pow(2UL, Dimension + 1); ++j) {
						matrix[i][j] = (matrix[i][j] * 2) % 3;
					}
				}

				for (size_t j = i + 1; j < matrix.size(); ++j) {
					unsigned int aij = matrix[j][i];
					for (size_t k = i; k < math::pow(2UL, Dimension + 1); ++k) {
						matrix[j][k] = (6 + matrix[j][k] - aij * matrix[i][k]) % 3;
					}
				}
			}

			if (rank != math::pow(2UL, Dimension + 1)) {
				for (size_t i = rank + 1; i < math::pow(2UL, Dimension + 1); ++i) {
					for (size_t k = 0; k < matrix.size(); ++k) {
						if (matrix[k][i] != 0) {
							++rank;
							break;
						}
					}
				}
			}

			return rank;
		}

		decltype(auto) computeHyperplanesFromVeldkampLines(const std::vector<std::bitset<NbrPoints>>& veldkampPoints,
		                                                   const std::vector<std::array<unsigned int, NbrPointsPerLine>>& pVLines) {
			constexpr size_t NewNbrPoints = math::pow(NbrPointsPerLine, Dimension + 1);

			std::vector<std::bitset<NewNbrPoints>> hyperplanes;

			// Compute the hyperplane of the next geometry using the veldkamp lines of the current geometry.
			for (size_t i = 0; i < pVLines.size(); ++i) {
				std::array<std::bitset<NbrPoints>, NbrPointsPerLine> hypers = getHyperplanesOfTheVeldkampLine(veldkampPoints, pVLines[i]);
				std::sort(hypers.begin(), hypers.end());

				do {
					std::bitset<NewNbrPoints> hyperplane;

					for (size_t j = 0; j < hypers.size(); ++j) {
						hyperplane |= copyBitset<NewNbrPoints>(hypers[j]) <<= (j * NbrPoints);
					}

					hyperplanes.push_back(std::move(hyperplane));
				} while (std::next_permutation(hypers.begin(), hypers.end()));
			}

			std::bitset<NewNbrPoints> fullLayout = copyBitset<NewNbrPoints>(std::bitset<NbrPoints>().flip());

			// Compute the missing hyperplanes by using 3 times the same hyperplane and the current full geometry.
			for (size_t i = 0; i < veldkampPoints.size(); ++i) {
				std::bitset<NewNbrPoints> hyperplane = copyBitset<NewNbrPoints>(veldkampPoints[i]);
				hyperplane |= copyBitset<NewNbrPoints>(veldkampPoints[i]) <<= NbrPoints;
				hyperplane |= copyBitset<NewNbrPoints>(veldkampPoints[i]) <<= 2 * NbrPoints;
				hyperplane |= fullLayout << 3 * NbrPoints;

				hyperplanes.push_back(std::move(hyperplane));

				hyperplane = copyBitset<NewNbrPoints>(veldkampPoints[i]);
				hyperplane |= copyBitset<NewNbrPoints>(veldkampPoints[i]) << NbrPoints;
				hyperplane |= fullLayout << 2 * NbrPoints;
				hyperplane |= copyBitset<NewNbrPoints>(veldkampPoints[i]) <<= 3 * NbrPoints;

				hyperplanes.push_back(std::move(hyperplane));

				hyperplane = copyBitset<NewNbrPoints>(veldkampPoints[i]);
				hyperplane |= fullLayout << NbrPoints;
				hyperplane |= copyBitset<NewNbrPoints>(veldkampPoints[i]) << 2 * NbrPoints;
				hyperplane |= copyBitset<NewNbrPoints>(veldkampPoints[i]) <<= 3 * NbrPoints;

				hyperplanes.push_back(std::move(hyperplane));

				hyperplane = fullLayout;
				hyperplane |= copyBitset<NewNbrPoints>(veldkampPoints[i]) <<= NbrPoints;
				hyperplane |= copyBitset<NewNbrPoints>(veldkampPoints[i]) <<= 2 * NbrPoints;
				hyperplane |= copyBitset<NewNbrPoints>(veldkampPoints[i]) <<= 3 * NbrPoints;

				hyperplanes.push_back(std::move(hyperplane));
			}

			return hyperplanes;
		}

		/**
		 * Returns the hyperplanes of the given veldkamp line.
		 * @param veldkampPoints list of all the veldkamp points of the current geometry.
		 * @param veldkampLine a veldkamp line.
		 * @return the hyperplanes of the given veldkamp line.
		 */
		static std::array<std::bitset<NbrPoints>, NbrPointsPerLine>
		getHyperplanesOfTheVeldkampLine(const std::vector<std::bitset<NbrPoints>>& veldkampPoints,
		                                const std::array<unsigned int, NbrPointsPerLine>& veldkampLine) {
			std::array<std::bitset<NbrPoints>, NbrPointsPerLine> hyperplanes;

			for (size_t i = 0; i < veldkampLine.size(); ++i) {
				hyperplanes[i] = veldkampPoints[veldkampLine[i]];
			}

			return hyperplanes;
		}

		decltype(auto) computeCartesianProduct() const noexcept {
			constexpr size_t NewNbrPoints = math::pow(NbrPointsPerLine, Dimension + 1);
			constexpr size_t NewNbrLines = math::pow(NbrPointsPerLine, Dimension) * (1 + Dimension);

			std::array<std::bitset<NewNbrPoints>, NewNbrLines> result;

			// Duplicates the current geometry to generate each layer of the cartesian product.
			std::generate(result.begin(), result.end(), [this, i = 0UL, j = 0UL]() mutable -> decltype(auto) {
				auto bitset = copyBitset<math::pow(NbrPointsPerLine, Dimension + 1)>(m_geometryLines[j]) <<=
						              static_cast<long unsigned int>(NbrPoints * i);
				++j;
				if (j % NbrLines == 0) {
					j = 0;
					++i;
				}

				return bitset;
			});

			// Computes the missing lines linking each layer.
			for (size_t i = 0; i < NbrPoints; ++i) {
				std::bitset<NewNbrPoints> line;
				for (size_t j = 0; j < NbrPointsPerLine; ++j) {
					line |= std::bitset<NewNbrPoints>(1) <<= (NbrPoints * j + i);
				}

				result[NbrLines * NbrPointsPerLine + i] = line;
			}

			return result;
		}

		decltype(auto) buildTensorPoints() const noexcept {
			constexpr size_t NewNbrPoints = math::pow(NbrPointsPerLine, Dimension + 1);

			std::array<std::array<unsigned int, TensorSize * 2>, NewNbrPoints> pts;

			for (size_t i = 0; i < NbrPointsPerLine; ++i) {
				for (size_t j = 0; j < NbrPoints; ++j) {
					for (unsigned int k = 0; k < TENSOR_2D[0].size(); ++k) {
						for (unsigned int l = 0; l < TensorSize; ++l) {
							// The mod operator is here because the coefficient in the associated space are {0, 1, 2}
							pts[i * NbrPoints + j][k * TensorSize + l] = TENSOR_2D[i][k] * m_geometryPoints[j][l] % 3;
						}
					}
				}
			}

			return pts;
		}

		std::vector<std::array<unsigned int, TensorSize>> buildMatrix(const std::bitset<NbrPoints>& veldkampPoint) const noexcept {
			std::vector<std::array<unsigned int, TensorSize>> matrix;

			for (size_t i = 0; i < NbrPoints; ++i) {
				if (veldkampPoint[i]) {
					matrix.push_back(m_geometryPoints[i]);
				}
			}

			return matrix;
		}

		template <bool OrderOfPoints>
		HyperplaneTableEntry getHyperplaneTableEntry(const std::bitset<NbrPoints>& hyperplane) const noexcept {
			HyperplaneTableEntry entry;
			entry.nbrPoints = static_cast<unsigned int>(hyperplane.count());

			std::vector<std::bitset<NbrPoints>> includedLines;
			for (const std::bitset<NbrPoints>& line : m_geometryLines) {
				if ((line & hyperplane) == line) {
					includedLines.push_back(line);
				}
			}

			entry.nbrLines = static_cast<unsigned int>(includedLines.size());

			if constexpr (OrderOfPoints) {
				if (entry.nbrLines == 0) {
					entry.pointsOfOrder[0] = entry.nbrPoints;
				} else {
					unsigned int pointOfOrder0 = entry.nbrPoints;
					for (size_t i = 0; i < NbrPoints; ++i) {
						unsigned int count = 0;
						if (hyperplane[i]) {
							for (size_t j = 0; j < includedLines.size(); ++j) {
								if (includedLines[j][i] != 0) {
									++count;
								}
							}

							if (count != 0) {
								++(entry.pointsOfOrder[count]);
								--pointOfOrder0;
							}
						}
					}

					if (pointOfOrder0 != 0) {
						entry.pointsOfOrder[0] = pointOfOrder0;
					}
				}
			}

			return entry;
		}

		template <bool OrderOfPoints>
		HyperplaneTableEntry getHyperplaneTableEntry(const std::bitset<NbrPoints>& hyperplane,
		                                             const std::vector<HyperplaneTableEntry>& precedent_table) const noexcept {
			HyperplaneTableEntry entry;

			if constexpr (!OrderOfPoints) {
				entry.nbrPoints = static_cast<unsigned int>(hyperplane.count());

				std::vector<std::bitset<NbrPoints>> includedLines;
				for (const std::bitset<NbrPoints>& line : m_geometryLines) {
					if ((line & hyperplane) == line) {
						includedLines.push_back(line);
					}
				}

				entry.nbrLines = static_cast<unsigned int>(includedLines.size());
			} else {
				entry = getHyperplaneTableEntry<OrderOfPoints>(hyperplane);
			}

			entry.subgeometries.resize(m_subGeometriesMasks.size());

			{
				std::size_t i = 0;
				for (const auto& direction_masks : m_subGeometriesMasks) {
					for (const auto& mask : direction_masks) {
						std::size_t nbr_points = (hyperplane & mask).count();

						std::vector<HyperplaneTableEntry>::const_iterator it = std::find_if(precedent_table.begin(), precedent_table.end(),
						                                                                    [&nbr_points](const HyperplaneTableEntry& e) {
							                                                                    return e.nbrPoints == nbr_points;
						                                                                    });

						if (it == precedent_table.cend()) {
							++(entry.subgeometries[i][-1]);
						} else {
							++(entry.subgeometries[i][static_cast<long long int>(std::distance(precedent_table.cbegin(), it))]);
						}
					}

					++i;
				}
			}

			return entry;
		}

		template <bool OrderOfPoints>
		std::vector<HyperplaneTableEntry> makeHyperplaneTable(const std::vector<std::bitset<NbrPoints>>& vPoints) const noexcept {
			std::vector<HyperplaneTableEntry> entries;

			for (const auto& vPoint : vPoints) {
				HyperplaneTableEntry entry = getHyperplaneTableEntry<OrderOfPoints>(vPoint);

				std::vector<HyperplaneTableEntry>::iterator it = std::find(entries.begin(), entries.end(), entry);
				if (it == entries.end()) {
					entry.count = 1;
					entries.push_back(entry);
				} else {
					++(it->count);
				}
			}

			return entries;
		}

		template <bool OrderOfPoints>
		std::vector<HyperplaneTableEntry> makeHyperplaneTable(const std::vector<std::bitset<NbrPoints>>& vPoints,
		                                                      const std::vector<HyperplaneTableEntry>& precedent_table) const noexcept {
			std::vector<HyperplaneTableEntry> entries;

			for (const auto& vPoint : vPoints) {
				HyperplaneTableEntry entry = getHyperplaneTableEntry<OrderOfPoints>(vPoint, precedent_table);

				// Check if entry already exist
				std::vector<HyperplaneTableEntry>::iterator it = std::find(entries.begin(), entries.end(), entry);
				if (it == entries.end()) {
					entry.count = 1;
					entries.push_back(entry);
				} else {
					++(it->count);
				}
			}

			return entries;
		}


		VeldkampLineTableEntry
		makeLinesTableEntry(bool isProjective, const std::array<unsigned int, NbrPointsPerLine>& line,
		                    const std::vector<std::bitset<NbrPoints>>& vPoints,
		                    const std::vector<HyperplaneTableEntry>& points_table) const noexcept {
			VeldkampLineTableEntry entry;
			entry.isProjective = isProjective;

			std::bitset<NbrPoints> kernel = vPoints[line[0]] & vPoints[line[1]];
			entry.coreNbrPoints = kernel.count();
			entry.coreNbrLines = 0;
			for (const std::bitset<NbrPoints>& geometryLine : m_geometryLines) {
				if ((kernel & geometryLine) == geometryLine) {
					++entry.coreNbrLines;
				}
			}


			for (unsigned int i = 0; i < NbrPointsPerLine; ++i) {
				const size_t nbr_points = vPoints[line[i]].count();
				std::vector<HyperplaneTableEntry>::const_iterator it = std::find_if(points_table.begin(), points_table.end(),
				                                                                    [&nbr_points](const HyperplaneTableEntry& e) {
					                                                                    return e.nbrPoints == nbr_points;
				                                                                    });

				if (it == points_table.end()) {
					IMPOSSIBLE;
				} else {
					++entry.pointsType[static_cast<long long int>(std::distance(points_table.cbegin(), it))];
				}
			}

			return entry;
		}

		std::vector<VeldkampLineTableEntry> makeVeldkampLinesTable(VeldkampLines<NbrPointsPerLine>& vLines,
		                                                                             const std::vector<std::bitset<NbrPoints>>& vPoints,
		                                                                             const std::vector<HyperplaneTableEntry>& points_table) const noexcept {
			static_assert(Dimension < 4, "Points type determination only work for Dimension < 4");

			const auto makeEntries = [&](std::vector<VeldkampLineTableEntry>& entries,
			                             const std::vector<std::array<unsigned int, NbrPointsPerLine>> lines, bool isProjective) {
				for (const std::array<unsigned int, NbrPointsPerLine>& line : lines) {
					VeldkampLineTableEntry entry = makeLinesTableEntry(isProjective, line, vPoints, points_table);

					typename std::vector<VeldkampLineTableEntry>::iterator
							it = std::find(entries.begin(), entries.end(), entry);
					if (it == entries.end()) {
						entry.count = 1;
						entries.push_back(entry);
					} else {
						++(it->count);
					}
				}
			};

			std::vector<VeldkampLineTableEntry> entries;
			makeEntries(entries, vLines.projectives, true);
			makeEntries(entries, vLines.exceptional, false);
			return entries;
		}

		std::vector<VeldkampLineTableEntryWithLines<Dimension, NbrPointsPerLine>> makeVeldkampLinesTableWithLines(
		  const VeldkampLines<NbrPointsPerLine>& vLines,
		  const std::vector<std::bitset<NbrPoints>>& vPoints,
		  const std::vector<HyperplaneTableEntry>& points_table) const noexcept {

			static_assert(Dimension < 4, "Points type determination only work for Dimension < 4");

			const auto makeEntries =
			  [&](std::vector<VeldkampLineTableEntryWithLines<Dimension, NbrPointsPerLine>>& entries,
			      const std::vector<std::array<unsigned int, NbrPointsPerLine>> lines,
			      bool isProjective) {

				for (const std::array<unsigned int, NbrPointsPerLine>& line : lines) {
					VeldkampLineTableEntry entry = makeLinesTableEntry(isProjective, line, vPoints, points_table);

					const ptrdiff_t pos = std::find_if(entries.begin(), entries.end(), [&entry](const VeldkampLineTableEntryWithLines<Dimension, NbrPointsPerLine>& oentry){
						return oentry.entry == entry;
					}) - entries.begin();
					if (pos >= static_cast<ptrdiff_t>(entries.size())) {
						entry.count = 1;
						VeldkampLineTableEntryWithLines<Dimension, NbrPointsPerLine> entrywl(entry);
						entrywl.lines.push_back(line);
						entries.push_back(std::move(entrywl));
					} else {
						using size_type = typename std::vector<VeldkampLineTableEntryWithLines<Dimension, NbrPointsPerLine>>::size_type;
						++(entries[static_cast<size_type>(pos)].entry.count);
						entries[static_cast<size_type>(pos)].lines.push_back(line);
					}
				}
			};

			std::vector<VeldkampLineTableEntryWithLines<Dimension, NbrPointsPerLine>> entries;
			makeEntries(entries, vLines.projectives, true);
			makeEntries(entries, vLines.exceptional, false);
			return entries;
		}

	private:

		template <typename T>
		void permutations(T index, T bits, T number, std::vector<std::bitset<NbrPoints>>& hyperplanes) const noexcept {
			if (index == 0) {
				if (bits == 0) {
					std::bitset<NbrPoints> potentialHyperplane(number);
					if (isHyperplane(potentialHyperplane)) {
						hyperplanes.push_back(std::move(potentialHyperplane));
					}
				}
				return;
			}

			if (index - 1 >= bits) {
				permutations(index - 1, bits, number, hyperplanes);
			}

			if (bits > 0) {
				permutations(index - 1, bits - 1, number | T((1 << (index - 1))), hyperplanes);
			}
		}

		void computeMasks() {
			if constexpr (Dimension == 1) { // no sub dimensions in dimension1
				return;
			}

			std::array<std::bitset<NbrPoints>, Dimension> gen_lines; // lines starting from 0, like a canonical base
			std::array<std::array<size_t, NbrPointsPerLine>, Dimension> gen_lines_indexes; // indexes of the points of the previous lines

			// Fill gen_lines and gen_lines_indexes
			for (size_t i = 0; i < gen_lines.size(); ++i) {
				for (size_t j = 0; j < NbrPointsPerLine; ++j) {

					gen_lines_indexes[i][j] = j * math::pow(NbrPointsPerLine, i);
					gen_lines[i][j * math::pow(NbrPointsPerLine, i)] = 1;
				}
			}

			// Generate masks:
			// In the loop:
			// - take dimension-1 lines
			// - A = first line
			// - For each lines (except first):
			//     - A += A shifted along the line
			// - A is a mask
			// - A shifted along the line not taken at first step generate other masks
			for (size_t ignored_line = 0;
			     ignored_line < Dimension;
			     ++ignored_line) { // to take dimension-1 lines, we choose an ignored line
				size_t line_index = static_cast<size_t>(!ignored_line); // current line: first non-ignored line
				std::bitset<NbrPoints> gen_line = gen_lines[line_index];

				// First mask
				m_subGeometriesMasks[ignored_line][0] = gen_line;
				for (size_t line_shift_index = 0; line_shift_index < Dimension; ++line_shift_index) {
					if (line_shift_index == line_index || line_shift_index == ignored_line) {
						continue;
					}

					for (size_t i = 1; i < NbrPointsPerLine; ++i) {
						m_subGeometriesMasks[ignored_line][0] |= (gen_line << gen_lines_indexes[line_shift_index][i]);
					}
					gen_line = m_subGeometriesMasks[ignored_line][0];
				}

				// Shifts of the first mask
				for (size_t submask = 1; submask < NbrPointsPerLine; ++submask) {
					m_subGeometriesMasks[ignored_line][submask] =
							m_subGeometriesMasks[ignored_line][submask - 1] << gen_lines_indexes[ignored_line][1];
				}
			}
		}

		std::array<std::bitset<NbrPoints>, NbrLines> m_geometryLines;
		std::array<std::array<unsigned int, TensorSize>, NbrPoints> m_geometryPoints;

		std::array<std::array<std::bitset<NbrPoints>, NbrPointsPerLine>, Dimension> m_subGeometriesMasks;
	};

	template <std::size_t NbrPointsPerLine>
	VeldkampLines<NbrPointsPerLine>::VeldkampLines(std::vector<std::array<unsigned int, NbrPointsPerLine>>&& exceptional_lines,
	                                               std::vector<std::array<unsigned int, NbrPointsPerLine>>&& projectives_lines) noexcept
			: exceptional(std::move(exceptional_lines))
			, projectives(std::move(projectives_lines)) {
	}
}

#endif //HYPERPLANEFINDER_POINTGEOMETRY_HPP
