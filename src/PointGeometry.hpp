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

#include <CombinationGenerator.hpp>

namespace math {
	template <typename T>
	inline constexpr T pow(T base, T exponent) {
		return exponent == 0 ? 1 : base * pow(base, exponent-1);
	}
}

namespace std {
	template<size_t N>
	bool operator<(const std::bitset<N>& x, const std::bitset<N>& y) {
		for (size_t i = N-1; i--;) {
			if (x[i] ^ y[i]) {
				return y[i];
			}
		}

		return false;
	}
}

namespace segre {
	struct Entry {
		Entry()
				: nbrPoints{0}
				, nbrLines{0}
				, pointsOfOrder{}
				, count{0} {
		}

		template <typename Map>
		bool map_compare (Map const &lhs, Map const &rhs) const {
			return lhs.size() == rhs.size()
			       && std::equal(lhs.begin(), lhs.end(),
			                     rhs.begin());
		}

		bool operator==(const Entry& entry) const {
			return nbrPoints == entry.nbrPoints &&
			       nbrLines == entry.nbrLines &&
			       map_compare(pointsOfOrder, entry.pointsOfOrder);
		}

		friend std::ostream& operator<<(std::ostream &os, const Entry& entry);

		unsigned int nbrPoints;
		unsigned int nbrLines;
		std::map<unsigned int, unsigned int> pointsOfOrder;
		size_t count;
	};

	template<size_t N1,size_t N2>
	inline std::bitset<N1> copyBitset(const std::bitset<N2>& bs2) {
		std::bitset<N1> bs1;
		for(size_t i = 0; i < N2; i++) {
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
				: m_lines(std::move(lines))
				, m_tPts(TENSOR_2D) {
		}

		explicit PointGeometry(std::array<std::bitset<NbrPoints>, NbrLines>&& lines,
		                       std::array<std::array<unsigned int, TensorSize>, NbrPoints>&& tensors) noexcept
				: m_lines(std::move(lines))
				, m_tPts(std::move(tensors)) {
		}

		std::vector<std::bitset<NbrPoints>> findHyperplanes() const noexcept {
			std::vector<std::bitset<NbrPoints>> hyperplanes;

			for (unsigned int j = 2; j < NbrPoints; ++j) {
				permutations<std::uint64_t>(NbrPoints, j, 0, hyperplanes);
			}

			return std::move(hyperplanes);
		}

		bool isHyperplane(const std::bitset<NbrPoints>& potentialHyperplane) const noexcept {
			for (const auto& line : m_lines) {
				auto intersection = line & potentialHyperplane;
				auto intersectionSize = intersection.count();
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

		std::pair<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>>
		findVeldkampLinesDimension4(const std::vector<std::bitset<NbrPoints>>& veldkampPoints) const noexcept {
			std::vector<std::vector<unsigned int>> supposedExceptional;
			std::vector<std::vector<unsigned int>> projectiveLines;

			CombinationGenerator gen;
			gen.initialize(static_cast<unsigned int>(veldkampPoints.size()), 2);

			while (!gen.isFinished()) {
				const std::vector<unsigned int>& currentCombination = gen.nextCombination();

				std::vector<std::pair<std::bitset<NbrPoints>, unsigned int>> sameCore;
				const std::bitset<NbrPoints>& h1 = veldkampPoints[currentCombination[0]];
				const std::bitset<NbrPoints>& h2 = veldkampPoints[currentCombination[1]];

				const std::bitset<NbrPoints> intersection12 = h1 & h2;

				for (size_t i = 0, n = veldkampPoints.size(); i < n; ++i) {
					const std::bitset<NbrPoints> intersection1i = h1 & veldkampPoints[i];
					const std::bitset<NbrPoints> intersection2i = h2 & veldkampPoints[i];

					if (intersection12 == intersection1i
							&& intersection12 == intersection2i
							&& intersection1i == intersection2i) {
						sameCore.emplace_back(veldkampPoints[i], i);
					}
				}

				CombinationGenerator gen2;
				gen2.initialize(static_cast<unsigned int>(sameCore.size()), 2);

				while (!gen2.isFinished()) {
					const std::vector<unsigned int>& currentCombination2 = gen2.nextCombination();

					if (sameCore[currentCombination2[0]].second > currentCombination[0] && sameCore[currentCombination2[1]].second > currentCombination[1]
							&& sameCore[currentCombination2[1]].second > currentCombination[0] && sameCore[currentCombination2[0]].second > currentCombination[1]){
						const std::bitset<NbrPoints>& ha = sameCore[currentCombination2[0]].first;
						const std::bitset<NbrPoints>& hb = sameCore[currentCombination2[1]].first;

						const std::bitset<NbrPoints> intersectionAB = ha & hb;

						if (intersection12 == intersectionAB) {
							if (sameCore.size() == 2) {
								projectiveLines.emplace_back(std::initializer_list<unsigned int>{currentCombination[0],
								                              currentCombination[1],
								                              sameCore[currentCombination2[0]].second,
								                              sameCore[currentCombination2[1]].second});
							} else {
								supposedExceptional.emplace_back(std::initializer_list<unsigned int>{currentCombination[0],
								                               currentCombination[1],
								                               sameCore[currentCombination2[0]].second,
								                               sameCore[currentCombination2[1]].second});
							}
						}
					}
				}
			}

			return std::make_pair(supposedExceptional, projectiveLines);
		}

		void distinguishVeldkampLines(std::pair<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>>& vLines,
		                              const std::vector<std::bitset<NbrPoints>>& vPoints,
		                              const PointGeometry<Dimension + 1, NbrPointsPerLine, math::pow(NbrPointsPerLine, Dimension) * (1 + Dimension)>& nextGeometry) const {
			std::vector<size_t> toRemove;
			constexpr size_t NewNbrPoints = math::pow(NbrPointsPerLine, Dimension + 1);

			for (size_t index = 0; index < vLines.first.size(); ++index) {
				std::bitset<NewNbrPoints> hyperplane;
				for (size_t i = 0; i < vLines.first[index].size(); ++i) {
					hyperplane |= copyBitset<NewNbrPoints>(vPoints[vLines.first[index][i]]) <<= (i * NbrPoints);
				}

				if (getRank(nextGeometry.buildMatrix(hyperplane)) < math::pow(2UL, Dimension + 1)) {
					toRemove.push_back(index);
				}
			}

			vLines.second.reserve(vLines.second.size() + toRemove.size());

			while (!toRemove.empty()) {
				vLines.second.push_back(vLines.first[toRemove.back()]);
				vLines.first.erase(std::next(vLines.first.begin(), static_cast<unsigned int>(toRemove.back())));

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

		decltype(auto) computeHyperplanes(const std::vector<std::bitset<NbrPoints>>& veldkampPoints, const std::vector<std::vector<unsigned int>>& pVLines) {
			constexpr size_t NewNbrPoints = math::pow(NbrPointsPerLine, Dimension + 1);

			std::vector<std::bitset<NewNbrPoints>> hyperplanes;

			for (size_t i = 0; i < pVLines.size(); ++i) {
				std::array<std::bitset<NbrPoints>, NbrPointsPerLine> hypers = getHyperplanes(veldkampPoints, pVLines[i]);
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

		static std::array<std::bitset<NbrPoints>, NbrPointsPerLine> getHyperplanes(const std::vector<std::bitset<NbrPoints>>& veldkampPoints,
		                                                 const std::vector<unsigned int>& veldkampLine) {
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

			std::generate(result.begin(), result.end(), [this, i = 0, j = 0] () mutable -> auto {
				auto bitset = copyBitset<NewNbrPoints>(m_lines[j]) <<= static_cast<long unsigned int>(NbrPoints * i);
				++j;
				if (j % NbrLines == 0) {
					j = 0;
					++i;
				}

				return bitset;
			});

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
							pts[i * NbrPoints + j][k * TensorSize + l] = TENSOR_2D[i][k] * m_tPts[j][l] % 3;
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
					matrix.push_back(m_tPts[i]);
				}
			}

			return matrix;
		}

		Entry getEntry(const std::bitset<NbrPoints>& hyperplane) const noexcept {
			Entry entry;
			entry.nbrPoints = static_cast<unsigned int>(hyperplane.count());

			std::vector<std::bitset<NbrPoints>> includedLines;
			for (const auto& line : m_lines) {
				if ((line & hyperplane) == line) {
					includedLines.push_back(line);
				}
			}

			entry.nbrLines = static_cast<unsigned int>(includedLines.size());

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

			return entry;
		}

		std::vector<Entry> makeTable(const std::vector<std::bitset<NbrPoints>>& vPoints) const noexcept {
			std::vector<Entry> entries;

			for (const auto& vPoint : vPoints) {
				Entry entry = getEntry(vPoint);

				std::vector<Entry>::iterator it = std::find(entries.begin(), entries.end(), entry);
				if (it == entries.end()) {
					entry.count = 1;
					entries.push_back(entry);
				} else {
					++(it->count);
				}
			}

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

			if (index-1 >= bits) {
				permutations(index-1, bits, number, hyperplanes);
			}

			if (bits > 0) {
				permutations(index-1, bits-1, number | T((1 << (index-1))), hyperplanes);
			}
		}

		std::array<std::bitset<NbrPoints>, NbrLines> m_lines;
		std::array<std::array<unsigned int, TensorSize>, NbrPoints> m_tPts;
	};

	std::ostream& operator<<(std::ostream &os, const Entry& entry) {
		os << "Entry{" <<
		   "Ps: " << entry.nbrPoints <<
		   ", Ls: " << entry.nbrLines <<
		   ", Order={";

		for(auto iterator = entry.pointsOfOrder.cbegin(); iterator != entry.pointsOfOrder.cend();) {
			os << iterator->first << "=" << iterator->second;
			if (++iterator != entry.pointsOfOrder.cend()) {
				os << ", ";
			}
		}

		os << "}, Crd: " << entry.count << '}';

		return os;
	}
}

#endif //HYPERPLANEFINDER_POINTGEOMETRY_HPP
