#include <iostream>
#include <iterator>
#include <utility>
#include <algorithm>
#include <vector>
#include <chrono>

#include <nlohmann/json.hpp>
#include <inja.hpp>

#include "PointGeometry.hpp"
#include "LatexPrinter.hpp"
#include "HyperplanesUtility.hpp"
#include "VeldkampLinesUtility.hpp"

using json = nlohmann::json;

constexpr size_t PPL = 4; // Points Per Lines
constexpr bool COMPUTE_AND_PRINT_POINTS_ORDER = false;
constexpr bool PRINT_SUBGEOMETRIES = true;

template<int N>
using VPoints = std::vector<std::bitset<math::pow(PPL,N)>>;

template<int>
using VLines = segre::VeldkampLines<PPL>;

VPoints<3> bruteForceD3Hyperplanes(const segre::PointGeometry<3, PPL, 48>& geometry3, const VPoints<2>& vPoints2){
	constexpr auto NbrPointsD3 = math::pow(PPL,3);
	constexpr auto NbrPointsD2 = math::pow(PPL,2);
	VPoints<3> vPoints3;

	for (const auto& vPoint : vPoints2) {
		std::bitset<NbrPointsD3> newH;

		newH |= segre::copyBitset<NbrPointsD3, NbrPointsD2>(vPoint);
		newH |= segre::copyBitset<NbrPointsD3, NbrPointsD2>(vPoint) << NbrPointsD2;
		newH |= segre::copyBitset<NbrPointsD3, NbrPointsD2>(vPoint) << (2 * NbrPointsD2);
		newH |= segre::copyBitset<NbrPointsD3, NbrPointsD2>(~vPoint) << (3 * NbrPointsD2);

		if (geometry3.isHyperplane(newH)) {
			vPoints3.push_back(newH);
		}
	}

	CombinationGenerator cg;
	cg.initialize(static_cast<unsigned int>(vPoints2.size()), PPL);

	while (!cg.isFinished()) {

		std::bitset<NbrPointsD3> newH;

		{
			unsigned int i = 0;
			for (unsigned int combination : cg.nextCombination()) {
				newH |= segre::copyBitset<NbrPointsD3, NbrPointsD2>(vPoints2[combination]) << (i * NbrPointsD2);
				++i;
			}

			if (geometry3.isHyperplane(newH)) {
				vPoints3.push_back(newH);
			}
		}

		if (cg.getNumLeft() % 100'000'000 == 0) {
			std::cout << cg.getNumLeft() << '\n';
		}
	}

	return vPoints3;
}

int main() {
	const auto time_start = std::chrono::system_clock::now();

	std::array<std::bitset<PPL>, 1> lines;
	lines[0] = std::bitset<PPL>(math::pow(2UL,PPL) - 1);

	segre::PointGeometry<1, PPL, 1> geometry1(std::move(lines));
	segre::PointGeometry<2, PPL, 8> geometry2(geometry1.computeCartesianProduct(), geometry1.buildTensorPoints());
	segre::PointGeometry<3, PPL, 48> geometry3(geometry2.computeCartesianProduct(), geometry2.buildTensorPoints());
	segre::PointGeometry<4, PPL, 256> geometry4(geometry3.computeCartesianProduct(), geometry3.buildTensorPoints());

	VPoints<2> vPoints2 = geometry2.findHyperplanesByBruteforce(); // brut force

	std::vector<segre::HyperplaneTableEntry> geometry2_hyp_table = geometry2.makeHyperplaneTable<COMPUTE_AND_PRINT_POINTS_ORDER>(vPoints2);
	std::sort(geometry2_hyp_table.begin(), geometry2_hyp_table.end(), [] (const segre::HyperplaneTableEntry& a, const segre::HyperplaneTableEntry& b) {
		return a.nbrPoints > b.nbrPoints;
	});

	VPoints<3> vPoints3 = bruteForceD3Hyperplanes(geometry3, vPoints2);

	std::vector<segre::HyperplaneTableEntry> geometry3_hyp_table = geometry3.makeHyperplaneTable<COMPUTE_AND_PRINT_POINTS_ORDER>(vPoints3);
	std::sort(geometry3_hyp_table.begin(), geometry3_hyp_table.end(), [] (const segre::HyperplaneTableEntry& a, const segre::HyperplaneTableEntry& b) {
		return a.nbrPoints > b.nbrPoints;
	});

	std::vector<segre::HyperplaneTableEntry> geometry3_hyp_table_with_points_order = geometry3.makeHyperplaneTable<true>(vPoints3);
	std::sort(geometry3_hyp_table_with_points_order.begin(), geometry3_hyp_table_with_points_order.end(), [] (const segre::HyperplaneTableEntry& a, const segre::HyperplaneTableEntry& b) {
		return a.nbrPoints > b.nbrPoints;
	});

	/*std::vector<segre::HyperplaneTableEntry> geometry3_hyp_table_with_subgeometries = geometry3.makeHyperplaneTable<COMPUTE_AND_PRINT_POINTS_ORDER>(vPoints3, geometry2_hyp_table);
	std::sort(geometry3_hyp_table_with_subgeometries.begin(), geometry3_hyp_table_with_subgeometries.end(), [] (const segre::HyperplaneTableEntry& a, const segre::HyperplaneTableEntry& b) {
		return a.nbrPoints > b.nbrPoints;
	});*/

	const auto time_end = std::chrono::system_clock::now();
	const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	std::cout << "Finished in " << static_cast<int>(elapsed.count()) << " seconds\n" << std::endl;

	std::cout << "\nDimension 2 points:\n";
	std::copy(geometry2_hyp_table.begin(), geometry2_hyp_table.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));

	std::cout << "\nDimension 3 points:\n";
	std::copy(geometry3_hyp_table.begin(), geometry3_hyp_table.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));

	LatexPrinter printer;
	printer.generateHyperplanesTable<COMPUTE_AND_PRINT_POINTS_ORDER,PRINT_SUBGEOMETRIES>(2, geometry2_hyp_table, 0);
	printer.generateHyperplanesTable<false,false>(3, geometry3_hyp_table, geometry2_hyp_table.size());
	printer.generateHyperplanesDiffTable<true, false>(3, geometry3_hyp_table, geometry3_hyp_table_with_points_order, geometry2_hyp_table.size());
	//printer.generateHyperplanesDiffTable<COMPUTE_AND_PRINT_POINTS_ORDER, true>(3, geometry3_hyp_table, geometry3_hyp_table_with_subgeometries, geometry2_hyp_table.size());

	return EXIT_SUCCESS;
}
