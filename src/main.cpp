#include <iostream>
#include <iterator>
#include <utility>
#include <algorithm>
#include <vector>
#include <chrono>

#include <math_utility.hpp>
#include <PointGeometry.hpp>
#include <LatexPrinter.hpp>
#include <EquationSolver.hpp>

#include <nlohmann/json.hpp>
#include <inja.hpp>

using json = nlohmann::json;

constexpr size_t PPL = 4; // Points Per Lines
constexpr bool CIMPUTE_AND_PRINT_POINTS_ORDER = true;
constexpr bool PRINT_SUBGEOMETRIES = true;

template <int N>
using VPoints = std::vector<std::bitset<math::pow(PPL,N)>>;

template<int>
using VLines = segre::VeldkampLines<PPL>;

template <int N>
using VLinesTableEntry = segre::VeldkampLineTableEntry<math::pow(PPL,N)>;

/*int main() {
	const auto time_start = std::chrono::system_clock::now();

	std::array<std::bitset<PPL>, 1> lines;
	lines[0] = std::bitset<PPL>(math::pow(2UL,PPL) - 1);

	segre::PointGeometry<1, PPL, 1> geometry1(std::move(lines));
	segre::PointGeometry<2, PPL, 8> geometry2(geometry1.computeCartesianProduct(), geometry1.buildTensorPoints());
	segre::PointGeometry<3, PPL, 48> geometry3(geometry2.computeCartesianProduct(), geometry2.buildTensorPoints());
	segre::PointGeometry<4, PPL, 256> geometry4(geometry3.computeCartesianProduct(), geometry3.buildTensorPoints());

	VPoints<2> vPoints2 = geometry2.findHyperplanesByBruteforce(); // brut force
	VLines<2> vLines2 = geometry2.computeVeldkampLines(vPoints2);
	geometry2.distinguishVeldkampLines(vLines2, vPoints2, geometry3);

	std::vector<segre::HyperplaneTableEntry> geometry2_hyp_table = geometry2.makeHyperplaneTable<CIMPUTE_AND_PRINT_POINTS_ORDER>(vPoints2);
	std::sort(geometry2_hyp_table.begin(), geometry2_hyp_table.end(), [] (const segre::HyperplaneTableEntry& a, const segre::HyperplaneTableEntry& b) {
		return a.nbrPoints > b.nbrPoints;
	});
	std::vector<VLinesTableEntry<2>> geometry2_lin_table = geometry2.makeVeldkampLinesTable(vLines2, vPoints2, geometry2_hyp_table);
	std::sort(geometry2_lin_table.begin(), geometry2_lin_table.end(), [](const VLinesTableEntry<2>& a, const VLinesTableEntry<2>& b){
		return std::make_tuple(a.isProjective, a.coreNbrPoints, a.coreNbrLines) < std::make_tuple(b.isProjective, b.coreNbrPoints, b.coreNbrLines);
	});

	VPoints<3> vPoints3 = geometry2.computeHyperplanesFromVeldkampLines(vPoints2, vLines2.projectives);
	VLines<3> vLines3 = geometry3.computeVeldkampLines(vPoints3);
	geometry3.distinguishVeldkampLines(vLines3, vPoints3, geometry4);

//	std::vector<segre::HyperplaneTableEntry> geometry3_hyp_table = geometry3.makeHyperplaneTable<CIMPUTE_AND_PRINT_POINTS_ORDER>(vPoints3, geometry2_hyp_table);
//	std::sort(geometry3_hyp_table.begin(), geometry3_hyp_table.end(), [] (const segre::HyperplaneTableEntry& a, const segre::HyperplaneTableEntry& b) {
//		return a.nbrPoints > b.nbrPoints;
//	});
//	std::vector<VLinesTableEntry<3>> geometry3_lin_table = geometry3.makeVeldkampLinesTable(vLines3, vPoints3, geometry3_hyp_table);
//	std::sort(geometry3_lin_table.begin(), geometry3_lin_table.end(), [](const VLinesTableEntry<3>& a, const VLinesTableEntry<3>& b){
//		return std::make_tuple(a.isProjective, a.coreNbrPoints, a.coreNbrLines) < std::make_tuple(b.isProjective, b.coreNbrPoints, b.coreNbrLines);
//	});

	VPoints<4> vPoints4 = geometry3.computeHyperplanesFromVeldkampLines(vPoints3, vLines3.projectives);
//	std::vector<segre::HyperplaneTableEntry> geometry4_hyp_table = geometry4.makeHyperplaneTable<CIMPUTE_AND_PRINT_POINTS_ORDER>(vPoints4, geometry3_hyp_table);
//
//	std::sort(geometry4_hyp_table.begin(), geometry4_hyp_table.end(), [] (const segre::HyperplaneTableEntry& a,
//	                                              const segre::HyperplaneTableEntry& b) {
//		return a.nbrPoints > b.nbrPoints;
//	});

	const auto time_end = std::chrono::system_clock::now();
	const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	std::cout << "Finished in " << static_cast<int>(elapsed.count()) << " seconds\n" << std::endl;

	std::size_t counter = 0;
	for (const auto& point : vPoints4) {
		std::string str = point.to_string();
		std::string rstr = std::string(str.rbegin(), str.rend());

		if (str == rstr) {
			++counter;
			//std::cout << str;
		}
	}
	std::cout << counter << std::endl;
//	std::cout << "\nDimension 2 points:\n";
//	std::copy(geometry2_hyp_table.begin(), geometry2_hyp_table.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));
//
//	std::cout << "\nDimension 2 lines:\n";
//	std::copy(geometry2_lin_table.begin(), geometry2_lin_table.end(), std::ostream_iterator<VLinesTableEntry<2>>(std::cout, "\n"));
//
//	std::cout << "\nDimension 3 points:\n";
//	std::copy(geometry3_hyp_table.begin(), geometry3_hyp_table.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));
//
//	std::cout << "\nDimension 3 lines:\n";
//	std::copy(geometry3_lin_table.begin(), geometry3_lin_table.end(), std::ostream_iterator<VLinesTableEntry<3>>(std::cout, "\n"));
//
//	std::cout << "\nDimension 4 points:\n";
//	std::copy(geometry4_hyp_table.begin(), geometry4_hyp_table.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));
//
//	LatexPrinter printer;
//	printer.generateHyperplanesTable<CIMPUTE_AND_PRINT_POINTS_ORDER,PRINT_SUBGEOMETRIES>(2, geometry2_hyp_table, 0);
//	printer.generateLinesTable(2, geometry2_lin_table, geometry2_hyp_table.size());
//	printer.generateHyperplanesTable<CIMPUTE_AND_PRINT_POINTS_ORDER,PRINT_SUBGEOMETRIES>(3, geometry3_hyp_table, geometry2_hyp_table.size());
//	printer.generateLinesTable(3, geometry3_lin_table, geometry3_hyp_table.size());
//	printer.generateHyperplanesTable<CIMPUTE_AND_PRINT_POINTS_ORDER,PRINT_SUBGEOMETRIES>(4, geometry4_hyp_table, geometry3_hyp_table.size());
//	printer.generateHyperplaneRepresentation<3, PPL>(geometry3_lin_table[14 - 1].core);

	return EXIT_SUCCESS;
}*/

int main() {
//	std::array<std::size_t, 16> orbit39{};
//	orbit39[0] = 0;
//	orbit39[1] = 0;
//	orbit39[2] = 0;
//	orbit39[3] = 1;
//	orbit39[4] = 0;
//	orbit39[5] = 1;
//	orbit39[6] = 1;
//	orbit39[7] = 0;
//	orbit39[8] = 1;
//	orbit39[9] = 0;
//	orbit39[10] = 0;
//	orbit39[11] = 0;
//	orbit39[12] = 0;
//	orbit39[13] = 1;
//	orbit39[14] = 2;
//	orbit39[15] = 1;

//	std::array<std::size_t, 16> orbit44{};
//	orbit44[0] = 0;
//	orbit44[1] = 0;
//	orbit44[2] = 0;
//	orbit44[3] = 1;
//	orbit44[4] = 0;
//	orbit44[5] = 1;
//	orbit44[6] = 1;
//	orbit44[7] = 0;
//	orbit44[8] = 1;
//	orbit44[9] = 0;
//	orbit44[10] = 0;
//	orbit44[11] = 0;
//	orbit44[12] = 1;
//	orbit44[13] = 0;
//	orbit44[14] = 2;
//	orbit44[15] = 2;

//	std::array<std::size_t, 16> orbit36{};
//	orbit36[0] = 0;
//	orbit36[1] = 0;
//	orbit36[2] = 0;
//	orbit36[3] = 1;
//	orbit36[4] = 0;
//	orbit36[5] = 1;
//	orbit36[6] = 1;
//	orbit36[7] = 0;
//	orbit36[8] = 1;
//	orbit36[9] = 0;
//	orbit36[10] = 0;
//	orbit36[11] = 0;
//	orbit36[12] = 0;
//	orbit36[13] = 0;
//	orbit36[14] = 1;
//	orbit36[15] = 1;

//	std::array<std::size_t, 16> orbit45{};
//	orbit45[0] = 0;
//	orbit45[1] = 0;
//	orbit45[2] = 0;
//	orbit45[3] = 1;
//	orbit45[4] = 0;
//	orbit45[5] = 1;
//	orbit45[6] = 1;
//	orbit45[7] = 0;
//	orbit45[8] = 1;
//	orbit45[9] = 0;
//	orbit45[10] = 0;
//	orbit45[11] = 0;
//	orbit45[12] = 1;
//	orbit45[13] = 1;
//	orbit45[14] = 2;
//	orbit45[15] = 0;

	std::array<std::size_t, 16> orbit37{};
	orbit37[0] = 0;
	orbit37[1] = 0;
	orbit37[2] = 0;
	orbit37[3] = 1;
	orbit37[4] = 0;
	orbit37[5] = 1;
	orbit37[6] = 1;
	orbit37[7] = 0;
	orbit37[8] = 1;
	orbit37[9] = 0;
	orbit37[10] = 0;
	orbit37[11] = 0;
	orbit37[12] = 0;
	orbit37[13] = 0;
	orbit37[14] = 2;
	orbit37[15] = 1;

//	std::array<std::size_t, 16> orbit42{};
//	orbit42[0] = 0;
//	orbit42[1] = 0;
//	orbit42[2] = 0;
//	orbit42[3] = 1;
//	orbit42[4] = 0;
//	orbit42[5] = 1;
//	orbit42[6] = 1;
//	orbit42[7] = 0;
//	orbit42[8] = 1;
//	orbit42[9] = 0;
//	orbit42[10] = 0;
//	orbit42[11] = 0;
//	orbit42[12] = 1;
//	orbit42[13] = 0;
//	orbit42[14] = 1;
//	orbit42[15] = 1;

//	std::array<std::size_t, 16> orbit30{};
//	orbit30[0] = 0;
//	orbit30[1] = 0;
//	orbit30[2] = 0;
//	orbit30[3] = 0;
//	orbit30[4] = 0;
//	orbit30[5] = 1;
//	orbit30[6] = 1;
//	orbit30[7] = 0;
//	orbit30[8] = 1;
//	orbit30[9] = 0;
//	orbit30[10] = 1;
//	orbit30[11] = 1;
//	orbit30[12] = 0;
//	orbit30[13] = 0;
//	orbit30[14] = 0;
//	orbit30[15] = 0;

//  std::array<std::size_t, 16> orbit49{};
//  orbit49[0] = 0;
//  orbit49[1] = 1;
//  orbit49[2] = 1;
//  orbit49[3] = 0;
//  orbit49[4] = 1;
//  orbit49[5] = 0;
//  orbit49[6] = 0;
//  orbit49[7] = 2;
//  orbit49[8] = 1;
//  orbit49[9] = 0;
//  orbit49[10] = 0;
//  orbit49[11] = 2;
//  orbit49[12] = 1;
//  orbit49[13] = 1;
//  orbit49[14] = 1;
//  orbit49[15] = 2;

//  std::array<std::size_t, 16> orbit47{};
//  orbit47[0] = 0;
//  orbit47[1] = 1;
//  orbit47[2] = 1;
//  orbit47[3] = 0;
//  orbit47[4] = 1;
//  orbit47[5] = 0;
//  orbit47[6] = 0;
//  orbit47[7] = 2;
//  orbit47[8] = 1;
//  orbit47[9] = 0;
//  orbit47[10] = 0;
//  orbit47[11] = 2;
//  orbit47[12] = 0;
//  orbit47[13] = 2;
//  orbit47[14] = 2;
//  orbit47[15] = 0;

//  std::array<std::size_t, 16> orbit{};
//  orbit[0] = 0;
//  orbit[1] = 0;
//  orbit[2] = 0;
//  orbit[3] = 1;
//  orbit[4] = 0;
//  orbit[5] = 1;
//  orbit[6] = 1;
//  orbit[7] = 0;
//  orbit[8] = 0;
//  orbit[9] = 1;
//  orbit[10] = 1;
//  orbit[11] = 0;
//  orbit[12] = 1;
//  orbit[13] = 0;
//  orbit[14] = 0;
//  orbit[15] = 1;

	auto solutions = segre::resolveEquation<4>(orbit37);
    std::cout << solutions.size() << '\n';

    for (auto solution : solutions) {
        for (auto couple : solution) {
            std::cout << '{' << +couple[0] << ", " << +couple[1] << "}, ";
        }

        std::cout << "\n";
    }

    auto hyperplane = segre::solutionsToHyperplane<4, 256>(solutions);
    //hyperplane = ~hyperplane;
    for (std::size_t i = 0; i < 256; ++i) {
        if (hyperplane[i]) {
            std::cout << i << " ";
        }
    }
	std::cout << '\n';

    std::string hyperplaneStr = hyperplane.to_string();
	std::cout << (hyperplaneStr == std::string(hyperplaneStr.rbegin(), hyperplaneStr.rend())) << '\n';
	std::cout << "count=" << hyperplane.count() << '\n';

    std::array<std::bitset<PPL>, 1> lines;
    lines[0] = std::bitset<PPL>(math::pow(2UL,PPL) - 1);

    segre::PointGeometry<1, PPL, 1> geometry1(std::move(lines));
    segre::PointGeometry<2, PPL, 8> geometry2(geometry1.computeCartesianProduct(), geometry1.buildTensorPoints());
    segre::PointGeometry<3, PPL, 48> geometry3(geometry2.computeCartesianProduct(), geometry2.buildTensorPoints());
    segre::PointGeometry<4, PPL, 256> geometry4(geometry3.computeCartesianProduct(), geometry3.buildTensorPoints());

    std::cout << std::boolalpha << geometry4.isHyperplane(hyperplane) << '\n';
    std::cout << geometry4.getHyperplaneTableEntry<true>(hyperplane) << '\n';

	LatexPrinter printer;
	printer.generateHyperplaneRepresentation<4, PPL>(hyperplane);

	for (std::size_t j = 0; j < PPL; ++j) {
		auto subHyperplanes = geometry4.getSubHyperplanesFromVL(geometry4.getVeldkampLineInBase4(hyperplane, hyperplane.count()), j);

		std::bitset<64> kernel{};
		kernel = ~kernel;

		for (const auto& subHyperplane : subHyperplanes) {
			kernel &= subHyperplane;
		}

		std::vector<std::pair<unsigned int, std::array<unsigned int, 3>>> pointsInBase4;

		for (std::size_t i = 0; i < kernel.size(); ++i) {
			if (kernel[i]) {
				pointsInBase4.emplace_back(i, geometry4.base10tobasePPL<3>(i));
			}
		}

		std::vector<std::vector<std::size_t>> distanceMatrix(pointsInBase4.size(), std::vector<std::size_t>(pointsInBase4.size(), 0));

		for (std::size_t pointIndex1 = 0; pointIndex1 < pointsInBase4.size(); ++pointIndex1) {
			for (std::size_t pointIndex2 = 0; pointIndex2 < pointsInBase4.size(); ++pointIndex2) {
				if (pointIndex1 == pointIndex2) {
					continue;
				}

				std::size_t distance = 0;

				for (std::size_t k = 0; k < 3; ++k) {
					if (pointsInBase4[pointIndex1].second[k] != pointsInBase4[pointIndex2].second[k]) {
						distance += 1;
					}
				}

				distanceMatrix[pointIndex1][pointIndex2] = distance;
			}
		}

		for (const auto& line : distanceMatrix) {
			for (const auto& e : line) {
				std::cout << e << " ";
			}
			std::cout << '\n';
		}

		printer.generateHyperplaneRepresentation<3, PPL>(kernel);
		std::cout << "\n\n";
	}

    return EXIT_SUCCESS;
}
