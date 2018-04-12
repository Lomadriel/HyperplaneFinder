#include <iostream>
#include <iterator>
#include <utility>
#include <algorithm>
#include <vector>
#include <chrono>

#include <PointGeometry.hpp>
#include <LatexPrinter.hpp>

#include <nlohmann/json.hpp>
#include <inja.hpp>

using json = nlohmann::json;

constexpr size_t PPL = 4; // Points Per Lines
constexpr bool CIMPUTE_AND_PRINT_POINTS_ORDER = true;
constexpr bool PRINT_SUBGEOMETRIES = true;

template<int N>
using VPoints = std::vector<std::bitset<math::pow(PPL,N)>>;

template<int>
using VLines = segre::VeldkampLines<PPL>;

int GRID() {
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
	std::vector<segre::VeldkampLineTableEntry> geometry2_lin_table = geometry2.makeVeldkampLinesTable(vLines2, vPoints2, geometry2_hyp_table);
	std::sort(geometry2_lin_table.begin(), geometry2_lin_table.end(), [](const segre::VeldkampLineTableEntry& a, const segre::VeldkampLineTableEntry& b){
		return std::make_tuple(a.isProjective, a.coreNbrPoints, a.coreNbrLines) < std::make_tuple(b.isProjective, b.coreNbrPoints, b.coreNbrLines);
	});

	VPoints<3> vPoints3 = geometry2.computeHyperplanesFromVeldkampLines(vPoints2, vLines2.projectives);
	VLines<3> vLines3 = geometry3.computeVeldkampLines(vPoints3);
	geometry3.distinguishVeldkampLines(vLines3, vPoints3, geometry4);

	std::vector<segre::HyperplaneTableEntry> geometry3_hyp_table = geometry3.makeHyperplaneTable<CIMPUTE_AND_PRINT_POINTS_ORDER>(vPoints3, geometry2_hyp_table);
	std::sort(geometry3_hyp_table.begin(), geometry3_hyp_table.end(), [] (const segre::HyperplaneTableEntry& a, const segre::HyperplaneTableEntry& b) {
		return a.nbrPoints > b.nbrPoints;
	});
	std::vector<segre::VeldkampLineTableEntry> geometry3_lin_table = geometry3.makeVeldkampLinesTable(vLines3, vPoints3, geometry3_hyp_table);
	std::sort(geometry3_lin_table.begin(), geometry3_lin_table.end(), [](const segre::VeldkampLineTableEntry& a, const segre::VeldkampLineTableEntry& b){
		return std::make_tuple(a.isProjective, a.coreNbrPoints, a.coreNbrLines) < std::make_tuple(b.isProjective, b.coreNbrPoints, b.coreNbrLines);
	});

	VPoints<4> vPoints4 = geometry3.computeHyperplanesFromVeldkampLines(vPoints3, vLines3.projectives);
	std::vector<segre::HyperplaneTableEntry> geometry4_hyp_table = geometry4.makeHyperplaneTable<CIMPUTE_AND_PRINT_POINTS_ORDER>(vPoints4, geometry3_hyp_table);

	std::sort(geometry4_hyp_table.begin(), geometry4_hyp_table.end(), [] (const segre::HyperplaneTableEntry& a,
	                                                                      const segre::HyperplaneTableEntry& b) {
		return a.nbrPoints > b.nbrPoints;
	});

	const auto time_end = std::chrono::system_clock::now();
	const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	std::cout << "Finished in " << static_cast<int>(elapsed.count()) << " seconds\n" << std::endl;

	std::cout << "\nDimension 2 points:\n";
	std::copy(geometry2_hyp_table.begin(), geometry2_hyp_table.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));

	std::cout << "\nDimension 2 lines:\n";
	std::copy(geometry2_lin_table.begin(), geometry2_lin_table.end(), std::ostream_iterator<segre::VeldkampLineTableEntry>(std::cout, "\n"));

	std::cout << "\nDimension 3 points:\n";
	std::copy(geometry3_hyp_table.begin(), geometry3_hyp_table.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));

	std::cout << "\nDimension 3 lines:\n";
	std::copy(geometry3_lin_table.begin(), geometry3_lin_table.end(), std::ostream_iterator<segre::VeldkampLineTableEntry>(std::cout, "\n"));

	std::cout << "\nDimension 4 points:\n";
	std::copy(geometry4_hyp_table.begin(), geometry4_hyp_table.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));

	LatexPrinter printer;
	printer.generateHyperplanesTable<CIMPUTE_AND_PRINT_POINTS_ORDER,PRINT_SUBGEOMETRIES>(2, geometry2_hyp_table, 0);
	printer.generateLinesTable(2, geometry2_lin_table, geometry2_hyp_table.size());
	printer.generateHyperplanesTable<CIMPUTE_AND_PRINT_POINTS_ORDER,PRINT_SUBGEOMETRIES>(3, geometry3_hyp_table, geometry2_hyp_table.size());
	printer.generateLinesTable(3, geometry3_lin_table, geometry3_hyp_table.size());
	printer.generateHyperplanesTable<CIMPUTE_AND_PRINT_POINTS_ORDER,PRINT_SUBGEOMETRIES>(4, geometry4_hyp_table, geometry3_hyp_table.size());

	return EXIT_SUCCESS;
}

int metod_stuff() {

	std::array<std::bitset<48>, 48> lines;
	lines[1 - 1][1 - 1] = lines[1 - 1][2 - 1] = lines[1 - 1][3 - 1] = 1;
	lines[2 - 1][4 - 1] = lines[2 - 1][5 - 1] = lines[2 - 1][6 - 1] = 1;
	lines[3 - 1][7 - 1] = lines[3 - 1][8 - 1] = lines[3 - 1][9 - 1] = 1;
	lines[4 - 1][10 - 1] = lines[4 - 1][11 - 1] = lines[4 - 1][12 - 1] = 1;
	lines[5 - 1][13 - 1] = lines[5 - 1][14 - 1] = lines[5 - 1][15 - 1] = 1;
	lines[6 - 1][16 - 1] = lines[6 - 1][17 - 1] = lines[6 - 1][18 - 1] = 1;
	lines[7 - 1][19 - 1] = lines[7 - 1][20 - 1] = lines[7 - 1][21 - 1] = 1;
	lines[8 - 1][22 - 1] = lines[8 - 1][23 - 1] = lines[8 - 1][24 - 1] = 1;
	lines[9 - 1][25 - 1] = lines[9 - 1][26 - 1] = lines[9 - 1][27 - 1] = 1;
	lines[10 - 1][28 - 1] = lines[10 - 1][29 - 1] = lines[10 - 1][30 - 1] = 1;
	lines[11 - 1][31 - 1] = lines[11 - 1][32 - 1] = lines[11 - 1][33 - 1] = 1;
	lines[12 - 1][34 - 1] = lines[12 - 1][35 - 1] = lines[12 - 1][36 - 1] = 1;
	lines[13 - 1][37 - 1] = lines[13 - 1][38 - 1] = lines[13 - 1][39 - 1] = 1;
	lines[14 - 1][40 - 1] = lines[14 - 1][41 - 1] = lines[14 - 1][42 - 1] = 1;
	lines[15 - 1][43 - 1] = lines[15 - 1][44 - 1] = lines[15 - 1][45 - 1] = 1;
	lines[16 - 1][46 - 1] = lines[16 - 1][47 - 1] = lines[16 - 1][48 - 1] = 1;
	lines[17 - 1][4 - 1] = lines[17 - 1][7 - 1] = lines[17 - 1][10 - 1] = 1;
	lines[18 - 1][1 - 1] = lines[18 - 1][8 - 1] = lines[18 - 1][11 - 1] = 1;
	lines[19 - 1][2 - 1] = lines[19 - 1][5 - 1] = lines[19 - 1][12 - 1] = 1;
	lines[20 - 1][3 - 1] = lines[20 - 1][6 - 1] = lines[20 - 1][9 - 1] = 1;
	lines[21 - 1][37 - 1] = lines[21 - 1][40 - 1] = lines[21 - 1][43 - 1] = 1;
	lines[22 - 1][38 - 1] = lines[22 - 1][41 - 1] = lines[22 - 1][46 - 1] = 1;
	lines[23 - 1][39 - 1] = lines[23 - 1][44 - 1] = lines[23 - 1][47 - 1] = 1;
	lines[24 - 1][42 - 1] = lines[24 - 1][45 - 1] = lines[24 - 1][48 - 1] = 1;
	lines[25 - 1][13 - 1] = lines[25 - 1][25 - 1] = lines[25 - 1][37 - 1] = 1;
	lines[26 - 1][1 - 1] = lines[26 - 1][14 - 1] = lines[26 - 1][38 - 1] = 1;
	lines[27 - 1][2 - 1] = lines[27 - 1][26 - 1] = lines[27 - 1][39 - 1] = 1;
	lines[28 - 1][3 - 1] = lines[28 - 1][15 - 1] = lines[28 - 1][27 - 1] = 1;
	lines[29 - 1][10 - 1] = lines[29 - 1][22 - 1] = lines[29 - 1][34 - 1] = 1;
	lines[30 - 1][11 - 1] = lines[30 - 1][35 - 1] = lines[30 - 1][46 - 1] = 1;
	lines[31 - 1][12 - 1] = lines[31 - 1][23 - 1] = lines[31 - 1][47 - 1] = 1;
	lines[32 - 1][24 - 1] = lines[32 - 1][36 - 1] = lines[32 - 1][48 - 1] = 1;
	lines[33 - 1][13 - 1] = lines[33 - 1][19 - 1] = lines[33 - 1][22 - 1] = 1;
	lines[34 - 1][14 - 1] = lines[34 - 1][16 - 1] = lines[34 - 1][20 - 1] = 1;
	lines[35 - 1][17 - 1] = lines[35 - 1][21 - 1] = lines[35 - 1][23 - 1] = 1;
	lines[36 - 1][15 - 1] = lines[36 - 1][18 - 1] = lines[36 - 1][24 - 1] = 1;
	lines[37 - 1][25 - 1] = lines[37 - 1][28 - 1] = lines[37 - 1][34 - 1] = 1;
	lines[38 - 1][29 - 1] = lines[38 - 1][31 - 1] = lines[38 - 1][35 - 1] = 1;
	lines[39 - 1][26 - 1] = lines[39 - 1][30 - 1] = lines[39 - 1][32 - 1] = 1;
	lines[40 - 1][27 - 1] = lines[40 - 1][33 - 1] = lines[40 - 1][36 - 1] = 1;
	lines[41 - 1][4 - 1] = lines[41 - 1][28 - 1] = lines[41 - 1][40 - 1] = 1;
	lines[42 - 1][16 - 1] = lines[42 - 1][29 - 1] = lines[42 - 1][41 - 1] = 1;
	lines[43 - 1][5 - 1] = lines[43 - 1][17 - 1] = lines[43 - 1][30 - 1] = 1;
	lines[44 - 1][6 - 1] = lines[44 - 1][18 - 1] = lines[44 - 1][42 - 1] = 1;
	lines[45 - 1][7 - 1] = lines[45 - 1][19 - 1] = lines[45 - 1][43 - 1] = 1;
	lines[46 - 1][8 - 1] = lines[46 - 1][20 - 1] = lines[46 - 1][31 - 1] = 1;
	lines[47 - 1][21 - 1] = lines[47 - 1][32 - 1] = lines[47 - 1][44 - 1] = 1;
	lines[48 - 1][9 - 1] = lines[48 - 1][33 - 1] = lines[48 - 1][45 - 1] = 1;

	segre::PointGeometry<2, 3, 48, 48> geometry2(std::move(lines), false);
	auto result = geometry2.findHyperplanesByBruteforce();
}

int main() {
	metod_stuff();
}
