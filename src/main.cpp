#include <iostream>
#include <iterator>
#include <utility>
#include <algorithm>
#include <vector>
#include <chrono>

#include <PointGeometry.hpp>

constexpr size_t PPL = 4; // Points Per Lines

template<int N>
using VPoints = std::vector<std::bitset<math::pow(PPL,N)>>;

template<int>
using VLines = segre::VeldkampLines<PPL>;

int main() {
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

	VPoints<3> vPoints3 = geometry2.computeHyperplanesFromVeldkampLines(vPoints2, vLines2.projectives);
	VLines<3> vLines3 = geometry3.computeVeldkampLines(vPoints3);
	geometry3.distinguishVeldkampLines(vLines3, vPoints3, geometry4);
	std::vector<segre::HyperplaneTableEntry> geometry3_table = geometry3.makeHyperplaneTable<false>(vPoints3);

	std::sort(geometry3_table.begin(), geometry3_table.end(), [] (const segre::HyperplaneTableEntry& a, const segre::HyperplaneTableEntry& b) {
		return a.nbrPoints > b.nbrPoints;
	});

	VPoints<4> vPoints4 = geometry3.computeHyperplanesFromVeldkampLines(vPoints3, vLines3.projectives);
	std::vector<segre::HyperplaneTableEntry> entries = geometry4.makeHyperplaneTable<false>(vPoints4, geometry3_table);

	std::sort(entries.begin(), entries.end(), [] (const segre::HyperplaneTableEntry& a,
	                                              const segre::HyperplaneTableEntry& b) {
		return a.nbrPoints > b.nbrPoints;
	});

	const auto time_end = std::chrono::system_clock::now();
	const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	std::cout << "Finished in " << static_cast<int>(elapsed.count()) << " seconds\n" << std::endl;

	std::cout << "Table 3\n\n";
	std::copy(geometry3_table.begin(), geometry3_table.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));

	std::cout << "Table 4\n";
	std::copy(entries.begin(), entries.end(), std::ostream_iterator<segre::HyperplaneTableEntry>(std::cout, "\n"));

	return EXIT_SUCCESS;
}
