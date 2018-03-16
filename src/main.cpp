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
using VLines = std::pair<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>>;

int main() {
	const auto time_start = std::chrono::system_clock::now();

	std::array<std::bitset<PPL>, 1> lines;
	lines[0] = std::bitset<PPL>(math::pow(2UL,PPL) - 1);

	segre::PointGeometry<1, PPL, 1> geometry1(std::move(lines));
	segre::PointGeometry<2, PPL, 8> geometry2(geometry1.computeCartesianProduct(), geometry1.buildTensorPoints());
	segre::PointGeometry<3, PPL, 48> geometry3(geometry2.computeCartesianProduct(), geometry2.buildTensorPoints());
	segre::PointGeometry<4, PPL, 256> geometry4(geometry3.computeCartesianProduct(), geometry3.buildTensorPoints());

	VPoints<2> vPoints2 = geometry2.findHyperplanes(); // brut force
	VLines<2> vLines2 = geometry2.findVeldkampLinesDimension4(vPoints2);
	geometry2.distinguishVeldkampLines(vLines2, vPoints2, geometry3);

	VPoints<3> vPoints3 = geometry2.computeHyperplanes(vPoints2, vLines2.second);
	VLines<3> vLines3 = geometry3.findVeldkampLinesDimension4(vPoints3);
	geometry3.distinguishVeldkampLines(vLines3, vPoints3, geometry4);

	VPoints<4> vPoints4 = geometry3.computeHyperplanes(vPoints3, vLines3.second);
	std::vector<segre::Entry> entries = geometry4.makeTable(vPoints4);
	//std::vector<segre::Entry> entries = geometry3.makeTable(vPoints3);

	std::sort(entries.begin(), entries.end(), [] (const segre::Entry& a,
	                                              const segre::Entry& b) {
		return a.nbrPoints > b.nbrPoints;
	});

	const auto time_end = std::chrono::system_clock::now();
	const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(time_end - time_start);
	std::cout << "Finished in " << static_cast<int>(elapsed.count()) << " seconds\n" << std::endl;

	std::copy(entries.begin(), entries.end(), std::ostream_iterator<segre::Entry>(std::cout, "\n"));

	return EXIT_SUCCESS;
}
