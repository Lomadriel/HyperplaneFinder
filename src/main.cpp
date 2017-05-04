#include <iostream>
#include <iterator>
#include <utility>
#include <algorithm>
#include <vector>

#include <PointGeometry.hpp>

int main() {
	std::array<std::bitset<4>, 1> lines;
	lines[0] = std::bitset<4>(15);

	segre::PointGeometry<1, 4, 1> geometry0(std::move(lines));
	segre::PointGeometry<2, 4, 8> geometry(geometry0.computeCartesianProduct(), geometry0.buildTensorPoints());
	segre::PointGeometry<3, 4, 48> geometry1(geometry.computeCartesianProduct(), geometry.buildTensorPoints());
	segre::PointGeometry<4, 4, 256> geometry2(geometry1.computeCartesianProduct(), geometry1.buildTensorPoints());

	auto vPoints0 = geometry.findHyperplanes();
	auto vLines = geometry.findVeldkampLinesDimension4(vPoints0);
	geometry.distinguishVeldkampLines(vLines, vPoints0, geometry1);
	auto vPoints = geometry.computeHyperplanes(vPoints0, vLines.second);

	auto vLines2 = geometry1.findVeldkampLinesDimension4(vPoints);
	geometry1.distinguishVeldkampLines(vLines2, vPoints, geometry2);

	auto vPoints2 = geometry1.computeHyperplanes(vPoints, vLines2.second);
	auto entries = geometry2.makeTable(vPoints2);
	//auto entries = geometry1.makeTable(vPoints);

	std::sort(entries.begin(), entries.end(), [] (const segre::Entry& a,
	                                              const segre::Entry& b) {
		return a.nbrPoints > b.nbrPoints;
	});

	std::copy(entries.begin(), entries.end(), std::ostream_iterator<segre::Entry>(std::cout, "\n"));

	return EXIT_SUCCESS;
}
