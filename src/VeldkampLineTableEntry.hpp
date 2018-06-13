#ifndef HYPERPLANEFINDER_VELDKAMPLINETABLEENTRY_HPP
#define HYPERPLANEFINDER_VELDKAMPLINETABLEENTRY_HPP


#include <array>
#include <cstddef>
#include <iostream>
#include <vector>
#include <map>

#include "math.hpp"

namespace segre {

	struct VeldkampLineTableEntry {

		VeldkampLineTableEntry()
		  : isProjective(false)
		    , coreNbrPoints(0)
		    , coreNbrLines(0)
		    , pointsType()
		    , count(0) {
		}

		bool operator==(const VeldkampLineTableEntry& entry) const {
			return isProjective == entry.isProjective
			       && coreNbrPoints == entry.coreNbrPoints
			       && coreNbrLines == entry.coreNbrLines
			       && pointsType == entry.pointsType;
		}

		friend std::ostream& operator<<(std::ostream& os, const VeldkampLineTableEntry& entry);

		bool isProjective;
		size_t coreNbrPoints;
		size_t coreNbrLines;
		std::map<long long int, std::size_t> pointsType;
		size_t count;
	};

	template<size_t NbrPointsPerLine>
	struct VeldkampLineTableEntryWithLines {

		VeldkampLineTableEntryWithLines()
		  : entry()
		  , lines() {
		}

		explicit VeldkampLineTableEntryWithLines(const VeldkampLineTableEntry& entry_)
		  : entry(entry_)
		  , lines() {
		}

		VeldkampLineTableEntryWithLines(const VeldkampLineTableEntry& entry_, const std::vector<std::array<unsigned int, NbrPointsPerLine>>& lines_)
		  : entry(entry_)
		  , lines(lines_) {
		}

		VeldkampLineTableEntry entry;
		std::vector<std::array<unsigned int, NbrPointsPerLine>> lines;
	};

	inline std::ostream& operator<<(std::ostream& os, const VeldkampLineTableEntry& entry) {
		os << "VeldkampLineEntry{"
		   << "Proj: " << std::boolalpha << entry.isProjective
		   << ", core{"
		   << "Ps: " << entry.coreNbrPoints
		   << ", Ls: " << entry.coreNbrLines
		   << "}, composition{";

		for(std::map<long long int, std::size_t>::const_iterator iterator = entry.pointsType.cbegin();
		    iterator != entry.pointsType.cend();) {
			os << "H" << iterator->first << ":" << iterator->second;
			if(++iterator != entry.pointsType.cend()) {
				os << ", ";
			}
		}

		os << "}, Crd: " << entry.count << '}';

		return os;
	}
}


#endif //HYPERPLANEFINDER_VELDKAMPLINETABLEENTRY_HPP
