#ifndef HYPERPLANEFINDER_HYPERPLANETABLEENTRY_HPP
#define HYPERPLANEFINDER_HYPERPLANETABLEENTRY_HPP


#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <set>

namespace segre{

	struct HyperplaneTableEntry {
		HyperplaneTableEntry()
		  : nbrPoints{0}
		  , nbrLines{0}
		  , pointsOfOrder{}
		  , subgeometries{}
		  , count{0}
		  , projective{true} {
		}

		template <typename Map>
		bool map_compare(Map const& lhs, Map const& rhs) const {
			return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
		}

		bool operator==(const HyperplaneTableEntry& oentry) const {
			if (projective != oentry.projective) {
				return false;
			}

			if(nbrPoints != oentry.nbrPoints)
				return false;

			if(nbrLines != oentry.nbrLines)
				return false;

			if(!map_compare(pointsOfOrder, oentry.pointsOfOrder))
				return false;

			// Compare subgeometries (order is not important) with those of the other entry
			// check if all subgeometries correspond to a subgeometry of the other entry
			std::set<const std::map<long long int, size_t>*> already_used;
			for(const std::map<long long int, size_t>& subgeometry : subgeometries){

				// Check if the subgeometry corresponds to a subgeometry of the other entry
				const std::vector<std::map<long long int, size_t>>::const_iterator itt = std::find_if(
				  oentry.subgeometries.cbegin(), oentry.subgeometries.cend(),
				  [&already_used, &subgeometry](const std::map<long long int, size_t>& osubgeometry){
                      // if equal and not already used in a correspondence
                      if(subgeometry == osubgeometry && already_used.count(&osubgeometry) == 0){
                          already_used.insert(&osubgeometry);
                          return true;
                      }
                      return false;
                  });
				if(itt == oentry.subgeometries.cend()){
					return false;
				}
			}

			return true;
		}

		friend std::ostream& operator<<(std::ostream& os, const HyperplaneTableEntry& entry);

		unsigned int nbrPoints;
		unsigned int nbrLines;
		std::map<unsigned int, unsigned int> pointsOfOrder;
		std::vector<std::map<long long int, std::size_t>> subgeometries;
		size_t count;
		bool projective;
	};

	inline std::ostream& operator<<(std::ostream& os, const HyperplaneTableEntry& entry) {
		os << "HyperplaneTableEntry{"
		   << "Ps: " << entry.nbrPoints
		   << ", Ls: " << entry.nbrLines
		   << ", Order={";

		for (auto iterator = entry.pointsOfOrder.cbegin(); iterator != entry.pointsOfOrder.cend();) {
			os << iterator->first << "=" << iterator->second;
			if (++iterator != entry.pointsOfOrder.cend()) {
				os << ", ";
			}
		}

		os << "}, SubGeometry={";

		for (std::size_t i = 0; i < entry.subgeometries.size(); ++i) {
			os << 'D' << i << "={";

			for (auto it = entry.subgeometries[i].cbegin(); it != entry.subgeometries[i].cend();) {
				if (it->first != -1) {
					os << 'H' << it->first << "=" << it->second;
				} else {
					os << "Full=" << it->second;
				}

				if (++it != entry.subgeometries[i].cend()) {
					os << ", ";
				}
			}

			os << "}";

			if (i + 1 != entry.subgeometries.size()) {
				os << ", ";
			}
		}

		/*for (auto iterator = entry.subgeometries.cbegin(); iterator != entry.subgeometries.cend();) {
			if (iterator->first != -1) {
				os << 'H' << iterator->first << "=" << iterator->second;
			} else {
				os << "Full=" << iterator->second;
			}

			if (++iterator != entry.subgeometries.cend()) {
				os << ", ";
			}
		}*/

		os << "}, Crd: " << entry.count << '}';

		return os;
	}
}


#endif //HYPERPLANEFINDER_HYPERPLANETABLEENTRY_HPP
