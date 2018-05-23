#ifndef HYPERPLANEFINDER_IMPOSSIBLE_HPP
#define HYPERPLANEFINDER_IMPOSSIBLE_HPP


#include <iostream>
#include <cstdlib>

#define IMPOSSIBLE std::cerr << __FILE__ << ":" << __LINE__ << " " << __func__ << ": Supposedly impossible branch taken, program stopped" ; std::exit(EXIT_FAILURE)


#endif //HYPERPLANEFINDER_IMPOSSIBLE_HPP
