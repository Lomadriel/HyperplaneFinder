#ifndef HYPERPLANEFINDER_PERMUTATIONGENERATOR_HPP
#define HYPERPLANEFINDER_PERMUTATIONGENERATOR_HPP


#include <array>

namespace{

	template<size_t nbrElements>
	std::array<unsigned int, nbrElements> generateInitialPermutation(){
		std::array<unsigned int, nbrElements> initial_permutation;
		for(unsigned int i = 0; i < nbrElements; ++i){
			initial_permutation[i] = i;
		}
		return initial_permutation;
	}

	template<size_t nbrElements>
	std::array<unsigned int, nbrElements - 1> withoutFirst(const std::array<unsigned int, nbrElements>& in_values){
		std::array<unsigned int, nbrElements - 1> out_values;
		std::copy(in_values.cbegin() + 1, in_values.cend(), out_values.begin());
		return out_values;
	}

	template<size_t permutationSize>
	struct Permutation;

	template<>
	struct Permutation<2>{

		Permutation(const std::array<unsigned int, 2>& values) noexcept
		  : permutation(values)
		  , ended(false) {

		}

		bool nextPermutation() {
			if(!ended){
				std::swap(permutation[0], permutation[1]);
				ended = true;
				return true;
			}
			return false;
		}

		std::array<unsigned int, 2> permutation;
		bool ended;
	};

	template<size_t permutationSize>
	struct Permutation{

		Permutation(const std::array<unsigned int, permutationSize>& values) noexcept
		  : subPermutation(withoutFirst(values))
		  , permutation(values)
		  , current(values[0])
		  , position(0) {

		}

		bool nextPermutation() {
			++position;
			if(position < permutationSize){
				std::swap(permutation[position - 1], permutation[position]);
				return true;
			}

			if(!subPermutation.nextPermutation()){
				return false;
			}

			position = 0;
			permutation[position] = current;
			std::copy(subPermutation.permutation.cbegin(), subPermutation.permutation.cend(), permutation.begin() + 1);
			return true;
		}

		Permutation<permutationSize - 1> subPermutation;
		std::array<unsigned int, permutationSize> permutation;
		const unsigned int current;
		unsigned int position;
	};

}

namespace segre{

	template<size_t nbrElements>
	class PermutationGenerator {

	public:

		PermutationGenerator() noexcept
		  : m_current_permutation(generateInitialPermutation<nbrElements>())
		  , m_finished(false) {

		}

		void reset() {
			m_current_permutation = generateInitialPermutation<nbrElements>();
			m_finished = false;
		}

		std::array<unsigned int, nbrElements> nextPermutation() {
			std::array<unsigned int, nbrElements> permutation = m_current_permutation.permutation;
			m_finished = !m_current_permutation.nextPermutation();
			return permutation;
		}

		bool isFinished() const {
			return m_finished;
		}

	private:

		Permutation<nbrElements> m_current_permutation;
		bool m_finished;
	};

}


#endif //HYPERPLANEFINDER_PERMUTATIONGENERATOR_HPP
