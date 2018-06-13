#ifndef HYPERPLANEFINDER_PERMUTATIONGENERATOR_HPP
#define HYPERPLANEFINDER_PERMUTATIONGENERATOR_HPP


#include <array>

namespace segre::detail {

	template<size_t nbrElements>
	std::array<unsigned int, nbrElements> generateInitialPermutation() {
		std::array<unsigned int, nbrElements> initial_permutation;
		for(unsigned int i = 0; i < nbrElements; ++i) {
			initial_permutation[i] = i;
		}
		return initial_permutation;
	}

	template<size_t nbrElements>
	std::array<unsigned int, nbrElements - 1> withoutFirst(const std::array<unsigned int, nbrElements>& in_values) {
		std::array<unsigned int, nbrElements - 1> out_values;
		std::copy(in_values.cbegin() + 1, in_values.cend(), out_values.begin());
		return out_values;
	}

	template<size_t permutationSize>
	struct Permutation;

	template<>
	struct Permutation<2> {

		explicit Permutation(const std::array<unsigned int, 2>& values) noexcept
		  : permutation(values)
		    , ended(false) {

		}

		Permutation(const Permutation<2>&) = delete;
		Permutation<2>& operator=(const Permutation<2>&) = delete;

		Permutation(Permutation<2>&&) noexcept = default;
		Permutation<2>& operator=(Permutation<2>&&) noexcept = default;

		void reset(const std::array<unsigned int, 2>& values) {
			permutation = values;
			ended = false;
		}

		bool nextPermutation() {
			if(!ended) {
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
	struct Permutation {

		explicit Permutation(const std::array<unsigned int, permutationSize>& values) noexcept
		  : subPermutation(withoutFirst(values))
		  , permutation(values)
		  , current(values[0])
		  , position(0) {

		}

		Permutation(const Permutation<permutationSize>&) = delete;
		Permutation<permutationSize>& operator=(const Permutation<permutationSize>&) = delete;

		Permutation(Permutation<permutationSize>&&) noexcept = default;
		Permutation<permutationSize>& operator=(Permutation<permutationSize>&&) noexcept = default;

		void reset(const std::array<unsigned int, permutationSize>& values) {
			subPermutation.reset(withoutFirst(values));
			permutation = values;
			current = values[0];
			position = 0;
		}

		bool nextPermutation() {
			++position;
			if(position < permutationSize) {
				std::swap(permutation[position - 1], permutation[position]);
				return true;
			}

			if(!subPermutation.nextPermutation()) {
				return false;
			}

			position = 0;
			permutation[position] = current;
			std::copy(subPermutation.permutation.cbegin(), subPermutation.permutation.cend(), permutation.begin() + 1);
			return true;
		}

		Permutation<permutationSize - 1> subPermutation;
		std::array<unsigned int, permutationSize> permutation;
		unsigned int current;
		unsigned int position;
	};

}

namespace segre {

	template<size_t nbrElements>
	class PermutationGenerator {

	public:

		PermutationGenerator() noexcept
		  : m_current_permutation(detail::generateInitialPermutation<nbrElements>())
		  , m_finished(false) {

		}

		PermutationGenerator(const PermutationGenerator<nbrElements>&) = delete;
		PermutationGenerator<nbrElements>& operator=(const PermutationGenerator<nbrElements>&) = delete;

		PermutationGenerator(PermutationGenerator<nbrElements>&&) noexcept = default;
		PermutationGenerator<nbrElements>& operator=(PermutationGenerator<nbrElements>&&) noexcept = default;

		void reset() {
			m_current_permutation.reset(detail::generateInitialPermutation<nbrElements>());
			m_finished = false;
		}

		std::array<unsigned int, nbrElements> nextPermutation() {
			const std::array<unsigned int, nbrElements> permutation = m_current_permutation.permutation;
			m_finished = !m_current_permutation.nextPermutation();
			return permutation;
		}

		bool isFinished() const {
			return m_finished;
		}

		static constexpr unsigned int getPermutationsNumber(){
			return math::factorial<nbrElements>;
		}

	private:

		detail::Permutation<nbrElements> m_current_permutation;
		bool m_finished;
	};

	template<size_t... nbrElements>
	class MultiPermutationGenerator {

	public:

		MultiPermutationGenerator() noexcept
		  : m_generators{}
		  , m_current_multiperm{}
		  , m_finished{false} {

			apply([](auto& perm, auto& generator) {
				perm = generator.nextPermutation();
			});
		}

		MultiPermutationGenerator(const MultiPermutationGenerator<nbrElements...>&) = delete;
		MultiPermutationGenerator<nbrElements...>& operator=(const MultiPermutationGenerator<nbrElements...>&) = delete;

		MultiPermutationGenerator(MultiPermutationGenerator<nbrElements...>&&) noexcept = default;
		MultiPermutationGenerator<nbrElements...>& operator=(MultiPermutationGenerator<nbrElements...>&&) noexcept = default;

		std::tuple<std::array<unsigned int, nbrElements>...> nextPermutation() {
			const std::tuple<std::array<unsigned int, nbrElements>...> permutations = m_current_multiperm;
			unsigned int i = 0;
			unsigned int to_increment = 0;
			apply([&](auto& perm, auto& generator) {
				if(i == to_increment) {
					if(generator.isFinished()) {
						if(i == sizeof...(nbrElements) - 1) {
							m_finished = true;
						}
						++to_increment;
						generator.reset();
					}
					perm = generator.nextPermutation();
				}
				++i;
			});
			return permutations;
		}

		bool isFinished() const {
			return m_finished;
		}

		static constexpr unsigned int getPermutationsNumber(){
			return (PermutationGenerator<nbrElements>::getPermutationsNumber() * ...);
		}

	private:

		template<typename Func>
		void apply(Func func) {
			apply_impl(func, std::make_index_sequence<sizeof...(nbrElements)>());
		}

		template<typename Func, size_t... Is>
		void apply_impl(Func func, std::index_sequence<Is...>) {
			(func(std::get<Is>(m_current_multiperm), std::get<Is>(m_generators)), ...);
		}

		std::tuple<PermutationGenerator<nbrElements>...> m_generators;
		std::tuple<std::array<unsigned int, nbrElements>...> m_current_multiperm;
		bool m_finished;
	};

}


#endif //HYPERPLANEFINDER_PERMUTATIONGENERATOR_HPP
