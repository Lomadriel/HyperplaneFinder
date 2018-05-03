#ifndef HYPERPLANEFINDER_INDEX_REPETITION_HPP
#define HYPERPLANEFINDER_INDEX_REPETITION_HPP


// Repetition holder
template <size_t... Ns> struct index_repetition {};

namespace {
	template <size_t rep, size_t... Ns> struct repetition_generator;

	// Recursion
	template <size_t rep, size_t N, size_t... Ns>
	struct repetition_generator<rep, N, Ns...>
	{
		using type = typename repetition_generator<rep - 1, N, N, Ns...>::type;
	};

	// Recursion end
	template <size_t N, size_t... Ns>
	struct repetition_generator<0, N, Ns...>
	{
		using type = index_repetition<Ns...>;
	};
}

// Make repetition
template <size_t rep, size_t N>
using make_index_repetition = typename repetition_generator<rep, N>::type;


#endif //HYPERPLANEFINDER_INDEX_REPETITION_HPP
