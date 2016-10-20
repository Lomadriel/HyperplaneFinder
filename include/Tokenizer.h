#ifndef HYPERPLANEFINDERD4_TOKENIZER_H
#define HYPERPLANEFINDERD4_TOKENIZER_H

#include <fstream>
#include <vector>

template <typename T>
void tokenize(const std::basic_string<T>& str,
                                   std::vector<std::basic_string<T>>& tokens,
                                   const std::basic_string<T>& delimiters,
                                   unsigned int limit = 0)
{
    size_t lastPos = str.find_first_not_of(delimiters, 0),
            pos = str.find_first_of(delimiters, lastPos);

    do {
        tokens.push_back(str.substr(lastPos, pos - lastPos));

        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);

    } while ((std::basic_string<T>::npos != pos
              || std::basic_string<T>::npos != lastPos)
             && (!limit || --limit));
}

#endif //HYPERPLANEFINDERD4_TOKENIZER_H
