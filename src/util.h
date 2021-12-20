

#ifndef UTIL_GENION
#define UTIL_GENION
#include <vector>
#include <string>
#include <array>
#include <iostream>
std::vector<std::string> rsplit(std::string str, std::string delim);

void get_fields(std::string line, std::vector<std::string> &v);

inline std::string strip_str(const std::string &inpt, const std::string &chrs)
{

    auto frst = inpt.find_first_not_of(chrs);
    auto last = inpt.find_last_not_of(chrs);
    if(frst == std::string::npos || last == std::string::npos){
        return "";
    }
    return inpt.substr(frst, last-frst+1);
}
inline void strip_for_each(std::vector<std::string> &vec, const std::string &chrs = " "){

    for( std::string &st : vec){
        st = strip_str(st, chrs);
    }
}

namespace std{
    template <typename T, size_t N>
    std::ostream &operator<< (std::ostream &os, const std::array<T,N> arr){
        
        for(const T& t: arr){
            os << t << "\t";
        }
        return os;
    }

    template <typename T>
    std::ostream &operator<< (std::ostream &os, const std::vector<T> arr){
        
        for(const T& t: arr){
            os << t << "\t";
        }
        return os;
    }

}
#endif

