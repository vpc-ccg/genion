

#ifndef UTIL_GENION
#define UTIL_GENION
#include <vector>
#include <string>
#include <array>
#include <iostream>
std::vector<std::string> rsplit(std::string str, std::string delim);

void get_fields(std::string line, std::vector<std::string> &v);

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

