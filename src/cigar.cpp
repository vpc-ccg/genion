#include <string>
#include <iostream>
#include <vector>

#include "cigar.h"

using std::string;
using std::vector;

#define CIGAR_CHARACTERS "DHIMNPSX="

cigar::cigar(const string &cgr){
    size_t prev = -1;
    size_t index = cgr.find_first_of(CIGAR_CHARACTERS);
    while(index != string::npos){
        int len = std::stoi(cgr.substr(1+prev,index));
        cigarray.push_back(std::make_pair(len,cgr[index]));
        prev = index;
        index = cgr.find_first_of(CIGAR_CHARACTERS,index + 1);
    }
}
decltype(cigar::cigarray.begin()) cigar::begin(){
    return cigarray.begin();
}

decltype(cigar::cigarray.end()) cigar::end(){
    return cigarray.end();
}

auto cigar::operator [](size_t index) const{
    return cigarray[index];
}
