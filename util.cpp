#include "util.h"
#include <sstream>
void get_fields(std::string line, std::vector<std::string> &v)
{
    v.clear();
    std::string s;
    std::istringstream sin(line);
    while(sin >> s)
        v.push_back(s);
}

std::vector<std::string> rsplit(std::string str, std::string delim){
    std::vector<std::string> splits;
    size_t p1 = 0;
    size_t p2 = 0;
    while((p2= str.find(delim,p1)) != std::string::npos){
        splits.push_back(str.substr(p1,p2-p1));
        p1 = p2+delim.size();
    }
    splits.push_back(str.substr(p1));
    return splits;
}
