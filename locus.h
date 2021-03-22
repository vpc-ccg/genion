#ifndef LOCUS_HEADER
#define LOCUS_HEADER

#include <string>


struct locus{
    std::string chr;
    int position;
    
    locus( const locus &other) : chr(other.chr), position(other.position){}
    locus( std::string  chr, int position) : chr(chr), position(position){}
    locus():chr(""),position(0){}

    bool operator>(const locus &other) const{
        if(chr > other.chr){
            return true;
        }else if (chr < other.chr){
            return false;
        }
        return position > other.position;
    }
    bool operator<(const locus &other) const{
        if(chr < other.chr){
            return true;
        }else if (chr > other.chr){
            return false;
        }
        return position < other.position;
    }
    bool operator==(const locus &other) const{
        return chr == other.chr && position == other.position;
    }
};



#endif
