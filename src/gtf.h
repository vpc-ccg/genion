#ifndef GTF_H
#define GTF_H
#include <map>
#include <string>
#include <vector>

#include "util.h"

class interval{
    public:
        int start;
        int end;

        interval():start(0),end(0){}
        interval(int start, int end): start(start), end(end){}
        interval(const interval &i1, const interval &i2): start(std::min(i1.start, i2.start)), end(std::max(i1.end, i2.end)){}
        int distance( const interval &other) const{
            if( other.start < this->start){
                return other.distance(*this);
            }
            return other.start - this->end;
        }
        int overlap(const interval &other)  const{
            if (other.end <= start){ //BEFORE
                return 0;
            }
            else if (other.start >= end){ //AFTER
                return 0;
            }
            else if (other.start >= start && other.end <= end){ // IN
                return other.end - other.start;
            }
            else if (other.start < start && other.end > end) { // AROUND
                return end - start;
            }
            else if (other.start < start && other.end < end && other.end > start){ // LEFT OVERLAP
                return other.end - start;
            }
            else if (other.start > start && other.start < end && other.end > end){ //RIGHT OVERLAP
                return end - other.start;
            }
            return 0;
        }
        friend int larger_interval( const interval &i1, const interval &i2){
            return std::max(i1.end - i1.start, i2.end - i2.start);
        }
        double reciprocal_overlap(const interval &other) const {
            return static_cast<double>(overlap(other)) / larger_interval(*this, other);
        }
        bool contains(int pos) const {
            return pos > start && pos < end;
        }
};


class ginterval: public interval{
    public:
        std::string chr;
        bool plus_strand;
        ginterval(): interval(0,0), chr(""), plus_strand(true){}
        ginterval(std::string chr, int start, int end, const std::string &strand): 
            interval(start, end), chr(chr), plus_strand(strand=="+"){}
        ginterval( const ginterval &g1, const ginterval &g2) : interval(g1,g2), chr(g1.chr), plus_strand(g1.plus_strand) {} // Merge constructor
        ginterval( const ginterval &g1) = default;
        virtual ~ginterval() {}
        int overlap( const ginterval &other) const{
            if( chr != other.chr){
                return 0;
            }
            return interval::overlap(other);
        } 
        double reciprocal(const ginterval &other) const {
            return static_cast<double>(overlap(other)) / larger_interval(*this, other);
        }
        bool operator<( const ginterval &other) const{
            if( other.chr != chr){
                return chr < other.chr;
            }
            if( start == other.start){
                return end < other.end;
            }
            return start < other.start;
            //}

        }
        bool operator==(const ginterval &other) const {
            return chr == other.chr && start == other.start && end == other.end && plus_strand == other.plus_strand;
        }
};

struct gtf: public ginterval{
    enum class entry_type{
        gene, transcript, exon,
        five_prime_utr, three_prime_utr,
        start_codon, stop_codon,
        CDS, Selenocysteine,
        other
    };
    static entry_type type_from_string( const std::string& type_str){

        if(type_str == "gene"){
            return entry_type::gene;
        }
        else if(type_str == "transcript"){
            return entry_type::transcript;
        }
        else if(type_str == "exon"){
            return entry_type::exon;
        }
        else if(type_str == "five_prime_utr"){
            return entry_type::five_prime_utr;
        }

        else if(type_str == "three_prime_utr"){
            return entry_type::three_prime_utr;
        }
        else if(type_str == "start_codon"){
            return entry_type::start_codon;
        }
        else if(type_str == "stop_codon"){
            return entry_type::stop_codon;
        }
        else if(type_str == "CDS"){
            return entry_type::CDS;
        }
        else if(type_str == "Selenocysteine"){
            return entry_type::Selenocysteine;
        }
        else{
            return entry_type::other;
        }
    }

    entry_type type;
    std::map<std::string, std::string> info;
    gtf( const std::string &gtf_line){
        std::vector<std::string> fields = rsplit(gtf_line, "\t");
        type = type_from_string(fields[2]);
        chr = fields[0];
        start = stoi(fields[3]);
        end   = stoi(fields[4]);
        plus_strand = (fields[6] == "+");
        
        std::string info_str{fields[8]};
        std::vector<std::string> info_vec = rsplit(info_str, ";");
        strip_for_each(info_vec , " ");

        for( auto iter = info_vec.begin(); iter != info_vec.end(); ++iter){
            if((*iter).size() <= 1){ continue;}
            std::string f{*iter};

            std::vector<std::string> fs = rsplit(f , " ");
            strip_for_each(fs, "\"");
            info[fs[0]] = fs[1];
        }
    }
};

#endif
