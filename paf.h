#ifndef PAF_HEADER
#define PAF_HEADER
#include <map>
#include <string>
#include <cstdint>
#include <vector>
#include <fstream>
#include <sstream>
#include <set>

#include <IITree.h>
#include "locus.h"
class interval{
    public:
        int start;
        int end;
        interval(int start, int end) : start(start), end(end){}
        friend std::ostream& operator<<(std::ostream& os, const interval& i);

        bool operator<(const interval &other) const{
            if(start == other.start){
                return end < other.end;
            }
            return start < other.start;
        }
        bool operator==(const interval &other) const{
            return start == other.start && end == other.end;
        }

        double reciprocal_overlap(const interval &other) const{
            double q = std::min(end,other.end) - std::max(start,other.start);
            double r1 = q/(end-start);
            double r2 = q/(other.end-other.start);
            double result = std::min(r1,r2);
            return std::max(result,0.0);
        }
};

namespace std{
    template <>
        struct hash<interval>{
            std::size_t operator()(const interval& k) const {
                using std::size_t;
                using std::hash;
                // Compute individual hash values for first,
                // second and third and combine them using XOR
                // and bit shifting:
                return ((hash<int>()(k.start)
                            ^ (hash<int>()(k.end) << 1)) >> 1);
            }
        };
}

std::ostream& operator<<(std::ostream& os, const interval& i);

class aligned_segment{
    public:

    std::string chr;
    interval query;
    interval tmplt;
    bool reverse_complemented;
    aligned_segment(std::string chr, interval query, interval tmplt, bool rc):
        chr(chr),
        query(query),
        tmplt(tmplt),
        reverse_complemented(rc){
    }
    aligned_segment(std::string chr, int start1, int end1, int start2, int end2, bool rc):
        chr(chr),
        query(start1,end1),
        tmplt(start2,end2),
        reverse_complemented(rc){}

        friend std::ostream& operator<<(std::ostream& os, const aligned_segment& i);
};

std::ostream& operator<<(std::ostream& os, const aligned_segment& i);


class exon{
    public:
        std::string chr;
        int start;
        int end;
        bool strand;
        std::string gene_id;
        std::string transcript_id;
        std::string exon_id;
        int exon_number;
        std::string gene_name;
        bool coding;
        bool operator<(const exon &other) const { 

            if(chr != other.chr){
                return chr < other.chr;
            }
            if(gene_id != other.gene_id)
                return (gene_id < other.gene_id);
            if(transcript_id != other.transcript_id)
                return (transcript_id < other.transcript_id);
            if(exon_id != other.exon_id)
                return exon_id < other.exon_id;
            if(start != other.start){
                return start < other.start;
            }
            if(end != other.end){
                return end < other.end;
            }
            return strand < other.strand;
        }

        bool operator==(const exon &other) const { 
            return (gene_id == other.gene_id) &&
                start == other.start &&
                end == other.end &&
                transcript_id ==  other.transcript_id &&
                exon_id == other.exon_id &&
                start == other.start &&
                end == other.end &&
                (strand == other.strand);
        }
        operator interval() const{
            return interval(start,end);
        }
        friend std::ostream& operator<<(std::ostream& os, const exon& i);
};

std::ostream& operator<<(std::ostream& os, const exon& i);


namespace std{
    template <>
        struct hash<exon>{
            std::size_t operator()(const exon& k) const{
                using std::size_t;
                using std::hash;
                using std::string;

                // Compute individual hash values for first,
                // second and third and combine them using XOR
                // and bit shifting:

                return ((hash<string>()(k.transcript_id)
                            ^ (hash<string>()(k.gene_id) << 1)) >> 1);
            }
        };
}


class paf_t{
    public:
        const static std::set<std::string> tag_list;

        std::string  qName;
        int32_t qLen;
        int32_t qStart;
        int32_t qEnd;
        char    strand;
        std::string  tName;

        int32_t tLen;
        int32_t tStart;
        int32_t tEnd;
        int32_t bMatch;
        int32_t bAlign;
        int32_t mapq;
        std::map<std::string, std::string> tags;

        //std::vector<aligned_segment> aligned_segments;
        std::vector<std::pair<aligned_segment,exon>> mappings;
        bool operator==(const paf_t &other) const{
            return this->tName  == other.tName;
        }
};

void str_split(std::string str, char delim, std::vector<std::string> &v);

bool compare_paf_query(const paf_t &p1, const paf_t &p2);
bool cosetmpare_paf_target(const paf_t &p1, const paf_t &p2);

class gene_t{
    public:
        std::string gene_name;
        std::string chr;
        int32_t chr_start;
        int32_t chr_end;
        int32_t rev_strand;
};

bool has_overlap(std::string t1, std::string t2, std::map<std::string, gene_t> &gene_info);

template <typename T>
T str2type(std::string str){
    T n;
    std::istringstream sin(str);
    sin >> n;
    return n;
}

template <typename T>
std::string type2str(T v){
    std::ostringstream sout;
    sout << v;
    return sout.str();
}


bool get_next_paf(std::ifstream &fin, paf_t &map, bool is_wg);
bool get_next_paf(std::ifstream &fin, paf_t &map, const IITree<locus,exon> &, bool is_wg);

enum class cigar_character_type{
    matched, onquery, ontemplate, hardclip, softclip, notcigar,
};

cigar_character_type what_is_this_cigar(char c);


void paf2mappings( const std::string &paf, std::vector<std::pair<aligned_segment,exon>> &mappings, int max_skip);
#endif
