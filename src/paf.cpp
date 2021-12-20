#include "paf.h"
#include "util.h"

#include <sstream>

#include "cigar.h"

#include <unordered_map>


std::ostream& operator<<(std::ostream& os, const exon& i){
    os << i.gene_id << " " << i.transcript_id << " " << i.exon_number << " " << i.start << " " << i.end;
    return os;
}

std::ostream& operator<<(std::ostream& os, const aligned_segment& i){
    os << i.chr << '\t' << i.query << '\t' << i.tmplt;
    return os;
}
std::ostream& operator<<(std::ostream& os, const interval& i){
    os << i.start << '\t' << i.end;
    return os;
}


cigar_character_type what_is_this_cigar(char c){
    switch(c){
        case 'M':
        case '=':
        case 'X':
            return cigar_character_type::matched;
        case 'D':
        case 'N':
            return cigar_character_type::ontemplate;
        case 'I':
        case 'P':
            return cigar_character_type::onquery;
        case 'H':
            return cigar_character_type::hardclip;
        case 'S':
            return cigar_character_type::softclip;
        default:
            return cigar_character_type::notcigar;
    }
}


//To Merge aligs from the same gene
void paf2mappings( const std::string &paf, std::vector<std::pair<aligned_segment,exon>> &mappings, int max_skip, const IITree<locus, exon> &exon_forest){


    std::istringstream ps(paf);

    std::string id;
    ps >> id;

    int rlen;
    int rstart;
    int rend;
    ps >> rlen;
    ps >> rstart;
    ps >> rend;

    std::string strand;
    ps >> strand;
    bool complemented = strand == "-";
    std::string chr;
    ps >> chr;
    if(chr.find("chr")!= std::string::npos){
        chr = chr.substr(3);
    }

    int template_len;
    int template_start;
    int template_end;
    ps >> template_len;
    ps >> template_start;
    ps >> template_end;

    int num_matches;
    int alig_block_len;
    int mapq;
    ps >> num_matches;
    ps >> alig_block_len;
    ps >> mapq;

    std::string field;
    ps >> field;

    while(!ps.eof() && field.substr(0,2) != "cg"){ ps >> field;}

    if( field.substr(0,2) != "cg"){
        std::cerr << "Paf line doesn't include alignment cigar! Exiting!." << std::endl;
        exit(-1);
    }
    cigar cig(field.substr(5));

    std::vector<aligned_segment> aligs;


    int st = template_start;
    int sq = rstart;
    int et, eq;
    if(complemented){
        sq = rend - 1;
    }

    
    for( auto pair : cig){
        int length = pair.first;
        char c = pair.second;
        cigar_character_type type = what_is_this_cigar(c);

        if( type == cigar_character_type::matched){
            et = st + length;

            if(complemented){
                eq = sq - length;
                aligs.emplace_back(chr,eq,sq-1,st,et-1,strand=="-");
            }
            else{
                eq = sq + length;
                aligs.emplace_back(chr,sq,eq-1,st,et-1,strand=="-");
            }

            sq = eq;
            st = et;
        }
        else if( type == cigar_character_type::ontemplate){
            st = st + length;
        }
        else if( type == cigar_character_type::onquery){
            if(complemented){
                sq = sq - length;
            }   
            else{
                sq = sq + length;
            }
        }
        else{

        }
    }


    std::vector<aligned_segment> merged;
    if( max_skip > 0){

        st = aligs[0].tmplt.start;
        et = aligs[0].tmplt.end;
        sq = aligs[0].query.start;
        eq = aligs[0].query.end;

        for(auto iter = std::begin(aligs); std::next(iter) != std::end(aligs); ++iter){

            if(complemented){

                if( std::next(iter)->tmplt.start - et < max_skip){
                    sq = std::next(iter)->query.start;
                    et = std::next(iter)->tmplt.end;
                
                }
                else{
                    merged.emplace_back(chr,sq,eq,st,et,strand=="-");
                    auto nxt = std::next(iter);
                    st = nxt->tmplt.start;
                    et = nxt->tmplt.end;
                    sq = nxt->query.start;
                    eq = nxt->query.end;
                }
            }
            else{

                if( std::next(iter)->tmplt.start - et < max_skip){
                    et = std::next(iter)->tmplt.end;
                    eq = std::next(iter)->query.end;
                }
                else{
                    merged.emplace_back(chr,sq,eq,st,et,strand=="-");
                    auto nxt = std::next(iter);
                    st = nxt->tmplt.start;
                    et = nxt->tmplt.end;
                    sq = nxt->query.start;
                    eq = nxt->query.end;
                }
            }
        }
        merged.emplace_back(chr,sq,eq,st,et,strand=="-");

    }
    else{
        merged = aligs;
    }

    std::vector<size_t> overlaps;
    for( const auto &alig : merged){
       
        locus start(chr, alig.tmplt.start);
        locus end(chr,   alig.tmplt.end);
        exon_forest.overlap(start,end,overlaps);
        for( const auto &ov : overlaps){

            exon e = exon_forest.data(ov);
            //if( e.strand == alig.reverse_complemented){
                mappings.push_back(std::make_pair(alig,e));
            //}
        }
        overlaps.clear();
    }
}


const std::set<std::string> paf_t::tag_list {
            "tp:A:",
            "cm:i:",
            "s1:i:",
            "s2:i:",
            "NM:i:",
            "MD:Z:",
            "AS:i:",
            "ms:i:",
            "nn:i:",
            "ts:A:",
            "cg:Z:",
            "cs:Z:",
            "dv:f:"
        };
bool compare_paf_query(const paf_t &p1, const paf_t &p2)
{
    return (p1.qStart < p2.qStart);
}

bool compare_paf_target(const paf_t &p1, const paf_t &p2)
{
    if(p1.tName < p2.tName)
    {
        return true;
    }
    else if(p1.tName > p2.tName)
    {
        return false;
    }
    else // equality
    {
        return (p1.tStart < p2.tStart);
    }
}

bool get_next_paf(std::ifstream &fin, paf_t &map, const  IITree<locus, exon> &exon_forest  ,bool is_wg)
{
    std::string line;
    if(getline(fin, line))
    {
        map.mappings.clear();
        if(is_wg){

            paf2mappings(line,  map.mappings, 10, exon_forest);
        }
        std::vector<std::string> fields;
        str_split(line, '\t', fields);
        map.qName = fields[0];
        map.qLen = str2type<int32_t>(fields[1]);
        map.qStart = str2type<int32_t>(fields[2]);
        map.qEnd = str2type<int32_t>(fields[3]);
        map.strand = fields[4][0];
        map.tName = fields[5].substr(0, 15);
        map.tLen = str2type<int32_t>(fields[6]);
        map.tStart = str2type<int32_t>(fields[7]);
        map.tEnd = str2type<int32_t>(fields[8]);
        map.bMatch = str2type<int32_t>(fields[9]);
        map.bAlign = str2type<int32_t>(fields[10]);
        map.mapq = str2type<int32_t>(fields[11]);
        for(unsigned j = 12; j < fields.size(); j++)
        {
            if(paf_t::tag_list.count(fields[j].substr(0, 5)) > 0)
            {
                map.tags[fields[j].substr(0, 2)] = fields[j].substr(5);
            }
        }
        return true;
    }
    else
    {
        return false;
    }
}



bool get_next_paf(std::ifstream &fin, paf_t &map, bool is_wg)
{
    std::string line;
    if(getline(fin, line))
    {
        if(is_wg){
//            paf2mappings(line, map.mappings, 2000);
        }
        std::vector<std::string> fields;
        str_split(line, '\t', fields);
        map.qName = fields[0];
        map.qLen = str2type<int32_t>(fields[1]);
        map.qStart = str2type<int32_t>(fields[2]);
        map.qEnd = str2type<int32_t>(fields[3]);
        map.strand = fields[4][0];
        map.tName = fields[5].substr(0, 15);
        map.tLen = str2type<int32_t>(fields[6]);
        map.tStart = str2type<int32_t>(fields[7]);
        map.tEnd = str2type<int32_t>(fields[8]);
        map.bMatch = str2type<int32_t>(fields[9]);
        map.bAlign = str2type<int32_t>(fields[10]);
        map.mapq = str2type<int32_t>(fields[11]);
        for(unsigned j = 12; j < fields.size(); j++)
        {
            if(paf_t::tag_list.count(fields[j].substr(0, 5)) > 0)
            {
                map.tags[fields[j].substr(0, 2)] = fields[j].substr(5);
            }
            // else it is comment
        }
        return true;
    }
    else
    {
        return false;
    }
}

void str_split(std::string str, char delim, std::vector<std::string> &v)
{
    v.clear();
    size_t p1 = 0;
    size_t p2 = 0;
    while((p2 = str.find(delim, p1)) != std::string::npos)
    {
        v.push_back(str.substr(p1, p2-p1));
        p1 = p2+1;
    }
    v.push_back(str.substr(p1));
}



// TODO: Make overlap less stringent
bool has_overlap(std::string t1, std::string t2, std::map<std::string, gene_t> &gene_info)
{
    if(gene_info[t1].chr != gene_info[t2].chr)
        return false;
    // return true; // only genes on different chromosoms do not overlap!!
    int32_t vicinity = 1000; // TODO: this might be removing some good candidates that are caused by short deletions
    if(gene_info[t1].chr_start < gene_info[t2].chr_start)
    {
        if(gene_info[t1].chr_end + vicinity < gene_info[t2].chr_start)
            return false;
        else
            return true;
    }
    else // gene_info[t1].chr_start >= gene_info[t2].chr_start
    {
        if(gene_info[t2].chr_end + vicinity < gene_info[t1].chr_start)
            return false;
        else
            return true;
    }
}


