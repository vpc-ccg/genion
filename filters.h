#ifndef FILTERS_H
#define FILTERS_H
#include <string>

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <iostream>
#include <zlib.h>

#include "paf.h"
#include "candidate.h"
#include "chainchain.h"
#include "locus.h"

#include <algorithm>
#include <IITree.h>


std::string get_prefix_suffix(const Candidate &cand);
//Hash function injection
namespace std{
    template<> struct hash<pair<string,string>>
    {
        typedef pair<string,string> argument_type;
        typedef std::size_t result_type;
        result_type operator()(argument_type const& s) const noexcept
        {
            result_type const h1 ( std::hash<std::string>{}(s.first) );
            result_type const h2 ( std::hash<std::string>{}(s.second) );
            return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion)
        }
    };
}



class AbstractFilter{
    public:
    virtual ~AbstractFilter() = 0;
    virtual bool operator()(const Candidate &cand) const = 0;

};


        
class StrandSwitchFilter : public AbstractFilter{
    std::map<std::string, gene_t> gene_info;

    public:
    
    StrandSwitchFilter(const std::map<std::string, gene_t> &gene_info):
        gene_info(gene_info) {}
    virtual bool operator()(const Candidate &cand) const;
    ~StrandSwitchFilter(){
    }
};


class WholeGenomeSelfAlignFilter: public AbstractFilter{
    std::unordered_set<std::pair<std::string,std::string>> homologs_g2g;
    bool homology_check(const std::string &prefix, const std::string &suffix) const;
    public:
    

    WholeGenomeSelfAlignFilter(std::unordered_set<std::pair<std::string,std::string>> homologs_g2g): homologs_g2g(homologs_g2g){}
    virtual bool operator()(const Candidate &cand) const;
    ~WholeGenomeSelfAlignFilter(){

    }
};
class TranscriptomeReferenceSelfAlignFilter: public AbstractFilter{
    std::unordered_set<std::pair<std::string,std::string>> homologs_g2g;
    bool homology_check(const std::string &prefix, const std::string &suffix) const;
    public:
    

    TranscriptomeReferenceSelfAlignFilter(std::unordered_set<std::pair<std::string,std::string>> homologs_g2g): homologs_g2g(homologs_g2g){}
    virtual bool operator()(const Candidate &cand) const;
    ~TranscriptomeReferenceSelfAlignFilter(){

    }
};
class HomologyFilter: public AbstractFilter{
    std::unordered_map<std::pair<std::string,std::string>,double> homolog2similarity;

    double min_idt;
    double min_end;
    
    double max_sim;
    bool homology_check(const std::string &prefix, const std::string &suffix) const;
    public:
    
    HomologyFilter(std::unordered_map<std::pair<std::string,std::string>,double> homolog2similarity, double min_idt, double min_end): homolog2similarity(homolog2similarity),min_idt(min_idt),min_end(min_end),max_sim(0){ }
    HomologyFilter(std::unordered_map<std::pair<std::string,std::string>,double> homolog2similarity, double min_idt, double min_end, double max_sim): homolog2similarity(homolog2similarity),min_idt(min_idt),min_end(min_end),max_sim(max_sim){}
    virtual bool operator()(const Candidate &cand) const;
    ~HomologyFilter(){

    }
};

class WholeGenomePrimaryAlignmentCountFilter : public AbstractFilter{
    unsigned min_count;
    public:
    WholeGenomePrimaryAlignmentCountFilter(int min_count) : min_count(min_count){
    }
    virtual bool operator()(const Candidate &cand) const override;
    ~WholeGenomePrimaryAlignmentCountFilter(){

    }
};

class TranscriptomePrimaryAlignmentCountFilter : public AbstractFilter{
    unsigned min_count;
    public:
    TranscriptomePrimaryAlignmentCountFilter(int min_count) : min_count(min_count){
    }
    virtual bool operator()(const Candidate &cand) const override;
    ~TranscriptomePrimaryAlignmentCountFilter(){

    }
};

class TranscriptomeProperPrefixFilter : public AbstractFilter{
    double min_end;
    public:
    TranscriptomeProperPrefixFilter(double min_end): min_end(min_end){
    }
    virtual bool operator()(const Candidate &cand) const override;
    ~TranscriptomeProperPrefixFilter(){

    }
};

class TranscriptomeProperSuffixFilter : public AbstractFilter{
    double min_end;

    public:
    TranscriptomeProperSuffixFilter(double min_end): min_end(min_end){
    }
    virtual bool operator()(const Candidate &cand) const override;
    ~TranscriptomeProperSuffixFilter(){

    }

};

class TranscriptomeStrictMidFilter : public AbstractFilter{

    int max_mid_dist;
    public:
    TranscriptomeStrictMidFilter(double max_distance_between_alignments): max_mid_dist(max_distance_between_alignments){
    }

    virtual bool operator()(const Candidate &cand) const override;

    ~TranscriptomeStrictMidFilter(){

    }
};

class TranscriptomeProperMidFilter : public AbstractFilter{

    double max_mid;
    public:
    TranscriptomeProperMidFilter(double max_mid): max_mid(max_mid){
    }

    virtual bool operator()(const Candidate &cand) const override;

    ~TranscriptomeProperMidFilter(){

    }
};

class TranscriptomePalindromeFilter : public AbstractFilter{
    std::map<std::string, gene_t> &gene_info;
    public:
    TranscriptomePalindromeFilter(std::map<std::string, gene_t> &gene_info):gene_info(gene_info){
    }

    virtual bool operator()(const Candidate &cand) const override;

    ~TranscriptomePalindromeFilter(){

    }
};


class TranscriptomeSecondaryPalindromeFilter : public AbstractFilter{

    double min_idt; // minimum identity of a "good" alignment
    double min_end;

    public:
    TranscriptomeSecondaryPalindromeFilter(double min_idt, double min_end) : min_idt(min_idt), min_end(min_end){
    }

    virtual bool operator()(const Candidate &cand) const override;
    ~TranscriptomeSecondaryPalindromeFilter(){
    }

};
//ChainFilterMinSegment::operator
class ChainFilterMinSegment : public AbstractFilter{
    IITree<locus, exon> exon_forest;
    int min_base_count;
    public:

    ChainFilterMinSegment ( const std::string &gtf_path, int min_base_count): exon_forest(build_exon_forest(gtf_path)),
         min_base_count(min_base_count){

    }
    ChainFilterMinSegment ( const IITree<locus, exon> &exon_forest,int min_base_count): exon_forest(exon_forest),
         min_base_count(min_base_count){

    }
    virtual bool operator()(const Candidate &cand) const override;
    ~ChainFilterMinSegment(){
    }
};
class ChainFilter : public AbstractFilter{
    IITree<locus, exon> exon_forest;
    public:

    ChainFilter ( const std::string &gtf_path): exon_forest(build_exon_forest(gtf_path)){

    }
    ChainFilter ( const IITree<locus, exon> &exon_forest): exon_forest(exon_forest){

    }
    virtual bool operator()(const Candidate &cand) const override;
    ~ChainFilter(){
    }
};

class QaDFilter : public AbstractFilter{

    public:
    QaDFilter (){
    }

    virtual bool operator()(const Candidate &cand) const override;
    ~QaDFilter(){
    }
};

class WGaTConsistentFilter : public AbstractFilter{

    public:
    WGaTConsistentFilter (){
    }

    virtual bool operator()(const Candidate &cand) const override;
    ~WGaTConsistentFilter(){

    }
};

#endif

