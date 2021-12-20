#include "filters.h"
#include "paf.h"
#include <array>

AbstractFilter::~AbstractFilter(){
}

std::string get_prefix_suffix(const Candidate &cand){
    return cand.transcriptome.prefix().tName + "\t" + cand.transcriptome.suffix().tName;
}

bool HomologyFilter::homology_check(const std::string &prefix, const std::string &suffix) const{

    auto similar = homolog2similarity.find(std::make_pair(prefix,suffix));
    if( similar != homolog2similarity.end()){
        if( similar->second > max_sim){
            return false;
        }
    }
    similar = homolog2similarity.find(std::make_pair(suffix,prefix));
    if( similar != homolog2similarity.end()){
        if( similar->second > max_sim){
            return false;
        }
    }
    return true;
}

bool StrandSwitchFilter::operator()(const Candidate &cand) const{

    paf_t pfix = cand.transcriptome.prefix();
    paf_t sfix = cand.transcriptome.suffix();
    if(pfix == sfix){return false;}
    auto prefix_gene_p = gene_info.find(pfix.tName);
    auto suffix_gene_p = gene_info.find(sfix.tName);

    if(prefix_gene_p == gene_info.end()){
        std::cerr << "Cannot find prefix gene " << pfix.tName << " in gene info" << std::endl;
        exit(-1);
    }
    if(suffix_gene_p == gene_info.end()){
        std::cerr << "Cannot find suffix gene " << sfix.tName << " in gene info" << std::endl;
        exit(-1);
    }

//    int prefix_strand = prefix_gene_p->second.rev_strand;
//    int suffix_strand = suffix_gene_p->second.rev_strand;

    return pfix.strand == sfix.strand; 
}

bool TranscriptomeReferenceSelfAlignFilter::operator()(const Candidate &cand) const{  
    paf_t pfix = cand.transcriptome.prefix();
    paf_t sfix = cand.transcriptome.suffix();
    if(pfix == sfix){return false;}
    auto similar = homologs_g2g.find(std::make_pair(cand.transcriptome.prefix().tName,cand.transcriptome.suffix().tName));
    if(similar != homologs_g2g.end()){
        return false;
    }
    return true;
}



bool HomologyFilter::operator()(const Candidate &cand) const{  
    paf_t pfix = cand.transcriptome.prefix();
    paf_t sfix = cand.transcriptome.suffix();
    if(pfix == sfix){return false;}
    if( !homology_check(cand.transcriptome.prefix().tName,cand.transcriptome.suffix().tName)){
        return false;
    }
    for(auto iter = cand.transcriptome.csecondary_begin(); iter != cand.transcriptome.csecondary_end(); iter++){
    
        if(cand.transcriptome.prefix().qStart < iter->qStart){ // it is a suffix
            if(((double)(iter->qEnd - iter->qStart) / (iter->qLen - iter->qStart)) >= min_end){ // it is large
                if(((double)iter->bMatch / iter->bAlign) >= min_idt){ // it is good
                    if(!homology_check(iter->tName,cand.transcriptome.prefix().tName)){
                        return false;
                    }        
                }
            }
        }else if(iter->qStart < cand.transcriptome.suffix().qStart){ // it is a suffix
            if(((double)(iter->qEnd - iter->qStart) / iter->qEnd) >= min_end){ // it is large
                if(((double)iter->bMatch / iter->bAlign) >= min_idt){ // it is good
                    if(!homology_check(iter->tName,cand.transcriptome.suffix().tName)){
                        return false;
                    }
                }
            }
        }

    }
    return true;
}


bool WholeGenomePrimaryAlignmentCountFilter::operator()(const Candidate &cand) const {
    paf_t pfix = cand.whole_genome.prefix();
    paf_t sfix = cand.whole_genome.suffix();
    if(pfix == sfix){return false;}
    return cand.whole_genome.primary_count()>=min_count;
}
bool TranscriptomePrimaryAlignmentCountFilter::operator()(const Candidate &cand) const {
    paf_t pfix = cand.transcriptome.prefix();
    paf_t sfix = cand.transcriptome.suffix();
    if(pfix == sfix){return false;}
    return cand.transcriptome.primary_count()>=min_count;
}

bool TranscriptomeProperPrefixFilter::operator()(const Candidate &cand) const {
    paf_t pfix = cand.transcriptome.prefix();
    paf_t sfix = cand.transcriptome.suffix();
    if(pfix == sfix){return false;}
    if(((double)(cand.transcriptome.prefix().qEnd - cand.transcriptome.prefix().qStart) / cand.transcriptome.prefix().qEnd) < min_end){
        return false;
    }
    return true;
}

bool TranscriptomeProperSuffixFilter::operator()(const Candidate &cand) const {
    paf_t pfix = cand.transcriptome.prefix();
    paf_t sfix = cand.transcriptome.suffix();
    if(pfix == sfix){return false;}
    if(((double)(cand.transcriptome.suffix().qEnd - cand.transcriptome.suffix().qStart) / (cand.transcriptome.suffix().qLen - cand.transcriptome.suffix().qStart)) < min_end){
        return false;
    }
    return true;
} 

bool TranscriptomeStrictMidFilter::operator()(const Candidate &cand) const {

    auto pfix = cand.transcriptome.prefix();
    auto sfix = cand.transcriptome.suffix();
    if(pfix == sfix){return false;}

    std::array<int,2> prefix_indices{{pfix.qStart, pfix.qEnd}};
    std::array<int,2> suffix_indices{{sfix.qStart, sfix.qEnd}};
    if( prefix_indices[0] > prefix_indices[1]){
        std::swap(prefix_indices[0],prefix_indices[1]);
    }
    if( suffix_indices[0] > suffix_indices[1]){
        std::swap(suffix_indices[0],suffix_indices[1]);
    }
    return std::abs(suffix_indices[0] - prefix_indices[1]) < max_mid_dist;
}

bool TranscriptomeProperMidFilter::operator()(const Candidate &cand) const {
    if(((double)(cand.transcriptome.suffix().qStart - cand.transcriptome.prefix().qEnd) / cand.transcriptome.suffix().qLen) > max_mid){
        return false;
    }
    return true;
}

bool TranscriptomePalindromeFilter::operator()(const Candidate &cand) const {
    auto pfix = cand.transcriptome.prefix();
    auto sfix = cand.transcriptome.suffix();
    if(pfix == sfix){return false;}
    if(cand.transcriptome.prefix().tName == cand.transcriptome.suffix().tName){
        return false;
    } 
    if(has_overlap(cand.transcriptome.prefix().tName, cand.transcriptome.suffix().tName, gene_info) == true){
        return false;
    }
    return true;
}

bool TranscriptomeSecondaryPalindromeFilter::operator()(const Candidate &cand) const {

    auto pfix = cand.transcriptome.prefix();
    auto sfix = cand.transcriptome.suffix();
    if(pfix == sfix){return false;}
    // is there a secondary suffix for our primary prefix?
    for( auto iter = cand.transcriptome.csecondary_begin(); iter != cand.transcriptome.csecondary_end(); iter++){
        if(cand.transcriptome.prefix().qStart < iter->qStart){ // it is a suffix
            if(((double)(iter->qEnd - iter->qStart) / (iter->qLen - iter->qStart)) >= min_end){ // it is large
                if(((double)iter->bMatch / iter->bAlign) >= min_idt){ // it is good
                    if(cand.transcriptome.prefix().tName == iter->tName){ // same gene
                        return false;
                    }
                }
            }
        }
    }

    // is there a secondary prefix for our primary suffix?
    for( auto iter = cand.transcriptome.csecondary_begin(); iter != cand.transcriptome.csecondary_end(); iter++){
        if(iter->qStart < cand.transcriptome.suffix().qStart){ // it is a suffix
            if(((double)(iter->qEnd - iter->qStart) / iter->qEnd) >= min_end){ // it is large
                if(((double)iter->bMatch / iter->bAlign) >= min_idt){ // it is good
                    if(cand.transcriptome.suffix().tName == iter->tName){ // same gene
                        return false;
                    }
                }
            }
        }
    }
    return true;
}



bool ChainFilter::operator()(const Candidate &cand) const{
//    if(cand.easy){
//        return false;
//    }
    std::unordered_set<std::string> genes;
    for(auto alig : cand.canonical){
        genes.insert( alig.second.gene_id);
    }
    return genes.size() > 1;
}


bool ChainFilterMinSegment::operator()(const Candidate &cand) const{
//    if(cand.easy){
//        return false;
//    }
    std::unordered_map<std::string,int> genes;
    for(auto alig : cand.canonical){
        //genes.insert( alig.second.gene_id);
        genes[alig.second.gene_id] += (alig.first.tmplt.end - alig.first.tmplt.start);
    }
    for(auto iter = genes.begin(); iter != genes.end(); ) {
        if (iter->second < min_base_count) {
              iter = genes.erase(iter);
        }
        else {
              ++iter;
        }
    }
    return genes.size() > 1;
}

bool WholeGenomeSelfAlignFilter::operator()(const Candidate &cand) const{  
    std::unordered_set<std::string> genes;
    for(auto alig : cand.canonical){
        genes.insert( alig.second.gene_id);
    }
    auto g1 = std::begin(genes);
    auto g2 = std::next(g1);
    auto similar = homologs_g2g.find(std::make_pair(*g1,*g2));
    if(similar != homologs_g2g.end()){
        return false;
    }
    return true;
}

/*
 * TODO
 * Checks if breakpoints for whole genome and transcriptome alignment are consisitent
 *
 */
bool WGaTConsistentFilter::operator()(const Candidate &cand) const{ 
    for( auto iter = cand.canonical.begin();
            std::next(iter) !=cand.canonical.end();
              ++iter){
        
    }
   return false; 
}

