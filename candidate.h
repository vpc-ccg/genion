#ifndef CANDIDATE_HEADER
#define CANDIDATE_HEADER
#include <algorithm>
#include "paf.h"
#include <vector>
#include "chainchain.h"
#include <algorithm>

class Candidate{
    public:
    class Alignment{
        private:
            std::vector<paf_t> primary_list;
            std::vector<paf_t> secondary_list;
        public:
            Alignment(){}
            Alignment(std::vector<paf_t> map_list) {
                for(paf_t paf : map_list){
                    if(paf.tags["tp"] == "P")
                    {
                        primary_list.push_back(paf);
                    }else{
                        secondary_list.push_back(paf);
                    }
                }
                std::sort(primary_list.begin(), primary_list.end(), compare_paf_query);
                std::sort(secondary_list.begin(), secondary_list.end(), compare_paf_query);
            }
            paf_t prefix() const{
                return primary_list.front();
            }
            paf_t suffix() const{
//                return primary_list.back();
                  return  primary_list[std::min(static_cast<int>(primary_list.size())-1,1)];
            }
            unsigned primary_count() const{
                return primary_list.size();
            }
            const auto csecondary_begin() const {
               return secondary_list.cbegin();
            }
            const auto csecondary_end() const {
               return secondary_list.cend();
            }
            const auto cprimary_begin() const {
                return primary_list.cbegin();
            }
            const auto cprimary_end() const {
                return primary_list.cend();
            }
    };
    Alignment transcriptome;
    Alignment whole_genome;
    std::vector<std::pair<aligned_segment, exon>> canonical;
    void print_chains(std::ofstream *chainer){
        for( auto pp : canonical){
        *chainer << "\t" << pp.first.tmplt.start << "\t" << pp.first.tmplt.end << "\t" <<
                        pp.first.chr << "\t" <<  pp.first.query.start << "\t" << pp.first.query.end << "\t" << 
                        pp.first.reverse_complemented << "\t" << 
                        pp.second.chr << "\t" << pp.second.start << "\t" << pp.second.end << "\t" << pp.second.strand  << "\t" <<
                        pp.second.gene_id << "\t" << pp.second.transcript_id << "\t" << pp.second.exon_number << "\n";
        }
    }
    bool easy;
    
    bool filtered;
    std::vector<bool> filter_stats;

    Candidate() : filtered(false){}
    Candidate(std::vector<paf_t> map_list): transcriptome(Alignment(map_list)), filtered(false){}
    Candidate(const Alignment &tra) : transcriptome(tra), filtered(false){}
    Candidate(const Alignment &tra, const Alignment &wg): transcriptome(tra), whole_genome(wg), filtered(false){}
    void set_wg(std::vector<paf_t> wg_map_list, bool add_supp_reads = false){
        whole_genome = wg_map_list;
        std::vector<std::pair<aligned_segment,exon>> mappings;
        for(auto iter =  whole_genome.cprimary_begin();
                 iter != whole_genome.cprimary_end();
               ++iter){
//            for( auto &mappa : iter->mappings){
//                std::cerr   << static_cast<int>(iter-whole_genome.cprimary_begin()) << "\t" << mappa.first << "\t-\t" << mappa.second << "\n";
//            }
            mappings.insert(mappings.begin(),iter->mappings.begin(),iter->mappings.end());
    //        aligs.insert(aligs.begin(),iter->alignments.begin(),iter->alignments.end());
        }
        if(add_supp_reads){
            for(auto iter = whole_genome.csecondary_begin();
                    iter != whole_genome.csecondary_end();
                    ++iter){
                mappings.insert(mappings.begin(),iter->mappings.begin(),iter->mappings.end());
            }
        }
        chainer c(mappings,25);

        std::vector<std::pair<aligned_segment, exon>> best_chain = c.chain();
        canonical = best_chain;
    }
    void add_filter_status(bool status){
        filtered = filtered || !status;
        filter_stats.push_back(status);
    }

    void print_filter_status_wg( std::ofstream &feature_table, const std::string &id){
        feature_table << id << "\t";
        std::unordered_set<std::string> genes;

        if( canonical.size() == 0){

            feature_table << "INTERGENIC::INTERGENIC";
            for(bool val: filter_stats){
                feature_table << "\t" <<  val;
            }
            feature_table << "\n";
            return;
        }
        for(auto alig : canonical){
            genes.insert( alig.second.gene_id);
        }
        auto g1 = std::begin(genes);
        auto g2 = std::next(g1);
        if(g2==genes.end()){
            g2 = g1;
        }
        feature_table << *g1 << "::" << *g2;
        for(bool val: filter_stats){
            feature_table << "\t" <<  val;
        }
        feature_table << "\n";
    }
    void print_filter_status( std::ofstream &feature_table, const std::string &id){
        feature_table << id << "\t";

        feature_table << this->transcriptome.prefix().tName << "::";
        feature_table << this->transcriptome.suffix().tName;
        for(bool val: filter_stats){
            feature_table << "\t" <<  val;
        }
        feature_table << "\n";
    }
};
#endif
