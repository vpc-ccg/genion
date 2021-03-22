#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <utility>

#include <cmath>

#include "util.h"
#include "cigar.h"
#include "chainchain.h"

#include <IITree.h>
using std::vector;
using std::set;
using std::pair;
using std::unordered_set;
using std::unordered_map;
using std::ostream;
using std::ifstream;


//9087.180525.9087        3605    75      3605    +       7       159345973       76036108        76280132        3271    3862    60      NM:i:591        ms:i:2081       AS:i:2045       nn:i:0  ts:A:+  tp:A:P  cm:i:532        s1:i:1940       s2:i:174        de:f:0.1026     rl:i:701        cg:Z:6M1I24M1D9M1I12M1I1M4D78M2I3M3D12M2D2M2D18M1D2M1I33M1D12M4D11M3D23M2I23M2D3M1D3M2I2M1I1M2I2M2D2M240272N18M1D7M2I13M3D1M1D3M3D26M1D4M1D6M1I1M1I26M2D3M1D7M1D9M3I10M3D40M3I7M1I8M1D1M1D13M1D17M1I24M1D3M2D17M2D13M3I7M2I42M2I8M1D13M1D7M3D4M1D11M2I31M1I16M2D11M2D8M1I13M4D28M2I10M1D22M3D33M2D25M1I12M1I4M1I13M1D40M3D26M3I40M2D1M2D3M5D32M1D7M2I27M2D17M3D38M6D6M1D9M1D44M1I7M7D6M1D6M2D8M1I15M1D18M10D19M1D15M4D25M1I2M1D10M1D9M1D11M2D3M1D25M3D6M2D12M1I15M1D3M1D42M2D10M1D29M1D32M5D10M5D22M6D7M2I18M1D10M2D3M1D2M2D19M1D13M3D2M1I12M3D5M1I16M1D6M1D9M2D13M13D12M2D9M1D3M2I46M2D28M1D14M1I18M1D2M1D5M1D2M1I19M1D8M1I7M3D1M1D5M1I25M3D23M1I13M1I9M2D10M1I7M1D9M1D14M1D8M2I6M1I26M1I14M14D6M2D59M2D5M1I9M2I3M1D20M1I15M5D6M5D4M1D1M1D15M4D41M2I14M4D4M2D21M4I8M1D14M3D8M1D32M1I29M2D39M3D19M1I10M1D3M1D9M1D16M2D4M1I19M1D4M2D55M1I32M2I6M5D4M1D20M1D3M1I9M1D5M3I26M1D41M1D27M4D41M2I25M1I31M1D12M1D7M4I10M1I11M1D19M1I4M3D26M1D4M1D12M1D4M1I6M2D4M1D6M4D34M1D32M1D15M5D31M1D6M1I2M1I5M1I34M3D20M2D16M2D8M1D18M2I27M1D15M1D9M1D34M3I3M2I66M5D14M1D32M5D29M1D39M2I1M1I4M1D6M2D28M   ENSG00000117385 ENSG00000135679


void cigar2intervals( int alignment_start, cigar cig, int max_na, vector<interval> &result){
    int start = alignment_start;
    int end;
    vector<interval> ivals;

    for( auto pair : cig){
        int len = pair.first;
        char  c = pair.second;
        cigar_character_type type = what_is_this_cigar(c);
        if( type == cigar_character_type::matched){
            end = start + len;
            ivals.emplace_back(start,end-1);
            start = end;
        }
        else if( type == cigar_character_type::ontemplate){
            start = start + len;
        }
        else if( type == cigar_character_type::onquery){
            // Nothing for now
        }
    }

    if( max_na > 0){
        vector<interval> merged;
        start = ivals[0].start;
        end = ivals[0].end;
        for(auto iter = std::begin(ivals); std::next(iter) != std::end(ivals); ++iter){
            if( std::next(iter)->start - iter->end < max_na){
                end = std::next(iter)->end;
            }
            else{
                merged.emplace_back(start,end-1);
                start = std::next(iter)->start;
                end = std::next(iter)->end;
            }
        }
        merged.emplace_back(start,end-1);
        result.insert(std::end(result),std::begin(merged),std::end(merged));
    }
    else{
        result.insert(std::end(result),std::begin(ivals),std::end(ivals));
    }
}

vector<interval> cigar2intervals( int alignment_start, cigar cig, int max_na){
    vector<interval> ivals;
    cigar2intervals( alignment_start, cig, max_na, ivals);
    return ivals;
}
//1       ensembl exon    69055   70108   .       +       .       gene_id "ENSG00000186092"; gene_version "6"; transcript_id "ENST00000335137"; transcript_version "4"; exon_number "1"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F5-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS30547"; exon_id "ENSE00002319515"; exon_version "2"; tag "basic"; transcript_support_level "NA (assigned to previous version 3)";
// 1       havana  transcript      65419   71585   .       +       .       gene_id "ENSG00000186092"; gene_version "6"; transcript_id "ENST00000641515"; transcript_version "2"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F5-202"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "basic";


IITree<locus, exon> build_exon_forest( const std::string &path_to_gtf){

    std::unordered_set<std::string> regular_chr= {"1","2","3","4","5","6","7","8","9","10",
        "11","12","13","14","15","16","17","18","19","20",
        "21","22","X","Y","MT"};
    for( auto r: regular_chr){
        regular_chr.insert("chr" + r);
    }
    IITree<locus, exon> exon_forest;
    std::set<std::string> lkup;

    bool keep_nc = true;
    bool keep_alt = false;
    bool keep_mito = false;

    ifstream gtf_stream( path_to_gtf);
    if(!gtf_stream.is_open()){
        std::cerr << "failed to open " << path_to_gtf << std::endl;
        exit(-1);
    }
    std::string line;

    while(getline(gtf_stream,line)){
        if(line[0]=='#'){ //Comment
            continue;
        }
        std::vector<std::string> tabs = rsplit(line,"\t");
        if(tabs[2] != "exon"){
            continue;
        }

        std::vector<std::string> fields = rsplit(tabs[8],";");

        bool coding = false;

        if(!keep_mito){
            if(tabs[0]=="MT"){
                continue;
            }
        }

        if(!keep_alt){
            if(regular_chr.find(tabs[0]) == regular_chr.end()){
                continue;
            }
        }


        for(auto iter = fields.begin(); iter != fields.end(); iter++){
            if( iter->find("transcript_biotype") != std::string::npos){
                std::string _id = iter->substr(iter->find("e ")+3);
                _id.pop_back();
                if( _id != "protein_coding"){
                    coding = true;
                    break;            
                }

            } 
        }


        if((keep_nc) || coding){
            exon ex;
            ex.coding = coding;
            ex.chr = tabs[0];
            if( ex.chr.find("chr") == 0){
                ex.chr = ex.chr.substr(3);
            }
            ex.start = std::stoi(tabs[3]);
            ex.end = std::stoi(tabs[4]);
            ex.strand = tabs[6] == "-";
            for(auto iter = fields.begin(); iter != fields.end(); iter++){

                if( iter->find("gene_id") != std::string::npos){
                    std::string _id = iter->substr(iter->find("d ")+3);
                    _id.pop_back();
                    if( _id.find(".") != std::string::npos){
                        size_t dot_pos = _id.find(".");
                        _id = _id.substr(0,dot_pos);
                    }
                    ex.gene_id = _id;
                } 
                if( iter->find("transcript_id") != std::string::npos){
                    std::string _id = iter->substr(iter->find("d ")+3);
                    _id.pop_back();
                    if( _id.find(".") != std::string::npos){
                        size_t dot_pos = _id.find(".");
                        _id = _id.substr(0,dot_pos);
                    }
                    ex.transcript_id = _id;
                } 

                if( iter->find("exon_number") != std::string::npos){
                    std::string _id;
                    if( iter->find("\"")!=std::string::npos){
                        _id = iter->substr(iter->find("r ")+3);
                        _id.pop_back();
                    }else{
                        _id = iter->substr(iter->find("r ")+2);
                    }
                    ex.exon_number = std::stoi(_id);
                }
                if( iter->find("gene_name") != std::string::npos){
                    std::string _id = iter->substr(iter->find("e ")+3);
                    _id.pop_back();
                    ex.gene_name = _id;
                }
            }

            locus start(ex.chr,ex.start);
            locus end(ex.chr,ex.end);
            exon_forest.add(start,end,ex);
        }
    }

    exon_forest.index();

    return exon_forest;
}

chainer::chainer( const vector<std::pair<aligned_segment,exon>> &read_mappings, int permissability): 
    mappings(read_mappings),
    permissability(permissability){

        P.resize(mappings.size(),-1);
        D.resize(mappings.size(),0.0);
        done.resize(mappings.size(),false);
    }

chainer::chainer( const vector<aligned_segment> &read_aligs, const IITree<locus, exon> &exon_tree, int permissability): 
    read_aligs(read_aligs), 
    permissability(permissability){
        std::cerr << read_aligs.size() << "\n";
        for( aligned_segment a : read_aligs){
            vector<size_t> overlaps;


            locus start(a.chr,a.tmplt.start);
            locus end(a.chr,a.tmplt.end);
            exon_tree.overlap(start, end, overlaps);
            for(size_t index : overlaps){
                exon ex = exon_tree.data(index);
                mappings.push_back(std::make_pair(a,ex));

            }

        }
        P.resize(mappings.size(),-1);
        D.resize(mappings.size(),0.0);
        done.resize(mappings.size(),false);
    }


void chainer::compute_fragment_score(size_t index){
    if(done[index]){
        return;
    }

    double max_score_value = -1;
    size_t max_parent_id    = -1;

    const pair<aligned_segment,exon> &val = mappings[index];
    double my_score = val.first.tmplt.reciprocal_overlap(val.second)  * (val.first.query.end - val.first.query.start);

    if(!mappings[index].second.coding){
        my_score = my_score -1;
    }
    for(size_t parent_id = 0; parent_id < mappings.size(); ++parent_id){
        if(!is_parent(index, parent_id)){
            continue;
        }


        compute_fragment_score(parent_id);

        bool same_gene = mappings[parent_id].second.gene_id == mappings[index].second.gene_id;
        bool same_tcipt = mappings[parent_id].second.transcript_id == mappings[index].second.transcript_id;
        double parent_score = D[parent_id];
        double mscore = my_score;
        if(!same_gene){
            mscore = mscore/2;
        }
        else if(!same_tcipt){
            mscore = 5*mscore/6;
        }
        /*
        if( mappings[parent_id].first.chr == mappings[index].first.chr){
            mscore = mscore / std::log2(std::abs(mappings[parent_id].first.tmplt.end - mappings[index].first.tmplt.start));
        }
        else{
            mscore = mscore / std::log2(300000000);
        }*/
        mscore += parent_score;
        if(mscore > max_score_value){
            max_score_value = mscore;

            max_parent_id = parent_id;
        }
    }

    D[index] = max_score_value;
    P[index] = max_parent_id;
    done[index] = true;

}

bool chainer::easy(){
    if(mappings.size() < 2){
        return true;
    }
    return false;
}

vector<pair<aligned_segment,exon>> chainer::chain(){
    double max_score = -1.0;
    size_t max_index =   -1;


    for(size_t i =0; i< mappings.size() ; ++i){

        compute_fragment_score(i);
        if( D[i] > max_score){
            max_score = (double) D[i];
            max_index = i;

        }

    }
    vector<pair<aligned_segment, exon>> optimal_chain;

    while(max_index != (size_t)-1){
        optimal_chain.push_back(mappings[max_index]);
        max_index = P[max_index];
    }

    std::reverse(optimal_chain.begin(),optimal_chain.end());
    return optimal_chain;
}



bool exon_overlaps( exon a, exon b, double ratio){
    interval ia = a;
    interval ib = b;
    return ia.reciprocal_overlap(ib) > ratio;
}

bool chainer::is_parent(size_t child, size_t prospective_parent ){

    if( child == prospective_parent){ return false;}
    
    auto &child_map = mappings[child];
    auto &parent_map = mappings[prospective_parent];

    const aligned_segment &cival = child_map.first;
    const aligned_segment &paval = parent_map.first;

    int permissed = std::min(permissability, (cival.query.end - cival.query.start)/4);

    if(cival.query.start + permissed <= paval.query.end){
        return false;
    }

    if(cival.reverse_complemented == paval.reverse_complemented && (child_map.second.gene_id == parent_map.second.gene_id) ){

        if(cival.reverse_complemented){
            if(cival.tmplt.end + permissed >= paval.tmplt.start){
                return false;
            }
            if((child_map.second).end + permissed >= (parent_map.second).start){
                return false;
            }
        }
        else{
            if(cival.tmplt.start + permissed <= paval.tmplt.end){
                return false;
            }
            if((child_map.second).start + permissed <= (parent_map.second).end){
                return false;
            }
        }
    }
    else{
        /*
        if( exon_overlaps(child_map.second, parent_map.second,0.25)){
            return child_map.second.exon_id == parent_map.second.exon_id;
        }
        */

    }
    return true;
}

