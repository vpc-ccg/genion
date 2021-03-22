
#include <utility>

#include <unordered_set>
#include <set>
#include <numeric>
#include <vector>
#include <tuple>
#include <functional>
#include <sstream>
#include <fstream>
#include <iostream>

#include <cmath>

#include <cxxopts.hpp>
#include <IITree.h>

#include "annotate.h"
#include "locus.h"
#include "util.h"

#include <iostream>
namespace annotate{



    std::ostream& operator<<(std::ostream& os, const locus &lc){
        os << lc.chr << "\t" << lc.position;
        return os;
    }
    enum class SEQDIR{ forward, reverse, unknown};
    cxxopts::ParseResult parse_args(int argc, char **argv ){
        try{
            cxxopts::Options *options = new cxxopts::Options(argv[0], "Fusion candidate annotation");

            options->add_options()


                ("i,input", "Output path of Genion filter stage", cxxopts::value<std::string>())
                ("o,output", "Output path of Genion annotation stage", cxxopts::value<std::string>())
                ("q,read_dir", "Read direction tabular", cxxopts::value<std::string>())
                ("s,minsupport", "min support to flag PASS", cxxopts::value<size_t>()->default_value("3"))
                ("f,minfin", "min fin to flag PASS", cxxopts::value<double>()->default_value("0.001"))
                ("d,duplications", "genomicSuperDups.txt, unzipped",cxxopts::value<std::string>())//can be found at http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz
                ("r,reference", "Reference path used in filter stage",cxxopts::value<std::string>())
                ("c,keep_non_coding", "Keep non coding genes", cxxopts::value<bool>()->default_value("false"))
                ("h,help", "Prints help")
                ;
            cxxopts::ParseResult result = options->parse(argc, argv);
            int ret = 0;
            if( result.count("h")){
                std::cerr << options->help() << std::endl;
                ret |=1;
            }

            if(!result.count("i")){
                std::cerr << "input is required" << std::endl;
                ret |=8;
            }
            if(!result.count("o")){
                std::cerr << "output is required" << std::endl;
                ret |=16;
            }

            if(!result.count("d")){
                std::cerr << "Duplication annotation is required" << std::endl;
                ret |=32;
            }
            if(!result.count("r")){
                std::cerr << "reference is required" << std::endl;
                ret |=64;
            }

            if(ret != 0){
                std::cerr << options->help() << std::endl;
                exit(-1);
            }
            return result;
        }
        catch (const cxxopts::OptionException& e)
        {
            std::cerr << "error parsing options: " << e.what() << std::endl;
            exit(1);
        }
    }

//585     chr1    10000   87112   chr15:101906152 0       -       chr15   101906152       101981189       75037   11764   1000    N/A     N/A     N/A     N/A     align_both/0009/both0046049     77880   71      3611  74269   73743   526     331     195     0.992918        0.991969        0.00711601      0.00711937 
    IITree<locus, std::tuple<std::string, int, int, double> > read_duplication_annotation(std::string path){
        IITree<locus, std::tuple<std::string, int, int, double> > duplications;

        std::ifstream dup_file(path);
        if(!dup_file.is_open()){
            std::cerr << "[ERROR] Cannot open file:" << path << std::endl;
            exit(-1);
        } 
        std::string line;
        
        std::string ch;
        int start;
        int end;

        std::string m_ch;
        int m_start;
        int m_end;

        double frac_match;
        while(std::getline(dup_file, line)){
            std::vector<std::string> fields = rsplit(line,"\t");
            ch = fields[1];
            start = stoi(fields[2]);
            end = stoi(fields[3]);

            m_ch = fields[7];
            m_start = stoi(fields[8]);
            m_end = stoi(fields[9]);

            frac_match = stof(fields[26]);

            if(ch.find("chr")!= std::string::npos){
                ch = ch.substr(3);
            }
            if(m_ch.find("chr")!= std::string::npos){
                m_ch = m_ch.substr(3);
            }
            locus s_s(ch,start);
            locus s_e(ch,end);
            std::cerr << m_ch << "\t" << m_start << "\t" << m_end << "\t" << s_s.chr << "\t" << s_s.position << "\t" << s_e.chr << "\t" << s_e.position << "\tDUP" << "\n";
            duplications.add(s_s, s_e, std::make_tuple(m_ch,m_start,m_end,frac_match));
        }
        dup_file.close();
        duplications.index();
        return duplications;
    }
  
    class interval{

        public:
        std::string chr;
        int start;
        int end;
        bool reverse_strand;

        interval( const std::string &chr, int start, int end, bool reverse_strand) :
            chr(chr),
            start(start),
            end(end),
            reverse_strand(reverse_strand) {}
        std::pair<locus, locus> as_loci()const {
            return std::make_pair(locus(chr,start),locus(chr,end));
        }
        
        interval () : chr(""), start(-1), end(-1), reverse_strand(0) {}
        
        std::pair< locus,locus> as_loci(){
            return std::make_pair(locus(chr,start),locus(chr,end));
        }
        bool overlaps(const interval &other) const{
            if( chr != other.chr){
                return false;
            }
            if(start > other.start){
                return start < other.end;
            }
            return end > other.start;
        }
    };
//chr1    HAVANA          gene    65419   71585   .       +       .       gene_id "ENSG00000186092.6"; gene_type "protein_coding"; gene_name "OR4F5"; level 2; hgnc_id "HGNC:14825"; havana_gene "OTTHUMG00000001094.4";
//   1    ensembl_havana  gene    65419   71585   .       +       .       gene_id "ENSG00000186092"; gene_version "6"; gene_name "OR4F5"; gene_source "ensembl_havana"; gene_biotype "protein_coding";`
    

    class gene{
        public:
        interval range;

        bool reverse_strand;
        std::string gene_id;
        std::string gene_name;
        std::string gene_type;
        bool coding;

        gene( interval range, const std::string &gene_id, const std::string &gene_name, 
                const std::string &gene_type) 
            : range(range),
              gene_id(gene_id),
              gene_name(gene_name),
              gene_type(gene_type),
              coding(gene_type=="protein_coding")
        {}
        gene( interval range, const std::string &gene_id, const std::string &gene_name, 
                const std::string &gene_type, bool coding) 
            : range(range),
              gene_id(gene_id),
              gene_name(gene_name),
              gene_type(gene_type),
              coding(coding)
        {}
        gene(){}
    };

    class exon{
        public:
        interval range;
        std::string gene_id;  
        std::string transcript_id;  
        int exon_no;

        exon( interval range, const std::string &gene_id, const std::string &transcript_id, int exon_no) :
            range(range),
            gene_id(gene_id),
            transcript_id(transcript_id),
            exon_no(exon_no) {}
    };

    class candidate_read{

        public:

        std::string read_id;
        std::vector<std::pair<interval, exon> > blocks;
        candidate_read(const std::string &rid) : read_id(rid) {}
        std::vector<int> first_exons;

        std::map<std::string, locus> get_breakpoints(bool direction) const {
            
            std::map<std::string, locus> bps;
            
            std::string first_gene = blocks[0].second.gene_id;
            for(auto &block : blocks){

                std::string gene_id = block.second.gene_id;
                bool is_first = (gene_id == first_gene) == direction;

                bool reverse = block.first.reverse_strand;
                auto loci = block.first.as_loci();
                
                    //std::cerr << "RT\t" << reverse << "\n";
                auto gene_ptr  = bps.find(gene_id);
                if ( gene_ptr == bps.end()){
                    if(reverse != is_first){
                        bps.emplace(gene_id,loci.first);
                    }
                    else{
                        bps.emplace(gene_id,loci.second);
                    }
                }
                else{
                    int pos = gene_ptr->second.position;

                    if(reverse){
                        if(is_first){
                            if(pos < loci.second.position){
                                bps[gene_id] = loci.second;
                            }
                        }
                        else{
                            if(pos > loci.first.position){
                                bps[gene_id] = loci.first;
                            }
                        }
                    }
                    else{
                        if(is_first){
                            if(pos > loci.first.position){
                                bps[gene_id] = loci.first;
                            }

                        }
                        else{
                            if(pos < loci.second.position){
                                bps[gene_id] = loci.second;
                            }

                        }
                    }

                    
                   /* 
                    if(reverse == is_first){
                        if(pos > loci.first.position){
                            bps[gene_id] = loci.first;
                        }
                    }
                    else{
                        if(pos < loci.second.position){
                            bps[gene_id] = loci.second;
                        }
                    }
*/
                }
            }

            return bps;
        }
        void add_block(const std::string &line){
            std::vector<std::string> fields = rsplit(line, "\t");
//        154633180       154633212       X       60      92      1       X       154632470       154633182       ENSG00000272681 ENST00000598177 1
            int start = stoi(fields[1]);
            int end   = stoi(fields[2]);
            std::string chr = fields[3];
            bool reverse_strand = fields[6] == "1";

            int ex_start = stoi(fields[8]);
            int ex_end   = stoi(fields[9]);
            
            bool ex_rev_strand = fields[10] == "1";
            std::string gene_id  = fields[11];
            std::string transcript_id = fields[12];
            int exon_no               = stoi(fields[13]);
            if(exon_no == 1){
                first_exons.push_back(blocks.size());
            }
            
            
            interval alig(chr,start,end, reverse_strand);

            interval expos(chr,ex_start,ex_end, ex_rev_strand);
            exon ex(expos, gene_id, transcript_id, exon_no);

            blocks.push_back(std::make_pair(alig,ex));
        }
    };

    class candidate_fusion{

        public:
        int fg_count = 0;
        int lg_count = 0;
        std::string name;
        std::vector<candidate_read> forward;
        std::vector<candidate_read> backward;

        std::vector<candidate_read> no_first;
        std::vector<candidate_read> multi_first;

        std::vector<std::pair<interval, interval> > duplications;
        std::vector<std::pair<gene, gene> > gene_overlaps;

        std::map<std::string, interval> fusion_gene_intervals(){
            std::map<std::string, std::string> chrs;
            std::map<std::string, int> mins; 
            std::map<std::string, int> maxs;
            std::map<std::string, bool> rev;
            for(const auto &v : {forward, backward, no_first, multi_first}){
                for(const auto &c : v){
                    for(const auto &i_e : c.blocks){
                        const interval &i = i_e.first;
                        const exon &e = i_e.second;
                        int mn = mins[e.gene_id];
                        int mx = maxs[e.gene_id];
                        if(i.start < mn || mn == 0){
                            mins[e.gene_id] = i.start;
                        }
                        if(i.end > mx){
                            maxs[e.gene_id] = i.end;
                        }
                        chrs[e.gene_id] = i.chr;
                        rev[e.gene_id] = i.reverse_strand;
                    }
                }
            } 
            std::map<std::string, interval> ivals;
            for(const auto &k_v : mins){
                const std::string &key = k_v.first;
                const int &mn = k_v.second;
                const int &mx = maxs[key];
                const std::string &chr = chrs[key];
                bool rs = rev[key];
                ivals.emplace(key, std::move(interval(chr,mn,mx,rs)));
            }
            return ivals;
        }


        candidate_fusion() {}
    };

    auto dash_fold(const std::string &a, const std::string &b){
        return std::move(a) + "::" + b;
    }
    class fusion_manager{
        public:
        std::unordered_map<std::string, candidate_fusion> fusions;
        std::unordered_map<std::string, int> gene_counts;

        
        fusion_manager() {}

        void add_read(const candidate_read &read, const std::unordered_map<std::string, gene> &gene_annot,
                std::unordered_map<std::string, int> &last_exons,
                std::unordered_map<std::string, SEQDIR> &directions){


            std::set<std::string> gene_ids;
            std::map<std::string,int> gene_order;
            std::set<std::string> transcript_ids;

            SEQDIR dir = directions[read.read_id];
            int index = 0;
            for(auto i_and_e : read.blocks){
                auto ite = gene_ids.find(i_and_e.second.gene_id);
                if(ite == gene_ids.end()){
                    gene_order[i_and_e.second.gene_id] = index;
                    ++index;
                }
                gene_ids.insert(i_and_e.second.gene_id);

                if(gene_annot.find(i_and_e.second.gene_id) == gene_annot.end()){
                    std::cerr << i_and_e.second.gene_id << " is not in annotation!\n";
                }
//                std::cerr << gene_annot.find(i_and_e.second.gene_id)->second.gene_name << std::endl;
            
            }
            bool first_good = true;
            bool last_good = true;
            for(auto i_and_e : read.blocks){
                std::string gid = i_and_e.second.gene_id;
                std::string tid = i_and_e.second.transcript_id;
                int last_exon = last_exons[tid];
                int order = gene_order[gid];

                if(dir == SEQDIR::forward){
                    if(order != 0){
                        if( i_and_e.second.exon_no >= last_exon - 1){
                            first_good = false;
                        }
                    }
                    else{
                        if( i_and_e.second.exon_no == 2){
                            last_good = false;
                        }

                    }
                }else if (dir == SEQDIR::reverse){
                    if(order == 0){
                        if( i_and_e.second.exon_no == 2){
                            last_good = false;
                        }

                    }
                    else{
                        if( i_and_e.second.exon_no >= last_exon - 1){
                            first_good = false;
                        }

                    }

                }else{
                    last_good = false;
                    first_good = false;
                }
         
            }

             
  //          std::cerr << "\n";
            
            std::string fusion_name = "";
            for(const std::string &id : gene_ids){
                fusion_name += gene_annot.find(id)->second.gene_name + "::";
            }

            fusion_name.pop_back();
            fusion_name.pop_back();

            std::string fusion_id = std::accumulate( std::next(std::begin(gene_ids)), std::end(gene_ids), *(std::begin(gene_ids)), dash_fold);

            for( std::string gid : gene_ids){
                gene_counts[gid]+=1;
            }
            auto &cand = fusions[fusion_id];
            if(first_good){
                cand.fg_count+=1;
            }
            if(last_good){
                cand.lg_count+=1;
            }
            cand.name = fusion_name;
            int last_first = - 1;
            if(read.first_exons.size() > 1){
                cand.multi_first.push_back(read);
//                std::cerr << "Multiple first " << fusion_id << "\n";
                return;
            }
            if(read.first_exons.size() == 0){
                cand.no_first.push_back(read);
//                std::cerr << "No first " << fusion_id << "\n";
                return;
            }
            last_first = read.first_exons.back();
            if(read.blocks[last_first].second.gene_id == *(gene_ids.rbegin())){
                cand.forward.push_back(read);
            }
            else{
                cand.backward.push_back(read);
            }
        }
    };

    bool make_gene(const std::string &line, gene &g){
        std::vector<std::string> tabs = rsplit(line, "\t");
        std::string ch(tabs[0]);
        if(ch.find("chr")!=std::string::npos){
            ch = ch.substr(3);
        }
        interval range( tabs[0], stoi(tabs[3]), stoi(tabs[4]), tabs[6]=="-");
        if(tabs[2] != "gene"){
            return false;
        }
        std::vector<std::string> fields = rsplit(tabs[8], ";");
        std::string gene_id, gene_name, gene_type;

        for(auto iter = fields.begin(); iter != fields.end(); iter++){

            if( iter->find("gene_id") != std::string::npos){
                std::string _id = iter->substr(iter->find("d ")+3);
                _id.pop_back();
                if( _id.find(".") != std::string::npos){
                    size_t dot_pos = _id.find(".");
                    _id = _id.substr(0,dot_pos);
                }
                gene_id = _id;
            } 

            if( iter->find("gene_name") != std::string::npos){
                std::string _id = iter->substr(iter->find("e ")+3);
                _id.pop_back();
                gene_name = _id;
            }
            if( iter->find("gene_biotype") != std::string::npos ||
                iter->find("gene_type") != std::string::npos 
              ){
                std::string _id = iter->substr(iter->find("e ")+3);
                _id.pop_back();
                gene_type = _id;
            }
        }
   //     std::cerr << gene_type << "\tTYPE\n";
        g = gene(range, gene_id, gene_name, gene_type);
        return true;
    }

// umap of transcript to last exon id
    std::unordered_map<std::string, int> read_last_exons(std::string gtf_path){
        std::ifstream gtf_file(gtf_path);
        if(!gtf_file.is_open()){
            std::cerr << "[ERROR] Cannot open file:" << gtf_path << std::endl;
            exit(-1);
        } 


        std::unordered_map<std::string, int> int_map;
        std::string line;
        while(std::getline(gtf_file, line)){
            if(line[0]=='#'){ //Comment
                continue;
            }
            std::vector<std::string> tabs = rsplit(line, "\t");
            if(tabs[2] != "exon"){
                continue;
            }

            std::vector<std::string> fields = rsplit(tabs[8], ";");
            std::string transcript_id = "-1";
            int exon_number = -1;
            for(auto iter = fields.begin(); iter != fields.end(); iter++){

                if( iter->find("transcript_id") != std::string::npos){
                    std::string _id = iter->substr(iter->find("d ")+3);
                    _id.pop_back();
                    if( _id.find(".") != std::string::npos){
                        size_t dot_pos = _id.find(".");
                        _id = _id.substr(0,dot_pos);
                    }
                    transcript_id = _id;

                }
                if( iter->find("exon_number") != std::string::npos){
                    std::string _id = iter->substr(iter->find("r ")+3);
                    _id.pop_back();
                    if( _id.find(".") != std::string::npos){
                        size_t dot_pos = _id.find(".");
                        _id = _id.substr(0,dot_pos);
                    }
                    exon_number = stoi(_id);
                        
                }
            }
            if( transcript_id == "-1"){
                std::cerr << "Transcript doesn't have transcript_id\n";
            }
            if( exon_number > int_map[transcript_id]){
                int_map[transcript_id] = exon_number;
            }
        }
        gtf_file.close();
//        std::cerr << "Read " << int_map.size() << " transcript annotations" << std::endl;
        return int_map;
    }
    template<class K, class V,template<class,class> class MAP>
    std::vector<std::pair<K, K>> get_key_pairs( const MAP<K,V> &map){
        std::vector<std::pair<K,K>> pairs;
        for(auto iter = std::begin(map); iter != std::end(map); ++iter){
            for(auto inner = std::next(iter); inner !=std::end(map); ++inner){
                pairs.emplace_back(iter->first, inner->first);
            }
        }
        return pairs;
    }

    std::unordered_map<std::string, gene> read_gene_annotation(std::string gtf_path){
        std::ifstream gtf_file(gtf_path);
        if(!gtf_file.is_open()){
            std::cerr << "[ERROR] Cannot open file:" << gtf_path << std::endl;
            exit(-1);
        } 

        std::unordered_map<std::string, gene> valid_set;
        std::string line;
        while(std::getline(gtf_file, line)){
            if(line[0]=='#'){ //Comment
                continue;
            }
            gene g;
            if(make_gene(line, g)){
                valid_set.emplace( g.gene_id, g);
            }
        }
        gtf_file.close();
//        std::cerr << "Read " << valid_set.size() << " gene annotations" << std::endl;
        return valid_set;
    }



    void annotate_duplications_and_overlaps(fusion_manager &fm, const std::unordered_map<std::string, gene> &gene_annot, const std::string &dup_path){

        IITree<locus, std::tuple<std::string, int, int, double> > duplications = read_duplication_annotation(dup_path);;
        std::vector< size_t> overlaps;
        for( auto &cand : fm.fusions){
            std::map<std::string, interval> ivals = cand.second.fusion_gene_intervals();
            auto key_pairs = get_key_pairs(ivals);
                
            //Duplication annotation
            for(const auto &key_pair : key_pairs){
                const std::string &f = key_pair.first;
                const std::string &s = key_pair.second;
                
                const interval &i = (ivals.find(f))->second;
                const auto loci = i.as_loci();

               // std::cerr << loci.first.chr << "\t" << loci.first.position << "\t" << loci.second.chr << "\t" << loci.second.position << "\t" << f << "\t" << s << "\n";

                duplications.overlap(loci.first,loci.second, overlaps);

                for(size_t d : overlaps){
                    const auto &dup = duplications.data(d);
//                    std::cerr << std::get<0>(dup) << "\t" <<  std::get<1>(dup) << "\t" << std::get<2>(dup) << "\n";
                    interval l(std::get<0>(dup), std::get<1>(dup), std::get<2>(dup), 0);
                    interval &r =( ivals.find(s))->second; 
                    if(l.overlaps(r)){
                        cand.second.duplications.emplace_back(i,l);
                    }
                }
                overlaps.clear(); 
            }
            //X

            //Gene overlap annotation
            for(const auto &key_pair : key_pairs){
                const gene &f = gene_annot.find(key_pair.first)->second;
                const gene &s = gene_annot.find(key_pair.second)->second;
                if(f.range.overlaps(s.range)){
                    cand.second.gene_overlaps.emplace_back(f,s); 
                }             
            }
            //X
        }
    }


    std::unordered_map<std::string, size_t> count_genes( const std::string &feature_table_path, bool all = false){
        
        std::unordered_map<std::string, size_t>  count_table;  
        
        std::ifstream feature_file(feature_table_path);
        
        std::string line;
        while(std::getline(feature_file, line)){
            std::vector<std::string> fields =  rsplit(line, "\t");
            if(!all && stoi(fields[2]) !=0){ // Split Alignment
                continue;
            }
            

            std::string gene_id1 = fields[1].substr(0,15);
            std::string gene_id2 = fields[1].substr(17);
//            std::cerr << gene_id1 << " - " << gene_id2 << std::endl;
            if(all){
                count_table[gene_id1]+=1;
                if(gene_id1 != gene_id2){
                    count_table[gene_id2]+=1;
                }
            }
            else{
                if(gene_id1 == gene_id2){
                    count_table[gene_id1]+=1;
                }
            }

        }
        return count_table;
    }


    std::unordered_map<std::string, SEQDIR>  read_read_directions(const std::string &path){
        std::unordered_map<std::string, SEQDIR> directions; 
        std::ifstream dir_file(path);
        std::string line;

        while(std::getline(dir_file, line)){
            std::vector<std::string> fields = rsplit(line,"\t");
            if(fields[1] =="NONE"){
                directions[fields[0]] = SEQDIR::unknown;
                continue;
            }
            if(fields[1] == "A" && stoi(fields[2]) > 50){
                directions[fields[0]] = SEQDIR::reverse;
            }
            else if (fields[1] == "T" && stoi(fields[2]) < 50){
                directions[fields[0]] = SEQDIR::forward;
            }
            else{
                directions[fields[0]] = SEQDIR::unknown;
            }
        }
        return directions;
    }
    int annotate_calls(int argc, char **argv){
        auto  opt = parse_args(argc, argv);

        size_t min_support = opt["minsupport"].as<size_t>();
        double min_fin_score = opt["minfin"].as<double>();

        std::string input_prefix(opt["input"].as<std::string>());
        std::string reference_path(opt["reference"].as<std::string>());

        bool filter_non_coding = !opt["c"].as<bool>();

        std::string gtf_path = reference_path + "/1.gtf";
        std::unordered_map<std::string, gene> gene_annot = read_gene_annotation(gtf_path);
        
        std::unordered_map<std::string, int> last_exons = read_last_exons(gtf_path);

        std::unordered_map<std::string, SEQDIR> read_directions = read_read_directions(opt["read_dir"].as<std::string>());

        std::string chains_path = input_prefix + "/chains.fixed.txt";
        std::string candidate_path = input_prefix + "/candidate-reads.list";

        std::vector<candidate_read> candidates;

        std::ifstream chain_file(chains_path);
        
        std::string line;
        while(std::getline(chain_file, line)){
            std::vector<std::string> fields =  rsplit(line, "\t");
            int block_count = stoi(fields[1]);
            candidate_read cr(fields[0]);
            for(int i = 0; i < block_count; ++i){
                std::getline(chain_file, line);
                cr.add_block(line);
            }
            candidates.push_back(std::move(cr));
        }

        chain_file.close();
        fusion_manager fm;
        for( auto &cand : candidates){
            fm.add_read(cand, gene_annot, last_exons, read_directions);
        }
        
        annotate_duplications_and_overlaps(fm, gene_annot, opt["duplications"].as<std::string>());
        
        std::string feature_table_path = input_prefix + "/feature_table.tsv";

        std::unordered_map<std::string, size_t> gene_counts = count_genes(feature_table_path,false);

        for( const auto &cand : fm.fusions){

            int total_count = cand.second.forward.size() + cand.second.backward.size() + cand.second.multi_first.size() + cand.second.no_first.size();
            int total_count_putative_full_length = cand.second.forward.size() + cand.second.backward.size();
            std::vector<std::string> genes = rsplit(cand.first,"::");

            bool coding_flag = false;
            if( filter_non_coding){

                for( const std::string &g : genes){
                    auto gptr = gene_annot.find(g);
                    if(gptr == gene_annot.end()){
                        std::cerr << "Gene " << g << " is not in annotation!\n";
                        continue;
                    }
                    if(gptr->second.coding == false){
                        coding_flag = true;
                        break;
                    }
                }
             //   if(coding_flag){
            //        break;
           //     }
            }
            double gene_count_sum = 0;
            std::string gene_count_string = "";
            std::string idf_string = "";
            double total_idf = 0;
            for(const std::string &gene : genes){
//                std::cerr  << "COUNT " << gene << "\t" << gene_counts[gene] << " for " << cand.first << "\n";
                gene_count_sum += gene_counts[gene];
                gene_count_string+= std::to_string(gene_counts[gene]) + ";";
                idf_string+= std::to_string(fm.gene_counts[gene]-total_count) + ";";
                total_idf+= fm.gene_counts[gene] - total_count;
            }
            
            
            double tfidf_score = total_count * std::log(fm.fusions.size()/(1+total_idf/2));
            double tfidf_score_full_len = total_count_putative_full_length * std::log(fm.fusions.size()/(1+total_idf/2));
            

            int tcpflnz;
            if(total_count == 0){
                tcpflnz = 1;
            }
            else{
                tcpflnz = total_count;
            }

            double fin_score = genes.size() * total_count / (gene_count_sum+1);

            std::string pass_fail_code = "";
            if(coding_flag){
                pass_fail_code += ":noncoding";
            }

            if(cand.second.gene_overlaps.size() > 0){
                pass_fail_code += ":overlaps";
            }
            if(cand.second.duplications.size() > 0){
                pass_fail_code += ":segdup";
            }

            if( fin_score < min_fin_score){
                pass_fail_code += ":lowfin";
            }

            //if( cand.second.forward.size() + cand.second.backward.size() < min_support){
            if( cand.second.forward.size() + cand.second.backward.size() 
                    + cand.second.multi_first.size() < min_support){
                pass_fail_code += ":lowsup";
            }
            if( pass_fail_code == ""){
                pass_fail_code = "PASS";
            }
            else{
                pass_fail_code = "FAIL" + pass_fail_code;
            }
            
            //#FusionID(Ensembl) Forward-Support Backward-Support Multi-First-Exon No-First-Exon Genes-Overlap Segmental-Duplication-Count FusionName(Symbol) FiN-Score Pass-Fail-Status total-normal-count fusion-count normal-counts proper-normal-count proper-FiN-Score total-other-fusion-count other-fusion-counts ffigf-score proper-ffigf-score A B Anorm Bnorm 
            std::cout << cand.first << "\t" << cand.second.forward.size() << "\t"
                << cand.second.backward.size()  << "\t" << cand.second.multi_first.size() << "\t" << cand.second.no_first.size()
                << "\t" <<  cand.second.gene_overlaps.size() << "\t" <<  cand.second.duplications.size()
                << "\t" << cand.second.name << "\t" << genes.size() * total_count / (gene_count_sum+1)
                << "\t" <<  pass_fail_code
                << "\t" << gene_count_sum << "\t" << total_count <<  "\t"  << gene_count_string << "\t"
                << total_count_putative_full_length << "\t" << genes.size() * total_count_putative_full_length / ( gene_count_sum + 1)
                << "\t" <<  total_idf << "\t" << idf_string << "\t" << tfidf_score << "\t" << tfidf_score_full_len
                << "\t" << cand.second.fg_count  << "\t" << cand.second.lg_count << "\t"
                << 1.0 * cand.second.fg_count / tcpflnz<< "\t"
                << 1.0 * cand.second.lg_count / tcpflnz << "\n";
        }
       
        std::string bp_file_path = opt["output"].as<std::string>() + "/breakpoints.tsv";

        std::ofstream bp_file(bp_file_path);


        for( const auto &cand : fm.fusions){
            std::string fusion_id = cand.first;
            std::map<std::string, std::vector<locus>> breakpoints;

            bool is_forward = true;
            for(const auto &ff : {cand.second.forward, cand.second.backward}){//, cand.second.no_first, cand.second.multi_first}){
                for(const auto &fus : ff){
                    //for(const auto &bp_pair : fus.get_breakpoints( &ff == &cand.second.forward)){
                    for(const auto &bp_pair : fus.get_breakpoints(is_forward)){


                        bp_file << fus.read_id <<"\t" << fusion_id << "\t" << bp_pair.first << "\t" << bp_pair.second << "\n";

                        breakpoints[bp_pair.first].push_back(bp_pair.second);
                    }
                }
                is_forward = false;
            }
            //for(const auto &bps : breakpoints){
            //    for(const auto &bp : bps.second){
            //        bp_file << fusion_id << "\t" << bps.first << "\t" << bp << "\n";
            //    }
            //}
        }
        bp_file.close();
        return 0;  
    }
       
} 
