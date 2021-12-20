#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
#include <utility>

#include <future>




#define FMT_STRING_ALIAS 1
#include <fmt/format.h>
#include <zlib.h>
#include <cxxopts.hpp>

#include "paf.h"
#include "filters.h"
#include "candidate.h"
#include "kseq.h"

#include "util.h"
#include "annotate.h"
#include "gtf.h"

#ifndef GENION_VERSION
#define GENION_VERSION "1.0.1"
#endif

using namespace std;
KSEQ_INIT(gzFile, gzread);

#if 1
void log_candidate(std::ofstream &, const std::string &, Candidate &){
#else
void log_candidate(std::ofstream &stream, const std::string &lastId, Candidate &cand){
//Make this part optional for LGBM

    stream << lastId << "\n";
    for(auto p = cand.transcriptome.cprimary_begin();
            p != cand.transcriptome.cprimary_end();
            ++p){
        stream << "\t" << p->tName;
    }
    stream << "\t;";
    for(auto p : cand.canonical){
        stream << "\t" << p.second.gene_id ;
    }
    stream << "\n";
#endif
}

bool apply_filter(Candidate &cand, AbstractFilter *filter){

    bool ret = cand.filtered;
    bool status = (*filter)(cand);
    cand.add_filter_status(status);
    return !ret && !status;
}



//Read dispatcher after WG alignment step
void dispatch_wgs_reads(const string &lr_path, 
        const unordered_map<string,
        unsigned> &id2index,
        vector<ofstream> &streams,
        const unordered_map<string,
        Candidate> &id2cand){

    gzFile file_ref = gzopen(lr_path.c_str(), "r");
    kseq_t *read_fq = kseq_init(file_ref);
    int val = 0;
    int cou = 0;
    while (kseq_read(read_fq) >= 0)
    {
        std::string id(read_fq->name.s);
        cou ++;

        if(id2cand.find(id) == id2cand.end()){
            continue;
        }
        if(id2index.find(id) == id2index.end()){

            continue;
        }

        val++;
        unsigned stream_index = id2index.find(id)->second;
        Candidate cand = id2cand.find(id)->second;
        streams[stream_index] << ">" << id << "\t" << cand.transcriptome.prefix().tName << "\t" << cand.transcriptome.suffix().tName << "\n";
        streams[stream_index] << read_fq->seq.s << "\n";
    }
    cerr << fmt::format("{}/{}",val,cou) << std::endl;
    kseq_destroy(read_fq);
    gzclose(file_ref);
}
//Read dispatcher before WG alignment step
void dispatch_reads(const string &lr_path, const unordered_map<string,unsigned> &id2index, vector<ofstream> &streams, ofstream &unmapped_stream, const unordered_map<string, Candidate> &id2cand){
    gzFile file_ref = gzopen(lr_path.c_str(), "r");
    kseq_t *read_fq = kseq_init(file_ref);
    int val = 0;
    int cou = 0;
    while (kseq_read(read_fq) >= 0)
    {
        std::string id(read_fq->name.s);
        cou ++;
        if(id2index.find(id) == id2index.end()){

            unmapped_stream << ">" << id  << "\n";
            unmapped_stream << read_fq->seq.s << "\n";
            continue;
        }
        if(id2cand.find(id) == id2cand.end()){
            unmapped_stream << ">" << id  << "\n";
            unmapped_stream << read_fq->seq.s << "\n";
            continue;
        }
        val++;
        unsigned stream_index = id2index.find(id)->second;
        Candidate cand = id2cand.find(id)->second;
        if(cand.transcriptome.primary_count() > 1){
            streams[stream_index] << ">" << id << "\t" 
                << cand.transcriptome.prefix().tName << "\t" << cand.transcriptome.suffix().tName << "\n";
        }
        else{
            streams[stream_index] << ">" << id << "\t" << cand.transcriptome.prefix().tName << "\t-" << "\n";
        }
        streams[stream_index] << read_fq->seq.s << "\n";
    }
    cerr << fmt::format("{}/{}",val,cou) << std::endl;
    kseq_destroy(read_fq);
    gzclose(file_ref);
}
auto init_gene_info_gtf( string gtf_path){

    // build gene_id -> gene_name map


    map<string, gene_t> gene_info;
    map<string, string> transcript2gene;
    std::ifstream gt_file(gtf_path);
    string buffer;
    buffer.reserve(1000);

    while(std::getline(gt_file, buffer)){
        if(buffer[0] == '#'){
            continue;
        }
        gtf entry(buffer);
        gene_t &current = gene_info[entry.info["gene_id"]];
        switch(entry.type){
            case gtf::entry_type::gene:
                current.gene_name = entry.info["gene_name"];
                current.chr       = entry.chr;
                current.chr_start = entry.start;
                current.chr_end   = entry.end;
                break;
            case gtf::entry_type::transcript:
                transcript2gene[ entry.info["transcript_id"]] = entry.info["gene_id"];
                break;
            case gtf::entry_type::exon:
                break;
            default: ;// Skip others
        }
    }

    return make_tuple(gene_info, transcript2gene);
}


map<string, gene_t> init_gene_info( string ref_path){

    // build gene_id -> gene_name map
    gzFile file_ref = gzopen(ref_path.c_str(), "r");
    if(file_ref == NULL)
    {
        cerr << "[ERROR] Cannot open file: " << ref_path << std::endl;
        exit(EXIT_FAILURE);
    }
    map<string, gene_t> gene_info;
    std::cerr << "Reading Transcriptome FASTA" << std::endl;

    kseq_t *seq_trans = kseq_init(file_ref);
    while (kseq_read(seq_trans) >= 0)
    {
        // cout<< "here" << endl;
        vector<string> comm, info, gene_fields;
        // 
        get_fields(seq_trans->comment.s, comm);
        str_split(comm[2], ':', info);
        string gene_id = info[1].substr(0,15);
        if(gene_info.count(gene_id) == 0) // This is O(N)
        {
            str_split(comm[1], ':', info);
            gene_info[gene_id].chr = info[2];
            gene_info[gene_id].chr_start = str2type<int32_t>(info[3]);
            gene_info[gene_id].chr_end = str2type<int32_t>(info[4]);
            gene_info[gene_id].rev_strand = (str2type<int32_t>(info[5]) < 0 ? 1 : 0);
            // 
            str_split(comm[5], ':', info);
            gene_info[gene_id].gene_name = info[1];
        }
        else // update boundaries if needed
        {
            str_split(comm[1], ':', info);
            int tmp_start = str2type<int32_t>(info[3]);
            int tmp_end = str2type<int32_t>(info[4]);
            if(tmp_start < gene_info[gene_id].chr_start)
                gene_info[gene_id].chr_start = tmp_start;
            if(tmp_end > gene_info[gene_id].chr_end)
                gene_info[gene_id].chr_end = tmp_end;
        }
    }
    kseq_destroy(seq_trans);
    gzclose(file_ref);
    return gene_info;
}

unordered_map<pair<string,string>,double> init_homolog_info( string csv_path){
    string line;
    ifstream file_h(csv_path);
    unordered_map<pair<string,string>,double> homologs;

    while(getline(file_h,line)){
        vector<string> fields = rsplit(line,",");
        homologs[make_pair(fields[0],fields[6])] = stod(fields[11]);
    }
    file_h.close();
    return homologs;
}

unordered_set<pair<string,string>> init_sa_homolog_info( string tsv_path, map<string, string> &transcript2genes){
    string line;
    ifstream file_h(tsv_path);
    unordered_set<pair<string,string>> homologs;

    while(getline(file_h,line)){
        vector<string> fields = rsplit(line,"\t");
        try{
            string g1 = transcript2genes.at(fields[0]);
            string g2 = transcript2genes.at(fields[1]);

            homologs.insert(make_pair(g1,g2));
            homologs.insert(make_pair(g2,g1));
        }
        catch(std::out_of_range &e){ //Some transcripts in the cdna are not in gtf, but not a problem for us.
 //           cerr << fields[0] << "\t" << fields[1] << "\tNot in Reference!\n";
        }
    }
    file_h.close();
    return homologs;
}
unordered_set<pair<string,string>> init_sa_homolog_info( string tsv_path){
    string line;
    ifstream file_h(tsv_path);
    unordered_set<pair<string,string>> homologs;

    while(getline(file_h,line)){
        vector<string> fields = rsplit(line,"\t");
        homologs.insert(make_pair(fields[0],fields[2]));
        homologs.insert(make_pair(fields[2],fields[0]));
    }
    file_h.close();
    return homologs;
}




int fusion_run(int argc, char **argv){

    cxxopts::Options options(argv[0], "GENe fusION");
    options.add_options()
        ("gtf", "GTF annotation path",cxxopts::value<std::string>())
        ("i,input", "Input fast{a,q} file",cxxopts::value<std::string>())
        ("g,gpaf", "Long read whole genom e alignment paf path",cxxopts::value<std::string>())
        ("d,duplications", "genomicSuperDups.txt, unzipped",cxxopts::value<std::string>())//can be found at http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz
        ("s,transcriptome-self-align", "Self align tsv",cxxopts::value<std::string>())
        ("o,output", "Output prefix for an existing path", cxxopts::value<std::string>())
        ("min-support", "Minimum read support for fusion calls", cxxopts::value<size_t>()->default_value("3"))
        ("max-rt-distance", "Maximum distance between genes for read-through events", cxxopts::value<int>()->default_value("500000"))
        ("max-rt-fin", "Maximum value of chimeric-count / normal-count for read-through events", cxxopts::value<double>()->default_value("0.2"))
        ("non-coding", "Allow non-coding genes and transcripts while calling gene fusions", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Prints help")
        ("v,version", "Prints version")
        ;
    cxxopts::ParseResult opt = options.parse(argc, argv);

    vector<string> mandatory_args {{"gtf","output", "gpaf", "duplications"}};


    if( opt.count("v")){
        std::cout << GENION_VERSION << "\n";
        return 0;
    }
    if( opt.count("h")){
        std::cerr << options.help({"","Mandatory"}) << std::endl;
        return 0;
    }
    int aindex = 1;
    int status = 0;
    for( string &arg : mandatory_args){
        if( opt.count(arg) == 0){
            std::cerr << arg << " is required!" << std::endl;
            status |=aindex;
        }
        aindex = aindex << 1;
    }


    if(status){
        std::cerr << options.help({"","Mandatory"}) << std::endl;
        exit(-1);
    }

    string output_prefix(opt["output"].as<std::string>());

    std::string gtf_path = opt["gtf"].as<std::string>();
    auto exon_forest = build_exon_forest(gtf_path);
    ChainFilterMinSegment chain_filter(exon_forest,50);


//    std::string cdna_path = fmt::format("{}/reference.cdna.fa.gz",opt["r"].as<std::string>());
    auto [ gene_info,transcript2gene] = init_gene_info_gtf(gtf_path);

    unordered_set<pair<string,string>> ref_self_align_info = init_sa_homolog_info(opt["s"].as<std::string>(), transcript2gene);
    WholeGenomeSelfAlignFilter wg_self_align_filter(std::move(ref_self_align_info));



    ifstream file_paf(opt["g"].as<std::string>());
    if(file_paf.is_open() == false){
        cerr << fmt::format( "[ERROR] Cannot open file: {}",opt["g"].as<std::string>()) << std::endl;
        return EXIT_FAILURE;
    }
    
    paf_t map1;

    vector<paf_t> map_list;

    std::string lastId = "-1";
    map_list.clear();
    std::array<int,3> cnts{{0,0,0}};

    //vector<Candidate> fusion_candidates;
    unordered_map<string,Candidate> fusion_candidates;
    std::unordered_map<string, size_t> gene_counts;
    while(get_next_paf(file_paf, map1, exon_forest, true)){
        if(map1.qName != lastId){
            if(map_list.size() > 0){
                Candidate cand{lastId};
                cand.set_wg(map_list,false);

                if(apply_filter(cand, &chain_filter)){
                    cnts[0]++;
                    std::unordered_set<string> genes;
                    for(auto &p : cand.canonical){
                        genes.insert(p.second.gene_id);
                    }
                    auto g1 = std::begin(genes);
                    if(g1!=genes.end()){
                        gene_counts[*g1]+=1;
                    }
                } 
                else if(apply_filter(cand, &wg_self_align_filter)){
                    cnts[1]++;
                }
                else{
                    cnts[2]++;
                    fusion_candidates[lastId] = cand;
                }
            }
            map_list.clear();
            lastId = map1.qName;
        }
        map_list.push_back(map1);
    }
    if(map_list.size() > 0)
    {
        Candidate cand{lastId};
        cand.set_wg(map_list,false);
        if(apply_filter(cand, &chain_filter)){
            cnts[0]++;
        } 
        else if(apply_filter(cand, &wg_self_align_filter)){
            cnts[1]++;
        }
        else{
            cnts[2]++;
            fusion_candidates[lastId] = cand;

        }
    }
    file_paf.close();


    vector<Candidate> cand_vec;
    double max_percent = 0.65;
    if(opt["input"].count() > 0 ){
        string fastq_filename {opt["input"].as<string>() };
        
        gzFile file_ref = gzopen(fastq_filename.c_str(), "r");
        kseq_t *read_fq = kseq_init(file_ref);

        while (kseq_read(read_fq) >= 0){
            std::string id(read_fq->name.s);
            
            auto iter = fusion_candidates.find(id);
            if(iter == fusion_candidates.end()){
                continue;
            }
            std::vector<std::pair<aligned_segment, exon>> new_canonical;
            std::set<string> new_genes;
            for(auto &seg : iter->second.canonical){
                int start = seg.first.query.start;
                int end = seg.first.query.end;
                map<char, int> base_counts;
                for(int i = start; i < end; ++i){
                    char c = read_fq->seq.s[i];
                    base_counts[c] += 1;
                }
                bool flag = true;
                int max_bc = (end - start) * max_percent;
                for( const std::pair<char, int> bc : base_counts){
                    if(bc.second > max_bc){
                        flag = true;
                        break;
                    }
                }
                if( flag){
                    new_canonical.push_back(seg);
                    new_genes.insert(seg.second.gene_id);
                }
            }
            if( new_canonical.size() > 0 && new_genes.size() > 1){
                iter->second.canonical = new_canonical;
                cand_vec.push_back(iter->second);
            }
        }
    }

    annotate::annotate_calls_direct(
            output_prefix,
            opt["gtf"].as<std::string>(),
            opt["duplications"].as<std::string>(),
            cand_vec,
            gene_counts,
            opt["min-support"].as<size_t>(),
            cnts[0], cnts[1]+cnts[2],
            opt["max-rt-distance"].as<int>(), opt["max-rt-fin"].as<double>(),
            !opt["non-coding"].as<bool>()
            );
    std::cerr << "Normal reads: " << cnts[0] << "\nHomologous chimeric reads: " << cnts[1] << "\nFusion candidate chimeric reads: " << cnts[2] << "\n";
    return 0;
}

int main(int argc, char **argv){
   return fusion_run(argc,argv);   
}
