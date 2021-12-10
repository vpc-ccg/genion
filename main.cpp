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


#include <filesystem>

#define FMT_STRING_ALIAS 1
#include <fmt/format.h>
#include <zlib.h>
#include <cxxopts.hpp>

#include "paf.h"
#include "filters.h"
#include "candidate.h"
#include "kseq.h"
#include "reference.h"
#include "util.h"
#include "annotate.h"

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

cxxopts::ParseResult parse_args(int argc, char **argv, bool is_wg = false){
    try{
        cxxopts::Options *options = new cxxopts::Options(argv[0], "Gene fusion");

        options->add_options()
            ("r,reference", "Reference path, see mkref",cxxopts::value<std::string>())
            ("p,tpaf", "Long read transcriptome alignment paf path",cxxopts::value<std::string>())
            ("g,gpaf", "Long read whole genom e alignment paf path",cxxopts::value<std::string>())
            ("m,homology", "Homolog gene pairs csv",cxxopts::value<std::string>())
            ("s,transcriptome-self-align", "Self align tsv",cxxopts::value<std::string>())
            ("prefix-filter", "Maximum number of unaligned prefix bases", cxxopts::value<int>())
            ("suffix-filter", "Maximum number of unaligned suffix bases", cxxopts::value<int>())
            ("mid-filter", "Maximum number of unaligned bases between fusion gene alignments", cxxopts::value<int>())
            ("no-strand-switch", "Don't allow strand switches")
            ("o,output", "Output prefix for an existing path", cxxopts::value<std::string>())
            ("t,threads", "Number of threads", cxxopts::value<unsigned>()->default_value("8"))
            ("e,force", "Force run, overwrites files in the output folder")
            ("keep-noncoding", "Keep non-coding exons")
            ("h,help", "Prints help")
            ;
        cxxopts::ParseResult result = options->parse(argc, argv);

        int status = 0;
        if( result.count("h")){
            std::cerr << options->help({"","Mandatory"}) << std::endl;
            exit(-1);

        }
        if(!result.count("r")){
            std::cerr << "reference is required" << std::endl; 
            status |=2;
        }
        if(!is_wg && !result.count("p")){
            std::cerr << "Transcriptome Paf alignment file is required" << std::endl;
            status |=4;
        }

        if(!result.count("o")){
            std::cerr << "Output prefix is required" << std::endl;           
            status |=32;
        }

        if(!result.count("g")){
            std::cerr << "Whole Genome Paf alignment file is required" << std::endl;
            status |=128;
        }

        if(status){
            std::cerr << options->help({"","Mandatory"}) << std::endl;
            exit(-1);
        }


        return result;

    }
    catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}

int fusion_run(int argc, char **argv){

    cxxopts::Options options(argv[0], "Gene fusion");
    options.add_options()
        ("r,reference", "Reference path, see mkref",cxxopts::value<std::string>())
        ("p,tpaf", "Long read transcriptome alignment paf path",cxxopts::value<std::string>())
        ("g,gpaf", "Long read whole genom e alignment paf path",cxxopts::value<std::string>())
        ("m,homology", "Homolog gene pairs csv",cxxopts::value<std::string>())
        ("d,duplications", "genomicSuperDups.txt, unzipped",cxxopts::value<std::string>())//can be found at http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz
        ("s,transcriptome-self-align", "Self align tsv",cxxopts::value<std::string>())
        ("o,output", "Output prefix for an existing path", cxxopts::value<std::string>())
        ("min-support", "Minimum read support for fusion calls", cxxopts::value<int>()->default_value("3"))
        ("max-rt-distance", "Maximum distance between genes for read-through events", cxxopts::value<int>()->default_value("500000"))
        ("max-rt-fin", "Maximum value of chimeric-count / normal-count for read-through events", cxxopts::value<double>()->default_value("0.2"))
        ("non-coding", "Allow non-coding genes and transcripts while calling gene fusions", cxxopts::value<bool>()->default_value("false"))
        ("prefix-filter", "Maximum number of unaligned prefix bases", cxxopts::value<int>())
        ("suffix-filter", "Maximum number of unaligned suffix bases", cxxopts::value<int>())
        ("mid-filter", "Maximum number of unaligned bases between fusion gene alignments", cxxopts::value<int>())
        ("no-strand-switch", "Don't allow strand switches")
        ("t,threads", "Number of threads", cxxopts::value<unsigned>()->default_value("8"))
        ("e,force", "Force run, overwrites files in the output folder")
        ("keep-noncoding", "Keep non-coding exons")
        ("h,help", "Prints help")
        ;
    cxxopts::ParseResult opt = options.parse(argc, argv);

    int status = 0;
    if( opt.count("h")){
        std::cerr << options.help({"","Mandatory"}) << std::endl;
        exit(-1);

    }
    if(!opt.count("r")){
        std::cerr << "reference is required" << std::endl; 
        status |=2;
    }

    if(!opt.count("o")){
        std::cerr << "Output prefix is required" << std::endl;           
        status |=32;
    }

    if(!opt.count("g")){
        std::cerr << "Whole Genome Paf alignment file is required" << std::endl;
        status |=128;
    }

    if(status){
        std::cerr << options.help({"","Mandatory"}) << std::endl;
        exit(-1);
    }

    string output_prefix(opt["output"].as<std::string>());

    std::string gtf_path = fmt::format("{}/1.gtf",opt["r"].as<std::string>());
    auto exon_forest = build_exon_forest(gtf_path);
    ChainFilterMinSegment chain_filter(exon_forest,50);

    if( !opt["force"].as<bool>()){
        if(std::filesystem::exists(output_prefix)){
            cerr << fmt::format("Path to {} already exists.\n",output_prefix);
            return -1;
        }
    }

    if(!std::filesystem::exists(output_prefix)){
        std::filesystem::create_directory(output_prefix);
    }

    std::filesystem::create_directory(output_prefix);
    std::string cdna_path = fmt::format("{}/reference.cdna.fa.gz",opt["r"].as<std::string>());
    map<string, gene_t> gene_info = init_gene_info(cdna_path);

    unordered_set<pair<string,string>> ref_self_align_info = init_sa_homolog_info(opt["s"].as<std::string>());
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
    vector<Candidate> fusion_candidates;
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
                    fusion_candidates.push_back(cand);
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
            fusion_candidates.push_back(cand);
        }
    }
    file_paf.close();

    annotate::annotate_calls_direct(
            output_prefix,
            opt["reference"].as<std::string>(),
            opt["duplications"].as<std::string>(),
            fusion_candidates,
            gene_counts,
            opt["min-support"].as<int>(),
            cnts[0], cnts[1]+cnts[2],
            opt["max-rt-distance"].as<int>(), opt["max-rt-fin"].as<double>(),
            !opt["non-coding"].as<bool>()
            );
    std::cerr << cnts[0] << "\t" << cnts[1] << "\t" << cnts[2] << "\n";
    return 0;
}


int fusion_filter_only_wg(int argc, char **argv){
    auto opt = parse_args(argc,argv,true);
    string output_prefix(opt["output"].as<std::string>());

    std::string gtf_path = fmt::format("{}/1.gtf",opt["r"].as<std::string>());
    auto exon_forest = build_exon_forest(gtf_path);
    ChainFilterMinSegment chain_filter(exon_forest,50);

    if( !opt["force"].as<bool>()){
        if(std::filesystem::exists(output_prefix)){
            cerr << fmt::format("Path to {} already exists.\n",output_prefix);
            return -1;
        }
    }

    if(!std::filesystem::exists(output_prefix)){
        std::filesystem::create_directory(output_prefix);
    }

    std::filesystem::create_directory(output_prefix);
    std::string cdna_path = fmt::format("{}/reference.cdna.fa.gz",opt["r"].as<std::string>());
    map<string, gene_t> gene_info = init_gene_info(cdna_path);

    unordered_set<pair<string,string>> ref_self_align_info = init_sa_homolog_info(opt["s"].as<std::string>());
    WholeGenomeSelfAlignFilter wg_self_align_filter(std::move(ref_self_align_info));

    vector<ofstream > streams;

    int cnt = 0;

    ofstream feature_table( fmt::format("{}/feature_table.tsv",output_prefix));
    streams.push_back(ofstream(fmt::format("{}/{:02d}-transcriptome-primary-alignment-count-filter-reads.list",output_prefix,cnt++))); 
    streams.push_back(ofstream(fmt::format("{}/{:02d}-wg-homology-filter-reads.list",output_prefix,cnt++))); 
    streams.push_back(ofstream(fmt::format("{}/candidate-reads.list",output_prefix)));

    ifstream file_paf(opt["g"].as<std::string>());
    if(file_paf.is_open() == false){
        cerr << fmt::format( "[ERROR] Cannot open file: {}",opt["g"].as<std::string>()) << std::endl;
        return EXIT_FAILURE;
    }
    
    paf_t map1;

    vector<paf_t> map_list;

    ofstream chains( fmt::format("{}/chains.txt",output_prefix));
    ofstream bad_chains( fmt::format("{}/bad_chains.txt",output_prefix));
    std::string lastId = "-1";
    map_list.clear();
    std::array<int,3> cnts{{0,0,0}};
    while(get_next_paf(file_paf, map1, exon_forest, true)){
        if(map1.qName != lastId){
            if(map_list.size() > 0){
                Candidate cand;
                cand.set_wg(map_list,false);

                ofstream *chainer;
                if(apply_filter(cand, &chain_filter)){
                    cnts[0]++;
                    log_candidate(streams[0], lastId, cand);
                    chainer = &bad_chains;
                } 
                else if(apply_filter(cand, &wg_self_align_filter)){
                    cnts[1]++;
                    log_candidate(streams[1], lastId, cand);
                    chainer = &bad_chains;
                }
                else{
                    cnts[2]++;
                    log_candidate(streams[2], lastId, cand);
                    chainer = &chains;
                }
                cand.print_filter_status_wg(feature_table,lastId);
                *chainer << lastId << "\t" << cand.canonical.size() << "\n";
                cand.print_chains(chainer);
            }
            map_list.clear();
            lastId = map1.qName;
        }
        map_list.push_back(map1);
    }
    if(map_list.size() > 0)
    {
        Candidate cand;
        cand.set_wg(map_list,false);

        ofstream *chainer;
        if(apply_filter(cand, &chain_filter)){
            log_candidate(streams[0], lastId, cand);
            chainer = &bad_chains;
        } 
        else if(apply_filter(cand, &wg_self_align_filter)){
            log_candidate(streams[1], lastId, cand);
            chainer = &bad_chains;
        }
        else{
            log_candidate(streams[2], lastId, cand);
            chainer = &chains;
        }
        cand.print_filter_status_wg(feature_table,lastId);
        *chainer << lastId << "\t" << cand.canonical.size() << "\n";
        cand.print_chains(chainer);
    }
    file_paf.close();
    for(unsigned i = 0; i < streams.size(); i++){
        streams[i].close();
    }

    std::cerr << cnts[0] << "\t" << cnts[1] << "\t" << cnts[2] << "\n";
    return 0;
}



int fusion_filter(int argc, char **argv){

    auto opt = parse_args(argc,argv);
    string output_prefix(opt["output"].as<std::string>());

    if( !opt["force"].as<bool>()){
        if(std::filesystem::exists(output_prefix)){
            cerr << fmt::format("Path to {} already exists.\n",output_prefix);
            return -1;
        }
    }

    if(!std::filesystem::exists(output_prefix)){
        std::filesystem::create_directory(output_prefix);
    }

    std::string cdna_path = fmt::format("{}/reference.cdna.fa.gz",opt["r"].as<std::string>());
    map<string, gene_t> gene_info = init_gene_info(cdna_path);

    std::cerr << "Reading Long Read paf file" << std::endl;
    // process paf file
    ifstream file_paf(opt["p"].as<std::string>());
    if(file_paf.is_open() == false)
    {
        cerr << "[ERROR] Cannot open file:" << opt["p"].as<std::string>() << std::endl;
        return EXIT_FAILURE;
    }

    std::string gtf_path = fmt::format("{}/1.gtf",opt["r"].as<std::string>());
    auto exon_forest = build_exon_forest(gtf_path);
    ChainFilter chain_filter(exon_forest);
    // get the first paf entry from the file

    paf_t map1;
    vector<paf_t> map_list;
    unordered_map<string,string> readId2cand;
    unordered_map<string,vector<string>> cand2readids;

    vector<AbstractFilter *> filters  {{     
        new TranscriptomePrimaryAlignmentCountFilter(2),      
        new TranscriptomePalindromeFilter(gene_info),
        new TranscriptomeSecondaryPalindromeFilter(0.7,0.8)
    }};
    for( const cxxopts::KeyValue &kv : opt.arguments()){
        if(kv.key() == "homology"){
            unordered_map<pair<string,string>,double> homolog_info = init_homolog_info(opt["m"].as<std::string>());
            filters.push_back(new HomologyFilter(std::move(homolog_info),0.7,0.8,0));
        }
        else if(kv.key() == "transcriptome-self-align"){
            unordered_set<pair<string,string>> ref_self_align_info = init_sa_homolog_info(opt["s"].as<std::string>());
            filters.push_back(new TranscriptomeReferenceSelfAlignFilter(std::move(ref_self_align_info)));
        }
        else if(kv.key() == "prefix-filter"){
            filters.push_back(new TranscriptomeProperPrefixFilter(opt["prefix-filter"].as<int>()));
        }
        else if(kv.key() == "suffix-filter"){
            filters.push_back(new TranscriptomeProperSuffixFilter(opt["suffix-filter"].as<int>()));
        }
        else if(kv.key() == "mid-filter"){
            filters.push_back(new TranscriptomeStrictMidFilter(opt["mid-filter"].as<int>()));
        }
        else if(kv.key() == "no-strand-switch"){
            filters.push_back(new StrandSwitchFilter(gene_info));
        }
    } 

    ofstream feature_table( fmt::format("{}/feature_table.tsv",output_prefix));

    unordered_map<string, Candidate> id2cand;

    vector<ofstream > streams;
    int cnt = 0;

    ofstream unmapped_fasta_stream( fmt::format("{}/{:02d}-unmapped-reads.list",output_prefix,cnt++));
    for(AbstractFilter *filter :filters){
        string _filter_id(typeid(*filter).name()+2);
        for( size_t index = 0; index < _filter_id.size(); index++){
            char _c = _filter_id[index];
            if(isupper(_c)){
                string _to;
                if( index == 0){   
                    _to = fmt::format("{}",(char)tolower(_c));
                }else{
                    _to = fmt::format("-{}",(char)tolower(_c));
                }
                _filter_id.replace(index,1,_to);
            }
        }
        streams.push_back(
                ofstream( fmt::format("{}/{:02d}-{}-reads.list",output_prefix,cnt++,_filter_id))
                );
    }

    string lastId = "-1";
    while(get_next_paf(file_paf, map1, false))
    {
        if(map1.qName != lastId)
        {
            if(map_list.size() > 0)
            {
                unsigned index;
                Candidate cand(map_list);
                for( index = 0; index < filters.size(); index++){
                    if(apply_filter(cand,filters[index])){
                        streams[index] << lastId << "\n";       
                    }
                }
                if(index == filters.size() && !cand.filtered){
                    id2cand[lastId] = cand;
                }
                else{
                    cand.print_filter_status(feature_table,lastId);
                }
            }
            map_list.clear();
            lastId = map1.qName;
        }
        map_list.push_back(map1);
    }
    if(map_list.size() > 0)
    {   
        Candidate cand(map_list);
        unsigned index;
        for(index =0 ; index < filters.size();index++){
            if(apply_filter(cand,filters[index])){
                streams[index] << map1.qName << "\n";       
            }
        }

        if(index == filters.size() && !cand.filtered){
            id2cand[map1.qName] = cand;
        }
        else{
            cand.print_filter_status(feature_table,lastId);
        }
    }

    file_paf.close();

    for(unsigned i = 0; i < streams.size(); i++){
        streams[i].close();
    }
    unmapped_fasta_stream.close();

    streams.clear();

    file_paf.open(opt["g"].as<std::string>());
    if(file_paf.is_open() == false)
    {
        cerr << fmt::format( "[ERROR] Cannot open file: {}",opt["g"].as<std::string>()) << std::endl;
        return EXIT_FAILURE;
    }

//    WholeGenomePrimaryAlignmentCountFilter wg_single_primary_filter(2);


    streams.push_back(ofstream(fmt::format("{}/{:02d}-wg-alignment-filter-reads.list",output_prefix,cnt++))); 
    streams.push_back(ofstream(fmt::format("{}/candidate-reads.list",output_prefix)));

    ofstream chains( fmt::format("{}/chains.txt",output_prefix));
    ofstream bad_chains( fmt::format("{}/bad_chains.txt",output_prefix));
    lastId = "-1";
    map_list.clear();
    while(get_next_paf(file_paf, map1, exon_forest, true)){
        if(map1.qName != lastId){
            if(map_list.size() > 0){
                auto candp = id2cand.find(lastId);
                if(candp == id2cand.end()){
                    map_list.clear();
                    lastId = map1.qName;
                    map_list.push_back(map1);
                    continue;
                }
                auto cand = candp->second;

                cand.set_wg(map_list,false);
                ofstream *chainer;
                
                if(apply_filter(cand, &chain_filter)){

                    streams[0] << lastId << "\n";
                    for(auto p = cand.transcriptome.cprimary_begin();
                            p != cand.transcriptome.cprimary_end();
                            ++p){
                        streams[0] << "\t" << p->tName;
                    }
                    streams[0] << "\t;";
                 //   << "\t" << cand.transcriptome.prefix().tName << "\t" <<
                   //     cand.transcriptome.suffix().tName;
                    for(auto p : cand.canonical){
                        streams[0] << "\t" << p.second.gene_id ;
                    }
                    streams[0] << "\n";
                    chainer = &bad_chains;
                }
                else{

                    streams[1] << lastId;
                    for(auto p = cand.transcriptome.cprimary_begin();
                            p != cand.transcriptome.cprimary_end();
                            ++p){
                        streams[1] << "\t" << p->tName;
                    }
                    streams[1] << "\t;";
                 //   << "\t" << cand.transcriptome.prefix().tName << "\t" <<
                   //     cand.transcriptome.suffix().tName;
                    for(auto p : cand.canonical){
                        streams[1] << "\t" << p.second.gene_id ;
                    }
                    streams[1] << "\n";
                    chainer = &chains;
                }

                *chainer << lastId << "\t" << cand.canonical.size() << "\n";
                for( auto pp : cand.canonical){
                    *chainer << "\t" << pp.first.tmplt.start << "\t" << pp.first.tmplt.end << "\t" <<
                        pp.first.chr << "\t" <<  pp.first.query.start << "\t" << pp.first.query.end << "\t" << 
                        pp.first.reverse_complemented << "\t" << 
                        pp.second.chr << "\t" << pp.second.start << "\t" << pp.second.end << "\t" << pp.second.strand  << "\t" <<
                        pp.second.gene_id << "\t" << pp.second.transcript_id << "\t" << pp.second.exon_number << "\n";
                }
            }
            map_list.clear();
            lastId = map1.qName;
        }
        map_list.push_back(map1);
    }
    if(map_list.size() > 0)
    {
        auto candp = id2cand.find(lastId);
        if( candp != id2cand.end()){
            auto cand = candp->second;
            cand.set_wg(map_list,false);
            if(apply_filter(cand,&chain_filter)){
                    streams[0] << lastId << "\n";
            }
            else{
                    streams[1] << lastId;
                    for(auto p = cand.transcriptome.cprimary_begin();
                            p != cand.transcriptome.cprimary_end();
                            ++p){
                        streams[1] << "\t" << p->tName;
                    }
                    streams[1] << "\t;";
                    for(auto p : cand.canonical){
                        streams[1] << "\t" << p.second.gene_id ;
                    }
                    streams[1] << "\n";
            }
        }
    }

    for(auto &pair: id2cand){
        pair.second.print_filter_status(feature_table,pair.first);
    }


    file_paf.close();
    for(unsigned i = 0; i < streams.size(); i++){
        streams[i].close();
    }

    streams.clear();

    id2cand.clear();

    return EXIT_SUCCESS;    
}

int unknown_command_exit(){
    std::cerr << fmt::format("\n{}\n\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n",
            "Usage:   fusion <command> [options]",
            "Commands:",
            "---------",
            "run",
            "filter",
            "wg-filter",
            "mkref",
            "annotate") << std::endl;
    return -1;
}

int main(int argc, char **argv){
    if(argc < 2){
        return unknown_command_exit();
    }
    string tool(argv[1]);


    if(tool == "run"){
        return fusion_run(argc-1,argv+1);   
    }
    if(tool == "filter"){
        return fusion_filter(argc-1,argv+1);   
    }
    if(tool == "wg-filter"){
        return fusion_filter_only_wg(argc-1,argv+1);   
    }
    if(tool == "mkref"){
        return reference::mkref(argc-1,argv+1);
    }
    if(tool == "annotate"){
        return annotate::annotate_calls(argc-1,argv+1);
    }
    if(tool == "intersect-gtf"){
        auto forest = build_exon_forest(argv[2]);
        string line;
        while(getline(cin, line)){
            vector<string> fields = rsplit(line,"\t");
            string ch = fields[0];
            int start = stoi(fields[1]);
            int end = stoi(fields[2]);
            vector<size_t> overlaps;
            forest.overlap(locus(ch,start),locus(ch,end),overlaps);
            for( size_t ov : overlaps){
                std::cout << ch << "\t" << start << "\t" << end << "\t" << forest.data(ov) << "\n";
            }
            overlaps.clear();
            return 0;
        }

    }
    return unknown_command_exit();
}
