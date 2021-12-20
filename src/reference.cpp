#include "reference.h"

#include <iostream>
#include <filesystem>
#include <fstream>
#include <unordered_set>


#include <fmt/format.h>
#include <fmt/ostream.h>
#include <cxxopts.hpp>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);

#include "util.h"

namespace reference{

    cxxopts::ParseResult parse_args(int argc, char **argv ){
        try{
            cxxopts::Options *options = new cxxopts::Options(argv[0], "Fusion reference building");

            options->add_options()
                ("c,cdna", "Transcriptome reference",cxxopts::value<std::string>())
                ("g,dna", "Genome reference",cxxopts::value<std::string>())
                ("f,gtf", "Input GTF file",cxxopts::value<std::string>())
                ("o,output", "Output prefix a folder in that path should exist", cxxopts::value<std::string>())
                ("h,help", "Prints help")
                ("m,keep-mitochondrial", "Keep mitochondrial transcripts")
                ("e,force", "Force run")
                ("p,keep-non-coding", "Keep non-coding transcripts")
                ("a,keep-alternative-contigs", "Keep genes from alternative contigs")
                ("T,TTTTTTTT", "Attach poly-A tail to transcripts")
                ;
            cxxopts::ParseResult result = options->parse(argc, argv);
            int ret = 0;
            if( result.count("h")){
                std::cerr << options->help({"","Mandatory"}) << std::endl;
                ret |=1;
            }
            if(!result.count("c")){
                std::cerr << "cDNA file is required" << std::endl;
                ret |=2;
            }
            if(!result.count("f")){
                std::cerr << "GTF file is required" << std::endl;
                ret |=4;
            }
            if(!result.count("o")){
                std::cerr << "Output path is required" << std::endl;
                ret |=8;
            }
            if(ret != 0){
                std::cerr << options->help({"","Mandatory"}) << std::endl;
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

    std::string copy_gtf(std::string prefix, std::string gtf){
        std::string new_gtf =  fmt::format("{}/1.gtf",prefix );

        std::string cmd;
        if(gtf.substr(gtf.find_last_of(".") +1) == "gz"){
            cmd = fmt::format("gzip -c -d {} > {}",gtf,new_gtf);
        }else{
            cmd = fmt::format("cp {} {}",gtf,new_gtf);
        }   
        int ret = std::system(cmd.c_str());
        if(ret){
            std::cerr << fmt::format("Command {} failed with error code {}",cmd,ret) << std::endl;
        }
        return new_gtf;
    }

    std::string get_gene_from_comment(std::string comment){
        
        std::vector<std::string> fields = rsplit(comment," ");
        for(auto iter = fields.begin(); iter != fields.end(); iter++){
            if( iter->find("gene") != std::string::npos){
                std::string _id = iter->substr(iter->find(":") + 1, iter->find_last_of("."));

                return _id;

            } 
        }
        return "";
    }

    std::string copy_dna(std::string prefix, std::string ref){
        std::string ext = "fa";
        if(ref.substr(ref.find_last_of(".") +1) == "gz"){
            ext = "fa.gz";
        }
        std::string new_ref = fmt::format("{}/0.dna.{}",prefix,ext);
        std::string cmd = fmt::format("cp {} {}", ref, new_ref );
        int ret = std::system(cmd.c_str());
        if(ret){
            std::cerr << fmt::format("Command {} failed with error code {}",cmd,ret) << std::endl;
        }
        return new_ref;
    }

    std::string copy_cdna(std::string prefix, std::string ref){
        std::string ext = "fa";
        if(ref.substr(ref.find_last_of(".") +1) == "gz"){
            ext = "fa.gz";
        }
        std::string new_ref = fmt::format("{}/0.cdna.{}",prefix,ext);
        std::string cmd = fmt::format("cp {} {}", ref, new_ref );
        int ret = std::system(cmd.c_str());
        if(ret){
            std::cerr << fmt::format("Command {} failed with error code {}",cmd,ret) << std::endl;
        }
        return new_ref;
    }
    std::unordered_set<std::string> init_valid_transcripts(std::string gtf_path, bool keep_mito, bool keep_nc, bool keep_alt){
        std::unordered_set<std::string> regular_chr= {"1","2","3","4","5","6","7","8","9","10",
                                                        "11","12","13","14","15","16","17","18","19","20",
                                                        "21","22","X","Y","MT"};
        std::ifstream gtf_file(gtf_path);
        if(!gtf_file.is_open()){
            std::cerr << "[ERROR] Cannot open file:" << gtf_path << std::endl;
            exit(-1);
        } 

        std::unordered_set<std::string> valid_set;
        std::string line;
        while(std::getline(gtf_file, line)){
            if(line[0]=='#'){ //Comment
                continue;
            }
            std::vector<std::string> tabs = rsplit(line,"\t");
            if(tabs[2] != "transcript"){
                continue;
            }

            std::vector<std::string> fields = rsplit(tabs[8],";");
            std::string chr = tabs[0];
            if(chr.find("chr") != std::string::npos){
                chr = chr.substr(3);
            }
            bool flag = true;
            
            if(!keep_mito){
                if(chr=="MT"){
                    continue;
                }
            }
            if(!keep_alt){
               if(regular_chr.find(chr) ==regular_chr.end()){
                    std::cout << chr << "\n";
                   continue;

               }
            }
            if(!keep_nc){
                for(auto iter = fields.begin(); iter != fields.end(); iter++){
                    if( iter->find("gene_biotype") != std::string::npos){
                        std::string _id = iter->substr(iter->find("e ")+3);
                        _id.pop_back();
                        

                        if( _id != "protein_coding"){
                            flag = false;
                            break;            
                        }

                    } 
                }
            }

            if(flag){
                for(auto iter = fields.begin(); iter != fields.end(); iter++){
                    if( iter->find("transcript_id") != std::string::npos){
                        std::string _id = iter->substr(iter->find("d ")+3);
                        _id.pop_back();
                        if( _id.find(".") != std::string::npos){
                            size_t dot_pos = _id.find(".");
                            _id = _id.substr(0,dot_pos);
                        }
                        valid_set.insert(_id);

                    } 
                }
            }
        }
        std::cerr << valid_set.size() << std::endl;
        return valid_set;
    }

    int mkref(int argc, char **argv){
        auto opt = parse_args(argc, argv);
        std::string output_prefix(opt["output"].as<std::string>());
        
        if( ! opt["force"].as<bool>()){
            if(std::filesystem::exists(output_prefix)){
                std::cerr << fmt::format("Path to {} already exists.\n",output_prefix);
                return -1;
            }
        }
        if(!std::filesystem::exists(output_prefix)){
            std::filesystem::create_directory(output_prefix);
        }



        std::string gtf_path = copy_gtf(output_prefix,opt["gtf"].as<std::string>());
        std::string cdna_path = opt["cdna"].as<std::string>();
//        std::string dna_path = copy_dna(output_prefix,opt["dna"].as<std::string>());


        std::unordered_set<std::string> valid_transcripts = init_valid_transcripts(gtf_path, opt["m"].count()>0, opt["p"].count() >0, opt["a"].count() > 0);

        gzFile cdna_file = gzopen(cdna_path.c_str(),"r");
        if(!cdna_file){
            std::cerr << fmt::format("Cannot open file: {}!",opt["cdna"].as<std::string>()) << std::endl;
            return -1;
        }
        kseq_t *seq_tc = kseq_init(cdna_file);
        std::string out_file(fmt::format("{}/reference.cdna.fa",output_prefix));
        std::ofstream formatted_fasta(out_file);
        while(kseq_read(seq_tc) >= 0){
            std::string id(seq_tc->name.s);
            std::string comment;
            if(id.find("|") != std::string::npos){
                std::vector<std::string> fields = rsplit(id, "|");
                for ( auto f: fields){
                    std::cerr << f << " || ";
                } 
                std::cerr << "\n";

                id = fields[0];
                comment = "gene:"  + fields[1];
                //>ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|
            }else{
                comment = seq_tc->comment.s;
            }
            std::string short_id = id.substr(0,id.find_last_of("."));
            auto tptr = valid_transcripts.find(short_id);
            if(tptr== valid_transcripts.end()){
    //            std::cout << id << std::endl;
    //

                continue;
            }


            std::string seq(seq_tc->seq.s);
            std::string gene = get_gene_from_comment(comment);

            formatted_fasta << fmt::format(">{}_{}\t{}\n",gene,id,comment);

            formatted_fasta << seq;
        
            if( opt["T"].count()>0){
                formatted_fasta << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
            }

            formatted_fasta <<  "\n";
        }
        formatted_fasta.close();
        gzclose(cdna_file);
        int ret_val;
        //gzip reference
        //int ret_val = std::system(fmt::format("gzip {}",out_file).c_str());
        //if(ret_val != 0){
        //    std::cerr << "Error while compressing the cdna reference!" << std::endl;
        //    return ret_val;
        //}
        
//        ret_val = std::system(fmt::format("gzip {}",dna_path).c_str());
//        if(ret_val != 0){
//            std::cerr << "Error while compressing the dna reference!" << std::endl;
//            return ret_val;
//        }
        
        ret_val = std::system(fmt::format("sum {} > {}/cdna.sum",cdna_path,output_prefix).c_str());
        if(ret_val != 0){
            std::cerr << "Error while calculating original cdna checksum!" << std::endl;
            return ret_val;
        }

        ret_val = std::system(fmt::format("gzip {}",out_file).c_str());
        if(ret_val != 0){
            std::cerr << "Error while compressing the cdna reference!" << std::endl;
            return ret_val;
        }
       
        return 0;
    }
}
