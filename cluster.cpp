#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "cluster.h"
#include "graph.h"

#include <cxxopts.hpp>

namespace clustering{

    void str_split(std::string str, std::string delim, std::vector<std::string> &v){
        v.clear();
        unsigned long p1 = 0;
        unsigned long p2 = 0;

        while((p2 = str.find(delim,p1)) != std::string::npos){
            v.push_back(str.substr(p1,p2-p1));
            p1 = p2 + 1;
        }
        v.push_back(str.substr(p1));
    }

    std::string tab = "\t";
    struct interval{
        std::string contig;
        int start;
        int end;
        interval(std::string ctg, int a, int b): contig(ctg),start(a),end(b){

        }

        bool overlaps(const interval &that, int range) const {
            if(this->contig != that.contig){ return false;}
            return ((that.start >= this->start - range) && (that.start <= this->end + range)) ||
                ((that.end >= this->start - range) && (that.end <= this->end + range));
        }
        friend std::ostream& operator << (std::ostream &os, const interval &i);
    };

    struct interval_pair{
        interval A;
        interval B;
        interval_pair(interval a, interval b) : A(a), B(b){

        }    
        bool overlaps( const interval_pair &that, int relax) const{
            return this->A.overlaps(that.A, relax) && this->B.overlaps(that.B, relax);
        }

        friend std::ostream& operator << (std::ostream &os, const interval_pair &i);
    };

    std::ostream& operator << (std::ostream &os, const interval &i){
        os << i.contig << '\t' << i.start << '\t' << i.end;
        return os;
    }

    std::ostream& operator << (std::ostream &os, const interval_pair &i){
        os << i.A << '\t' << i.B;
        return os;
    }

    cxxopts::ParseResult parse_args(int argc, char **argv ){
        try{
            cxxopts::Options *options = new cxxopts::Options(argv[0], "Gene fusion short read clustering");

            options->add_options()
                ("s,sam", "Input short read SAM file",cxxopts::value<std::string>())
                ("f,gtf", "Input GTF file",cxxopts::value<std::string>())
                ("o,output", "Output prefix a folder in that path should exist", cxxopts::value<std::string>())
                //            ("t,threads", "Number of threads", cxxopts::value<unsigned>()->default_value("8"))
                ("h,help", "Prints help")
                ;
            cxxopts::ParseResult result = options->parse(argc, argv);

            if( result.count("h")){

                std::cerr << options->help({"","Mandatory"}) << std::endl;
                exit(0);
            }
            if(!result.count("s")){
                std::cerr << "SAM file is required" << std::endl;
                exit(-1);
            }
            if(!result.count("f")){
                std::cerr << "GTF file is required" << std::endl;
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

    int fusion_cluster(int argc, char **argv){
        auto opt = parse_args(argc, argv);

        std::ifstream samfile( opt["sam"].as<std::string>());
        if(! samfile.is_open()){
            std::cerr << "Couldn't open " << opt["sam"].as<std::string>() << std::endl;
            return -1;
        }

        std::string line;
        std::vector<std::string> fields;
        std::vector<interval_pair> interval_pairs;

        std::vector<std::string> sa_field;
        while( getline(samfile,line)){
            str_split(line,"\t",fields);
            for( std::string part : fields){
                if(part[0] == 'S' && part[1] == 'A' && part[2] == ':' && part[3] == 'Z'){
                    int len = fields[9].length();
                    str_split(part,",",sa_field);
                    int sign = ((sa_field[2] == "+")?1:-1);
                    interval_pairs.push_back(
                            interval_pair(interval(fields[2],stoi(fields[3]),stoi(fields[3])+len),interval(sa_field[0].substr(5),stoi(sa_field[1]),stoi(sa_field[1])+(sign*(50-len)))));
                    //                        {{fields[2],stoi(fields[3]),stoi(fields[3])+stoi(fields[4])},
                    //                        {sa_field[0].substr(5),stoi(sa_field[1]),stoi(sa_field[1])+(sign*30)}}
                    //                        );
                    break;
                }
            }

        }
        samfile.close();

        std::unordered_map<std::string,int> geneLen;
        std::ifstream gtffile(opt["gtf"].as<std::string>());
        while( getline(gtffile,line)){
            if( line[0] == '#') {continue;} // Comment
            str_split(line,"\t",fields);
            if( fields[2].find("gene") == std::string::npos){ continue;} // Not Gene

            int len =  stoi(fields[4]) - stoi(fields[3]);
            str_split(fields[8]," ", fields);

            std::string gene = fields[1].substr(1);

            gene.pop_back();
            gene.pop_back();

            geneLen[gene] = len;
        }
        gtffile.close();
        interval_pairs.erase(std::remove_if(interval_pairs.begin(),interval_pairs.end(), [&geneLen] (const interval_pair &a){
                    std::vector<std::string> fields;
                    str_split(a.A.contig,"::",fields);

                    int size = geneLen[fields[0]];
                    return (a.A.contig != a.B.contig)|| a.A.end > size || a.B.start < size + 100;
                    }),interval_pairs.end());
        qcgraph<interval_pair> graph(0.5,0.6,40);

        for( interval_pair &i: interval_pairs){  
            graph.add_vertex(&i);
            //std::cout << i << std::endl;
        }
        int edge_count = 0;
        for( auto i1 = interval_pairs.begin(); i1 != interval_pairs.end(); i1++){  
            for( auto i2 = std::next(i1); i2 != interval_pairs.end(); i2++){  
                if(i1->overlaps(*i2,350)){
                    edge_count++;
                    graph.add_edge(*i1,*i2);   
                }
            }
        }

        std::cerr << "Node count: " << interval_pairs.size()<< std::endl;


        std::cerr << "Edge count: " << edge_count << std::endl;
        std::vector<std::vector<interval_pair *>> clusters = graph.find_clusters();

        std::ofstream cmd("igvcommands.batch");
        cmd << "new\n";
        cmd << "genome cand.fa\n";
        cmd << "load lr.bam\n";
        cmd << "load sr.bam\n";
        cmd << "snapshotDirectory snapshots\n";

        int count = 0;
        for( std::vector<interval_pair *> vv : clusters){


            if( vv.size() < 2) continue;

            std::cout << vv.size() << std::endl;
            std::cout << "Cluster: " << count ++ << "\n";
            //std::exit(-1);

            for( interval_pair *i : vv){
                std::cout << *i << "\n";
                if( i->B.end - i->A.start <44000){
                    cmd << "goto " << i->A.contig << ":" << i->A.start - 2500 << "-" << i->B.end + 2500 << "\n";
                    cmd << "snapshot " << count  << std::endl;
                }
                else{

                    cmd << "goto " << i->A.contig << ":" << i->A.start - 15000 << "-" << i->A.end + 15000 << "\n";
                    cmd << "snapshot " << count  << "-A" <<std::endl;
                    cmd << "goto " << i->A.contig << ":" << i->B.start - 15000 << "-" << i->B.end + 15000 << "\n";
                    cmd << "snapshot " << count  << "-B" <<  std::endl;
                }

            }
        }
        cmd.close();
        return 0;
    }
}
