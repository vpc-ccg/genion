#ifndef CHAIN_CHAIN
#define CHAIN_CHAIN
#include <iostream>
#include <algorithm>
#include <vector>

#include "paf.h"
#include "cigar.h"

#include <IITree.h>
#include "locus.h"
class chainer{

    std::vector<aligned_segment> read_aligs;
    //    set<exon> exon_set;
    //    vector<exon> exons;

    std::vector<std::pair<aligned_segment,exon>> mappings;
    int permissability;
    std::vector<bool> done;
    std::vector<size_t> P;
    std::vector<double> D;
    //    unordered_map<interval, set<exon>> interval2exons;


    void compute_fragment_score(size_t index);

    bool is_parent(size_t child, size_t prospective_parent );
    public:
    chainer( const std::vector<aligned_segment> &read_aligs, const IITree<locus, exon> &exon_tree, int permissability);
    chainer( const std::vector<std::pair<aligned_segment,exon>> &read_mappings, int permissability);

    bool easy();

    std::vector<std::pair<aligned_segment,exon>> chain();

};


std::vector<interval> cigar2intervals( int alignment_start, cigar cig, int max_na);
IITree<locus, exon> build_exon_forest(const std::string &path_to_gtf);



#endif
