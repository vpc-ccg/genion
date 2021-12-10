#ifndef GENION_ANNOTATE
#define GENION_ANNOTATE

#include "candidate.h"
#include <vector>
#include <unordered_map>
#include <string>

namespace annotate{
    int annotate_calls(int argc, char **argv);
    int annotate_calls_direct(
            const std::string &output_path,
            const std::string &reference_path,
            const std::string &duplications,
            const std::vector<Candidate> &candidates,
            const std::unordered_map<std::string, size_t> gene_counts,
            int min_support,
            int total_normal_count, int total_chimer_count,
            int maxrtdistance,      double maxrtfin,
            bool only_coding
    );
}
#endif
