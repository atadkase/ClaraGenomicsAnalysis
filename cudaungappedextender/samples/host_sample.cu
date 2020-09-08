/*
* Copyright 2020 NVIDIA CORPORATION.
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include <claraparabricks/genomeworks/cudaungappedextender/cudaungappedextender.hpp>
#include <claraparabricks/genomeworks/io/fasta_parser.hpp>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string>
#include <vector>

using namespace claraparabricks::genomeworks;
using namespace claraparabricks::genomeworks::cudaungappedextender;

int main(int argc, char* argv[])
{
    const int32_t input_xdrop = 910;
    const bool input_no_entropy = false;
    const int32_t score_threshold = 3000;

    // Define an alphabet for the sequences to be processed
    Alphabet alphabet  = make_alphabet("acgt");

    // Fasta query and target files
    std::string target_file_path = "../data/example.fa";
    std::unique_ptr<io::FastaParser> fasta_parser_target =
        io::create_kseq_fasta_parser(alphabet, target_file_path, 0, false);
    // Assumes that only one sequence is present per file
    SequenceVector target_sequences = fasta_parser_target->get_sequence_by_id(0);

    std::string query_file_path = "../data/example.fa";
    std::unique_ptr<io::FastaParser> fasta_parser_query =
        io::create_kseq_fasta_parser(alphabet, query_file_path, 0, false);
    // Assumes that only one sequence is present per file
    SequenceVector query_sequences = fasta_parser_query->get_sequence_by_id(0);

    // CSV SeedPairs file - Each row -> query_position_in_read_,
    // target_position_in_read_
    std::string seed_pairs_file_path = "../data/example_seed_pairs.csv";

    std::vector<SeedPair> h_seed_pairs;
    // Following function loops through all seed_pairs in the SeedPairs csv and returns
    // results in
    // the passed vector
    parse_SeedPairs(seed_pairs_file_path, h_seed_pairs);

    // Following sections TBD based on encoding
    ScoreMatrix score_matrix(a);
    score_matrix('a','a') = score;
    score_matrix('a','c') = score;
    score_matrix('a','g') = score;
    score_matrix('a','t') = score;
    score_matrix('c','c') = score;
    score_matrix('c','g') = score;
    score_matrix('c','t') = score;
    score_matrix('g','g') = score;
    score_matrix('g','t') = score;
    score_matrix('t','t') = score;

    // Create a stream for async use
    CudaStream stream0 = make_cuda_stream();
    // Create an ungapped extender object
    std::unique_ptr<UngappedExtender> ungapped_extender =
        std::make_unique<UngappedExtender>(0, score_matrix, input_xdrop,
                                           input_no_entropy, stream0.get());
    // Launch the ungapped extender host function
    int32_t query_idx  = 0;
    int32_t target_idx = 0;
    ungapped_extender->extend_async(query_sequences, query_idx, target_sequences, target_idx, score_threshold, h_seed_pairs);

    // Wait for ungapped extender to finish
    ungapped_extender->sync();

    // Get results
    const std::vector<ScoredSegmentPair>& segments =
        ungapped_extender->get_scored_segment_pairs();
    int32_t i = 0;
    for (const auto& segment : segments)
    {
        std::cout << "Segment: " << i << "Length: " << segment.length
                  << "Score: " << segment.score << std::endl;
        std::cout << "Position in query: "
                  << segment.seed_pair.query_position_in_read << std::endl;
        std::cout << "Position in target: "
                  << segment.seed_pair.target_position_in_read << std::endl;
        i++;
    }

    return 0;
}
