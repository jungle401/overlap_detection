#pragma once

#include <biovoltron/file_io/fasta.hpp>

class ReadsManip {
 private:
  static biovoltron::istring get_hpc(const biovoltron::istring_view& ori_seq) {
    auto seq = biovoltron::istring();
    if (ori_seq[0] == 4) {
      /* n_table.push_back({0, 1}); */
    }
    else {
      /* pos.emplace_back(0); */
      /* npos.emplace_back(n_table.size()); */
      seq += ori_seq[0];
      /* run_len.emplace_back(1); */
    }

    for (auto i = uint32_t{1}; i < ori_seq.size(); i++) {
      constexpr int N = 4;
      if (ori_seq[i] == N) {
        if (ori_seq[i - 1] == N) {
          /* n_table.back()[1]++; */
        }
        else {
          /* n_table.push_back({i, i + 1}); */
        }
      } else {
        if (ori_seq[i - 1] != ori_seq[i]) {
          /* if (seq.size() % POS_INTV == 0) { */
            /* pos.emplace_back(i); */
            /* npos.emplace_back(n_table.size()); */
          /* } */
          seq += ori_seq[i];
          /* run_len.emplace_back(1); */
        } else {
          /* run_len.back()++; */
        }
      }
    }
    return seq;
  }
 public:
  static void homopolymer_compression (std::vector<biovoltron::FastaRecord<true>>& reads) {
    for (auto& r : reads) {
      r.seq = get_hpc(r.seq);
    }
  }
};
