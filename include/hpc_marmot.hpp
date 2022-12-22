#pragma once

#include <biovoltron/utility/archive/serializer.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <fstream>

namespace biovoltron {

class HPC {
 public:
  const int POS_INTV = 64;

  istring seq;
  std::vector<uint8_t> run_len;
  std::vector<uint32_t> pos;
  std::vector<uint16_t> npos;
  std::vector<std::array<uint32_t, 2>> n_table;

  void load(std::ifstream &fin) {
    int POS_INTV;
    fin.read(reinterpret_cast<char*>(&POS_INTV), sizeof(POS_INTV));
    assert(POS_INTV == HPC::POS_INTV);

    Serializer::load(fin, seq);
    Serializer::load(fin, run_len);
    Serializer::load(fin, pos);
    Serializer::load(fin, npos);
    Serializer::load(fin, n_table);
  }

  void save(std::ofstream &fout) {
    fout.write(reinterpret_cast<const char*>(&POS_INTV), sizeof(POS_INTV));
    Serializer::save(fout, seq);
    Serializer::save(fout, run_len);
    Serializer::save(fout, pos);
    Serializer::save(fout, npos);
    Serializer::save(fout, n_table);
  }
  
  void build(istring_view ori_seq) {
    if (ori_seq[0] == 4)
      n_table.push_back({0, 1});
    else {
      pos.emplace_back(0);
      npos.emplace_back(n_table.size());
      seq += ori_seq[0];
      run_len.emplace_back(1);
    }

    for (auto i = uint32_t{1}; i < ori_seq.size(); i++) {
      constexpr int N = 4;
      if (ori_seq[i] == N) {
        if (ori_seq[i - 1] == N)
          n_table.back()[1]++;
        else
          n_table.push_back({i, i + 1});
      } else {
        if (ori_seq[i - 1] != ori_seq[i]) {
          if (seq.size() % POS_INTV == 0) {
            pos.emplace_back(i);
            npos.emplace_back(n_table.size());
          }
          seq += ori_seq[i];
          run_len.emplace_back(1);
        } else {
          run_len.back()++;
        }
      }
    }
  }

  // i: hpc position
  // return origin position
  auto compute_pos(size_t i) {
    // 靠近他左邊的位置 pos_beg
    const auto pos_beg = i / POS_INTV;
    if (POS_INTV == 1 or i % POS_INTV == 0)
      return pos[pos_beg];

    // 靠近他右邊的位置 pos_end
    const auto pos_end = pos_beg + 1;
    auto beg = pos_beg * POS_INTV;
    auto end = pos_end * POS_INTV;
    //  beg--i------- end
    if (pos_end == pos.size() or i - beg < end - i) {
      // backward search ------>
      auto cnt = uint32_t{};
      auto it = n_table.begin() + npos[pos_beg];
      for (; beg < i; beg++) {
        cnt += run_len[beg];
        if (it != n_table.end() and pos[pos_beg] + cnt == (*it)[0]) {
          cnt += (*it)[1] - (*it)[0];
          it++;
        }
        // TODO: 可以不用一直查 n_table, 如果 npos[pos_beg] npos[pos_end] 不一樣才要查 n_table
      }

      return pos[pos_beg] + cnt;
    } else {
      auto it = n_table.begin() + npos[pos_end];
      // forward search
      if (it == n_table.begin())
        it = n_table.end();
      else
        it--;

      auto cnt = uint32_t{};
      for (; i < end; end--) {
        if (it != n_table.end() and pos[pos_end] - cnt == (*it)[1]) {
          cnt += (*it)[1] - (*it)[0];
          if (it == n_table.begin())
            it = n_table.end();
          else
            it--;
        }
        cnt += run_len[end - 1];
      }

      return pos[pos_end] - cnt;
    }
  }
};
}
