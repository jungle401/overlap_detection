#pragma once
#include <algorithm>
#include <iostream>
#include <fstream>
#include <anchor_map.hpp>

using Anchor = std::vector<ReadOffset>;

enum class Strand {
  FORWARD,
  REVERSE,
};

class Detector {
 private:

  std::vector<std::pair<int, int>> dp;
  void dump(std::filesystem::path p, int r1len, int r2len) {
    auto fout = std::ofstream{p};
    assert(fout.is_open());
    fout << r1len << 'x';
    fout << r2len << '\n';
    for (const auto& [x, y] : dp) {
      fout << x << ',';
      fout << y << '\n';
    }
  }

  struct Dot {
    Dot(){}
    Dot (uint32_t _x, uint32_t _y) : x(_x), y(_y) { }
    uint32_t x;
    uint32_t y;
  };

  struct Segment {
    Segment (Dot _head, Dot _tail) : head(_head), tail(_tail) { }
    Dot head, tail;
  };

  struct Chain {
    Chain (Dot head_) : head(head_), tail(head_) {}
    Dot head, tail;
    uint32_t num_segs;
    float score = 0.f;
    /* to-be-abort */
    uint16_t cur_seg_len = 0;
  };

  struct ReferenceReadsScores {
    std::vector<std::vector<Chain>> forward_chains;
    std::vector<std::vector<Chain>> reverse_chains;
    ReferenceReadsScores (uint32_t reference_reads_size) {
      forward_chains.resize(reference_reads_size);
      reverse_chains.resize(reference_reads_size);
    }
  };

  template<Strand strand>
  typename std::enable_if<strand == Strand::FORWARD>::type scoring_basewise (
      ReferenceReadsScores& reference_reads_scores, 
      uint32_t query_read_index,
      auto shift_factor, 
      auto k,
      uint32_t base_index,
      const Anchor* anchor) {

    auto it_anchor = std::upper_bound(anchor->begin(), anchor->end(), query_read_index, [] (uint32_t query_read_index, const ReadOffset& rhs) {
        return query_read_index < rhs.read;
        });
    for (; it_anchor != anchor->end(); it_anchor = std::next(it_anchor)) {
      auto& read_index = it_anchor->read;
      auto read_index_rrs = read_index - query_read_index - 1;
      auto& offset = it_anchor->offset;
      auto append_to_some_chain = false;

      for (auto& chain : reference_reads_scores.forward_chains[read_index_rrs]) {
        if (base_index == chain.tail.x + 1 && offset == chain.tail.y + 1) {
          chain.cur_seg_len++;
          chain.tail.x += 1;
          chain.tail.y += 1;
          append_to_some_chain = true;
          break;
        } 
        if (offset <= chain.tail.y) {
          continue;
        }
        auto vert_dist = base_index - chain.tail.x;
        auto hori_dist = offset - chain.tail.y;
        auto& short_dist = vert_dist;
        auto& long_dist = hori_dist;
        if (short_dist > long_dist) { std::swap(short_dist, long_dist); }
        if (shift_factor * short_dist > long_dist) {
          /* account the last segment and append this pivot point as new segment to the chain. */
          chain.score += -std::log(std::pow(0.85, 2 * (k + chain.cur_seg_len)));
          chain.cur_seg_len = 0;
          chain.tail = Dot(base_index, offset);
          append_to_some_chain = true;
          /* break; */
        }
      }
      if (!append_to_some_chain) {
        /* start a new chain if nowhere to append */
        reference_reads_scores.forward_chains[read_index_rrs].push_back(Chain(Dot(base_index, offset)));
      }
    }
  };

  template<Strand strand>
  typename std::enable_if<strand == Strand::REVERSE>::type scoring_basewise (
      ReferenceReadsScores& reference_reads_scores, 
      uint32_t query_read_index,
      auto shift_factor, 
      auto k,
      uint32_t base_index,
      const Anchor* anchor) {

    auto it_anchor = std::upper_bound(anchor->begin(), anchor->end(), query_read_index, [] (uint32_t query_read_index, const ReadOffset& rhs) {
        return query_read_index < rhs.read;
        });

    for (; it_anchor != anchor->end(); it_anchor = std::next(it_anchor)) {
      auto& read_index = it_anchor->read;
      auto read_index_rrs = read_index - query_read_index - 1;
      auto& offset = it_anchor->offset;
      auto append_to_some_chain = false;

      for (auto& chain : reference_reads_scores.reverse_chains[read_index_rrs]) {
        if (base_index == chain.tail.x + 1 && offset == chain.tail.y - 1) {
          chain.cur_seg_len++;
          chain.tail.x += 1;
          chain.tail.y -= 1;
          append_to_some_chain = true;
          /* break; */
        } 
        if (chain.tail.y <= offset) {
          continue;
        }
        auto vert_dist = base_index - chain.tail.x;
        auto hori_dist = chain.tail.y - offset;
        auto& short_dist = vert_dist;
        auto& long_dist = hori_dist;
        if (short_dist > long_dist) { std::swap(short_dist, long_dist); }
        if (shift_factor * short_dist > long_dist) {
          /* account the last segment and append this pivot point as new segment to the chain. */
          chain.score += -std::log(std::pow(0.85, 2 * (k + chain.cur_seg_len)));
          chain.cur_seg_len = 0;
          chain.tail = Dot(base_index, offset);
          append_to_some_chain = true;
          continue;
          /* break; */
        }
      }
      if (!append_to_some_chain) {
        /* start a new chain if nowhere to append */
        reference_reads_scores.reverse_chains[read_index_rrs].push_back(Chain(Dot(base_index, offset)));
      }
    }
  };

  void scoring(
      InputParams& params, 
      AnchorsMap& anchors_map,
      uint32_t query_read_index,
      ReferenceReadsScores& reference_reads_scores,
      std::vector<biovoltron::FastaRecord<true>>& reads,
      std::vector<uint32_t>& ovg_cur_query_to_references) {
    auto shift_factor = 1.5f;
    auto& query_read = reads[query_read_index];
    auto k = params.kmer_size;
    const auto& anchors = anchors_map.get_anchors();
    const auto& blocked_loci = anchors_map.get_blocked_loci();
    const auto& blocked_query_loci = blocked_loci[query_read_index];

    for (uint32_t base_index = 0; base_index < query_read.seq.size() - k + 1; base_index++) {
      if (blocked_query_loci[base_index]) continue;
      /* if (blocked_query_loci.at( base_index + k )) { */
      /*   base_index += k; */
      /*   /1* TODO *1/ */
      /*   /1* auto additional_skip_dist = uint8_t(); *1/ */
      /*   continue; */
      /* } */
      auto seed_fw = biovoltron::istring_view(query_read.seq.c_str() + base_index, k);
      auto seed_rc = biovoltron::Codec::rev_comp(seed_fw);
      auto seed_fw_enc = biovoltron::Codec::hash(seed_fw);
      auto seed_rc_enc = biovoltron::Codec::hash(seed_rc);
      if (anchors.contains(seed_fw_enc)) {
        scoring_basewise<Strand::FORWARD>(reference_reads_scores, query_read_index, shift_factor, k, base_index, anchors.at(seed_fw_enc));
      }
      if (anchors.contains(seed_rc_enc)) {
        scoring_basewise<Strand::REVERSE>(reference_reads_scores, query_read_index, shift_factor, k, base_index, anchors.at(seed_rc_enc));
      }
    }
    for (uint32_t ref_read_index = 0; ref_read_index < reference_reads_scores.forward_chains.size(); ref_read_index++) {
      auto reported = false;
      for (auto& chain : reference_reads_scores.forward_chains[ref_read_index]) {
        if (std::min(chain.tail.x - chain.head.x, chain.tail.y - chain.head.y) < 200) {
          continue;
        }
        chain.score += -std::log(std::pow(0.85, 2 * (k + chain.cur_seg_len)));
        if (chain.score > params.thres_bin_score) {
          auto ref_read_id = ref_read_index + query_read_index + 1;
          ovg_cur_query_to_references.push_back(ref_read_id);
          break;
        }
      }
      if (reported) continue;
      for (auto& chain : reference_reads_scores.reverse_chains[ref_read_index]) {
        if (std::min(chain.tail.x - chain.head.x, chain.head.y - chain.tail.y) < 200) {
          continue;
        }
        chain.score += -std::log(std::pow(0.85, 2 * (k + chain.cur_seg_len)));
        if (chain.score > params.thres_bin_score) {
          auto ref_read_id = ref_read_index + query_read_index + 1;
          ovg_cur_query_to_references.push_back(ref_read_id);
          break;
        }
      }
    }
  }

 public:
  void detection (InputParams& params, AnchorsMap& anchors_map, std::vector<biovoltron::FastaRecord<true>>& reads) {
    Timer t("detection start");
    auto ovg = std::unordered_map<uint32_t, std::vector<uint32_t>>();
    /* auto& blocked_loci = anchors_map.get_blocked_loci(); */
    /* auto& anchors = anchors_map.get_anchors(); */
#pragma omp parallel for schedule(dynamic, 1)
    for (uint32_t query_read_index = 0; query_read_index < reads.size() - 1; query_read_index++) {
      auto ovg_cur_query_to_references = std::vector<uint32_t>();
      auto reference_reads_size = reads.size() - query_read_index - 1;
      auto reference_reads_scores = ReferenceReadsScores(reference_reads_size);
      scoring(params, anchors_map, query_read_index, reference_reads_scores, reads, ovg_cur_query_to_references);
#pragma omp critical
      {
        ovg.insert({query_read_index, ovg_cur_query_to_references});
      }
    }
    t.end();
    t.start("write ovg");
    auto fout = std::ofstream{std::filesystem::path(PROJECT_SOURCE_DIR) / "output/res.ovl"};
    assert(fout.is_open());
    for (const auto& [r1, r2s] : ovg) {
      for (const auto& r2 : r2s) {
        fout << r1 << ',';
        fout << r2 << '\n';
      }
    }
    t.end();
  }
};
