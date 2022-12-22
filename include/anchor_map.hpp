#pragma once
#include <algorithm>
#include <filesystem>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>

#include <options.hpp>
#include <utils.hpp>

#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/algo/align/exact_match/fm_index.hpp>


/* struct Read_offsets { */
/*   Read_offsets (uint32_t read_) : read(read_) {} */
/*   uint32_t read; */
/*   std::vector<uint32_t> offsets; */
/* }; */

/* using Anchor = std::vector<Read_offsets>; */

struct ReadOffset {
  ReadOffset(uint32_t offset_) : offset(offset_) {}
  uint32_t read;
  uint32_t offset;
};

using Anchor = std::vector<ReadOffset>;

class AnchorsMap {
 private:
  std::unordered_map<uint64_t, Anchor*> anchors;
  std::vector<std::vector<bool>> blocked_loci;

   void set_concat_reads (const InputParams& params, const std::vector<biovoltron::FastaRecord<true>>& reads, biovoltron::istring& concat_reads) {
    for (auto r : reads) {
      concat_reads += r.seq;
    }
    /* kind of compensation for some design flaw of this fm-index as I know */
    for (int i = 0; i < params.kmer_size; i++) {
      concat_reads.push_back(0);
    }
  }

  biovoltron::FMIndex<> get_fmi (const InputParams& params, const std::vector<biovoltron::FastaRecord<true>>& reads, biovoltron::istring& concat_reads) {
    Timer t("get_fmi");
    set_concat_reads(params, reads, concat_reads);
    auto fmi = biovoltron::FMIndex<>();
    auto fmi_file = params.fmi_dir / (std::string("batch_0") + ".fmi");
    if (params.fmi_src == "build") {
      if (!std::filesystem::exists(params.fmi_dir)) {
        std::filesystem::create_directory(params.fmi_dir);
      }
      fmi.build(concat_reads);
      auto fout = std::ofstream(fmi_file, std::ios::binary);
      assert(fout.is_open());
      fmi.save(fout);
    } else if (params.fmi_src == "load") {
      auto fin = std::ifstream(fmi_file, std::ios::binary);
      fmi.load(fin);
    } else {
      std::cout << "param error" << std::endl;
      exit(0);
    }
    t.end();
    return fmi;
  }

  std::vector<uint32_t> get_reads_start_offset_on_the_concatnated (const std::vector<biovoltron::FastaRecord<true>>& reads) {
    auto res = std::vector<uint32_t>{0};
    for (const auto& r : reads) {
      res.push_back(res.back() + r.seq.length());
    }
    return res;
  }

  std::unordered_set<uint64_t> get_encoded_seeds (const std::unordered_set<std::string>& blocked_seeds) {
    Timer t("encode blocked seeds");
    auto res = std::unordered_set<uint64_t>();
    for (const auto& s : blocked_seeds) {
      auto istr = biovoltron::Codec::to_istring(s);
      auto enc = biovoltron::Codec::hash(istr);
      res.insert(enc);
    }
    t.end();
    return res;
  }

 public:
  void build (const InputParams& params, const std::vector<biovoltron::FastaRecord<true>>& reads, const std::unordered_set<std::string>& blocked_seeds) {
    auto concat_reads = biovoltron::istring();
    auto fmi = get_fmi(params, reads, concat_reads);
    auto reads_start_offsets = get_reads_start_offset_on_the_concatnated(reads);
    auto length_concat_reads = reads_start_offsets.back();
    log(length_concat_reads, "length_concat_reads, ");
    auto short_seed_size = uint16_t(12);
    auto short_seeds_amount = uint32_t(std::pow(4, short_seed_size));
    auto blocked_loci_on_concat_reads = std::vector<bool>(length_concat_reads, false);
    auto checked_loci_on_concat_reads = std::vector<bool>(length_concat_reads, false);
    auto blocked_seeds_enc = get_encoded_seeds(blocked_seeds);

    using seed_enc_anchor = std::pair<uint64_t, Anchor*>;
    auto short_seed_wrt_anchor_list = std::vector<std::vector<seed_enc_anchor>*>(short_seeds_amount);

    Timer t("parallelly gather seeds");
#pragma omp parallel for
    for (uint32_t short_seed_enc = 0; short_seed_enc < short_seeds_amount; short_seed_enc++) {
      auto sseed = biovoltron::Codec::rhash(short_seed_enc, short_seed_size);
      auto [beg, end, ofsts] = fmi.get_range(sseed, 0);
      auto short_seed_offsets = fmi.get_offsets(beg, end);
      auto& space = params.dist_hori_filter;
      auto sseed_hit_amount = end - beg;
      if (sseed_hit_amount < 2) {
        continue;
      }
      auto level_seed_offsets_entry = new std::vector<seed_enc_anchor>();
      short_seed_wrt_anchor_list[short_seed_enc] = level_seed_offsets_entry;
      for (const auto& sseed_offset : short_seed_offsets) {
        if (checked_loci_on_concat_reads[sseed_offset]) {
          continue;
        }
        auto seed = biovoltron::istring_view(concat_reads.c_str() + sseed_offset, params.kmer_size);
        auto seed_enc = biovoltron::Codec::hash(seed);
        auto [b, e, o] = fmi.get_range(seed, 0);
        auto seed_hit_amount = e - b;
        /* should be a conservative estimation */
        /* auto maxSeedHitNum = reads.size(); */
        /* TODO */
        /* eliminate this constraint, use the previous line. */
        auto maxSeedHitNum = 100;
        if (seed_hit_amount < 2 || seed_hit_amount > maxSeedHitNum) {
          continue;
        }
        auto seed_offsets = fmi.get_offsets(b, e);
        if (blocked_seeds_enc.contains(seed_enc)) {
          for (const auto& seed_offset : seed_offsets) {
            checked_loci_on_concat_reads[seed_offset] = true;
            blocked_loci_on_concat_reads[seed_offset] = true;
          }
          continue;
        }
        auto anchor = new Anchor();
        short_seed_wrt_anchor_list[short_seed_enc]->push_back(seed_enc_anchor(seed_enc, std::move(anchor)));
        for (const auto& seed_offset : seed_offsets) {
          checked_loci_on_concat_reads[seed_offset] = true;
        }
        std::sort(seed_offsets.begin(), seed_offsets.end());
        /* see next only policy */
        /* auto& cur_anchor = short_seed_wrt_anchor_list[short_seed_enc]->back().second; */
        /* for (auto it = seed_offsets.begin(), next_it = std::next(it); */
        /*     next_it != seed_offsets.end(); */
        /*     it = next_it, next_it = std::next(next_it)) { */
        /*   if (*it + space > *next_it) { */
        /*     cur_anchor->push_back(ReadOffset(*it)); */
        /*   } else { */
        /*     /1* block the two and jump additional one iteration *1/ */
        /*     blocked_loci_on_concat_reads[*it] = true; */
        /*     blocked_loci_on_concat_reads[*next_it] = true; */
        /*     it = next_it; */
        /*     next_it = std::next(next_it); */
        /*     if (next_it == seed_offsets.end()) break; */
        /*   } */
        /* } */
        /* see next only policy */
        /* see prev and next policy */
        auto& cur_anchor = short_seed_wrt_anchor_list[short_seed_enc]->back().second;
        if (seed_offsets.size() == 2) {
          if (seed_offsets[0] + space < seed_offsets[1]) {
            cur_anchor->push_back(ReadOffset(seed_offsets[0]));
            cur_anchor->push_back(ReadOffset(seed_offsets[1]));
          } else {
            blocked_loci_on_concat_reads[seed_offsets[0]] = true;
            blocked_loci_on_concat_reads[seed_offsets[1]] = true;
          }
        } else {
          if (seed_offsets[0] + space < seed_offsets[1]) {
            cur_anchor->push_back(ReadOffset(seed_offsets[0]));
          } else {
            blocked_loci_on_concat_reads[seed_offsets[0]] = true;
          }
          for (int i = 1; i < seed_offsets.size() - 1; i++) {
            if (seed_offsets[i - 1] + space < seed_offsets[i] && 
                seed_offsets[i] + space < seed_offsets[i + 1]) {
              cur_anchor->push_back(seed_offsets[i]);
            } else {
              blocked_loci_on_concat_reads[seed_offsets[i]] = true;
            }
          }
          auto N = seed_offsets.size();
          if (seed_offsets[N - 2] + space < seed_offsets[N - 1]) {
            cur_anchor->push_back(seed_offsets[N - 1]);
          } else {
            blocked_loci_on_concat_reads[N - 1] = true;
          }
        }
        /* see prev and next policy */
      }
    }
    t.end();

    t.start("delimit anchors by reads");
    auto asdf = std::atomic<int>();
#pragma omp parallel for
    for (uint32_t short_seed_enc = 0; short_seed_enc < short_seeds_amount; short_seed_enc++) {
      /* skip invalid entries */
      if (short_seed_wrt_anchor_list[short_seed_enc] == nullptr) { continue; }

      for (const auto& [seed_enc, anchor] : *short_seed_wrt_anchor_list[short_seed_enc]) {
        auto it_read_start_loci = reads_start_offsets.begin();
        for (auto& read_offset : *anchor) {
          if (read_offset.offset >= *it_read_start_loci) {
            it_read_start_loci = std::upper_bound(it_read_start_loci, reads_start_offsets.end(), read_offset.offset);
            read_offset.read = std::distance(reads_start_offsets.begin(), it_read_start_loci) - 1;
            assert(read_offset.offset >= *std::prev(it_read_start_loci));
            read_offset.offset = read_offset.offset - *std::prev(it_read_start_loci);
          } else {
            read_offset.read = std::distance(reads_start_offsets.begin(), it_read_start_loci) - 1;
            assert(read_offset.offset >= *std::prev(it_read_start_loci));
            read_offset.offset = read_offset.offset - *std::prev(it_read_start_loci);
          }
        }
      }
    }
    t.end();

    t.start("merge into one");
    /* merge into one mapping table */
    for (uint32_t short_seed_enc = 0; short_seed_enc < short_seeds_amount; short_seed_enc++) {
      if (!short_seed_wrt_anchor_list[short_seed_enc]) {
        continue;
      }
      for (auto& [seed, anchor] : *short_seed_wrt_anchor_list[short_seed_enc]) {
        if (anchor->size() == 0) continue;
        this->anchors.insert({seed, anchor});
      }
    }
    log(this->anchors.size(), "this->anchors.size()");
    t.end();

    /* handle blocked from 1d to 2d */
    blocked_loci.resize(reads.size());
    auto it_rso = reads_start_offsets.begin();
    auto next_it_rso = std::next(it_rso);
    for (auto& r : blocked_loci) {
      r = std::vector<bool>(blocked_loci_on_concat_reads.begin() + *it_rso, blocked_loci_on_concat_reads.begin() + *next_it_rso);
      it_rso = next_it_rso;
      next_it_rso = std::next(next_it_rso);
    }
  }

  const std::unordered_map<uint64_t, Anchor*>& get_anchors () {
    return this->anchors;
  }

  const std::vector<std::vector<bool>>& get_blocked_loci () {
    return this->blocked_loci;
  }
};

