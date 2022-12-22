#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/algo/align/exact_match/fm_index.hpp>
#include <cmath>
#include <ranges>
#include <algorithm>

#include <utils.hpp>
#include <options.hpp>
#include <hpc_marmot.hpp>
#include <old_detect.hpp>

class Detection {
 public:
  std::ifstream fin;  // read in fasta reads file
  std::ofstream fout; // overlap pairs results
  InputParams params;
  struct Pivot_anchor {
    uint32_t index;
    uint32_t read;
    Pivot_anchor (uint32_t _index, uint32_t _read)
      : index(_index), read(_read) {}
  };

  struct Anchor {
    std::vector<uint32_t> offsets;
    std::vector<Pivot_anchor> pivot_anchors;
    Anchor() {
      offsets = std::vector<uint32_t>{};
      pivot_anchors = std::vector<Pivot_anchor>{};
    }
  };

  struct Dot_plot_meta {
    uint32_t id_last_appended_chain = UINT32_MAX;
    float last_chain_score;
    uint32_t times_switch_chain = 0;
  };

  Detection(InputParams& _params) 
    : params(_params)
  {
    fin = std::ifstream(params.reads_fasta);
    std::filesystem::create_directory(params.output_dir);
    fout = std::ofstream{params.output_dir / "ovlp.ovl"};
    assert(fout);
  }
  // thisexp
  std::unordered_map<uint32_t, std::vector<uint32_t>> ovg_gt;
  std::unordered_map<uint32_t, std::vector<uint32_t>> ovg_control;
  void set_ovg(std::filesystem::path ovg_gt_fpath, std::unordered_map<uint32_t, std::vector<uint32_t>>& ovg) {
    auto fin = std::ifstream{ovg_gt_fpath};
    assert(fin);
    auto str = std::string();
    auto acc = size_t(0);
    while (std::getline(fin, str, '\t')) {
      auto r1 = std::stoul(str);
      ovg[r1] = std::vector<uint32_t>();
      std::getline(fin, str);
      auto ss = std::stringstream(str);
      while (std::getline(ss, str, '\t')) {
        auto r2 = std::stoul(str);
        ovg[r1].emplace_back(std::move(r2));
      }
      acc += ovg[r1].size();
    }
    std::cout << "ovlp size()" << '\t';
    std::cout << acc << '\n';
  }

  biovoltron::istring get_hpc(const biovoltron::istring_view& ori_seq) {
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

  biovoltron::istring get_hpc(biovoltron::istring_view ori_seq, uint32_t l) {
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

    auto c = uint32_t(0);
    /* auto assigned = false; */
    for (auto i = uint32_t{1}; i < ori_seq.size(); i++) {
      constexpr int N = 4;
      if (ori_seq[i] == N) {
        // skip
      } else {
        if (ori_seq[i - 1] != ori_seq[i]) {
          if (c < l) {
            for (int i = 0; i < c; i++) {
              seq += seq.back();
            }
          }
          seq += ori_seq[i];
          c = 0;
        } else {
          c++;
        }
      }
    }
    return seq;
  }

  auto get_fmi(biovoltron::istring& ref) {
    static auto num_batch = uint16_t(0);
    auto fmi = biovoltron::FMIndex<>{};
    if (params.fmi_src == "build") {
      // check if fmi src dir exists
      if (!std::filesystem::exists(params.fmi_dir)) {
        std::filesystem::create_directory(params.fmi_dir);
      }
      // build and save
      // append A with length kmer_size, kind of unsolved bug for fm-index
      for (int i = 0; i < params.kmer_size; i++) {
        ref.push_back(0);
      }
      // ================================== CAUTION ==================================
      // Do TGS long reads has Ns? Comment this line for now, if error may due to this. 
      // fm-index only support 'ACTG', so we need to change 'N' to 'ACGT'
      // std::ranges::transform(ref, ref.begin(), [](auto& c) { return c % 4; });
      // ================================== CAUTION ==================================

      fmi.build(ref);
      auto fout = std::ofstream(params.fmi_dir / ("batch" + std::to_string(num_batch) + ".fmi"), std::ios::binary);
      fmi.save(fout);
    } else if (params.fmi_src == "load") {
      // check if fmi src dir exists
      if (!std::filesystem::exists(params.fmi_dir)) {
        std::cerr << "dir for fm-index files not exists! exit.\n";
        exit(0);
      }
      // load
      auto fmi_thisBatch = params.fmi_dir / ("batch" + std::to_string(num_batch) + ".fmi");
      if (!std::filesystem::exists(fmi_thisBatch)) {
        std::cerr << "fmi file for this batch not exists! exit.\n";
        exit(0);
      }
      auto fin = std::ifstream(fmi_thisBatch, std::ios::binary);
      assert(fin);
      fmi.load(fin);
    } else {
      std::cerr << "fmi_src option error! exit.\n";
      exit(0);
    }
    return fmi;
  }

  auto get_anchorIntervalsOnReads(
      biovoltron::FMIndex<>& fmi, 
      std::vector<biovoltron::FastaRecord<true>>& treads,
      /* std::vector<Anchor*>& anchors, */
      std::unordered_map<size_t, std::unique_ptr<Anchor>>& anchors,
      std::vector<std::vector<bool>>& valid_base_on_querys,
      biovoltron::istring& ref) {
    Timer timer;

    auto unordered_map_policy = [&] () {
      timer.start("long seed");
      // build intervals_template_reads
      auto acc_lengthOfTemplateReads = std::vector<uint32_t>{0};
      for (const auto& tread : treads) {
        acc_lengthOfTemplateReads.push_back(acc_lengthOfTemplateReads.back() + tread.seq.length());
      }
      auto get_templated_reads_total_length = [&treads] () {
        auto res = uint32_t(0);
        for (const auto& r : treads) {
          res += r.seq.size();
        }
        return res;
      };
      auto tlen = get_templated_reads_total_length();
      auto short_seed_size = uint16_t(9);
      auto short_seeds_amount = size_t(std::pow(4, short_seed_size));
      log(short_seeds_amount, "short_seeds_amount");
      /* auto actual_short_seed_amount = size_t(0); */
      /* pair of [ encoded long seed, offsets on tempalte reads ] */
      using long_anchor = std::pair<size_t, std::unique_ptr<Anchor>>;
      /* catgorize long seeds by their prefix */
      auto prefix_cat_anchors = std::vector<std::vector<long_anchor>>(short_seeds_amount);
      auto checked_loci = std::vector<bool>(tlen, false);
#pragma omp parallel for
      for (size_t seed_enc = 0; seed_enc < short_seeds_amount; seed_enc++) {
        auto sseed = biovoltron::Codec::rhash(seed_enc, short_seed_size);
        auto [beg, end, ofsts] = fmi.get_range(sseed, 0);
        auto short_seed_offsets = fmi.get_offsets(beg, end);
        auto checked_long_seed = std::unordered_set<int64_t>();
        for (const auto& offset : short_seed_offsets) {
          if (!checked_loci[offset]) {
            auto lseed = biovoltron::istring_view(ref.c_str() + offset, params.kmer_size);
            auto elseed = biovoltron::Codec::hash(lseed);
            auto lseed_exists = (checked_long_seed.find(elseed) != checked_long_seed.end());
            if (lseed_exists) {
              continue;
            } else {
              checked_long_seed.insert(elseed);
            }
            auto [b, e, o] = fmi.get_range(lseed, 0);
            if (e - b > 1 && e - b < params.maxNumSeedAnchors) {
              auto ofst_lsd = fmi.get_offsets(b, e);
              /* auto acr = new Anchor(); */
              auto acr = std::make_unique<Anchor>();
              std::sort(ofst_lsd.begin(), ofst_lsd.end());
              // =========== see prev and next policy ===========
              if (ofst_lsd[0] + params.dist_hori_filter < ofst_lsd[1]) {
                checked_loci[ofst_lsd[0]] = true;
                acr->offsets.emplace_back(std::move(ofst_lsd[0]));
              }
              for (int i = 1; i < ofst_lsd.size() - 1; i++) {
                if (ofst_lsd[i - 1] + params.dist_hori_filter < ofst_lsd[i] && 
                    ofst_lsd[i] + params.dist_hori_filter < ofst_lsd[i + 1]) {
                  checked_loci[ofst_lsd[i]] = true;
                  acr->offsets.emplace_back(std::move(ofst_lsd[i]));
                } else {
                }
              }
              if (ofst_lsd[ofst_lsd.size() - 2] + params.dist_hori_filter < ofst_lsd.back()) {
                checked_loci[ofst_lsd.back()] = true;
                acr->offsets.emplace_back(std::move(ofst_lsd.back()));
              }
              // =========== see prev and next policy ===========
              auto lseed_enc = biovoltron::Codec::hash(lseed);
              prefix_cat_anchors[seed_enc].emplace_back(std::move(long_anchor(lseed_enc, std::move(acr))));
            }
          }
        }
      }
      timer.end();
      /* ==================== */
#define comt false
#if comt
        for (size_t seed = 1; seed < anchors.size(); seed++) {
          // if (anchors[seed]->offsets.size() > params.minNumSeedAnchors) {
          if (anchors[seed]) {
            auto it_read = acc_lengthOfTemplateReads.begin();
            auto it_acr = anchors[seed]->offsets.begin();
            auto searcher = it_read;
            while (it_acr != anchors[seed]->offsets.end()) {
              searcher = std::upper_bound(it_read, acc_lengthOfTemplateReads.end(), *it_acr);
              if (searcher != acc_lengthOfTemplateReads.end()) {
                anchors[seed]->pivot_anchors.emplace_back(std::move(Pivot_anchor(std::distance(anchors[seed]->offsets.begin(), it_acr), std::distance(acc_lengthOfTemplateReads.begin(), searcher) - 1)));
                it_read = searcher;
                if (it_read == acc_lengthOfTemplateReads.end()) {
                  break;
                }
                while (it_acr != anchors[seed]->offsets.end() && *it_acr < *it_read) {
                  it_acr += 1;
                }
              } else {
                anchors[seed]->pivot_anchors.emplace_back(std::move(Pivot_anchor (std::distance(anchors[seed]->offsets.begin(), it_acr), std::distance(acc_lengthOfTemplateReads.begin(), searcher) - 1)));
                break;
              }
            }
          }
        }
#endif
      /* ==================== */
      timer.start("delimit");
#pragma omp parallel for
      for (size_t seed_enc = 0; seed_enc < short_seeds_amount; seed_enc++) {
        for (auto& [lseed_enc, acr] : prefix_cat_anchors[seed_enc]) {
          auto it_read = acc_lengthOfTemplateReads.begin();
          auto it_acr = acr->offsets.begin();
          auto searcher = it_read;
          while (it_acr != acr->offsets.end()) {
            searcher = std::upper_bound(it_read, acc_lengthOfTemplateReads.end(), *it_acr);
            if (searcher != acc_lengthOfTemplateReads.end()) {
              acr->pivot_anchors.emplace_back(Pivot_anchor(std::distance(acr->offsets.begin(), it_acr), std::distance(acc_lengthOfTemplateReads.begin(), searcher) - 1));
              it_read = searcher;
              if (it_read == acc_lengthOfTemplateReads.end()) {
                break;
              }
              while (it_acr != acr->offsets.end() && *it_acr < *it_read) {
                it_acr += 1;
              }
            } else {
              acr->pivot_anchors.emplace_back(Pivot_anchor(std::distance(acr->offsets.begin(), it_acr), std::distance(acc_lengthOfTemplateReads.begin(), searcher) - 1));
              break;
            }
          }
        }
      }
      timer.end();
      timer.start("merge");
      for (size_t seed_enc = 0; seed_enc < short_seeds_amount; seed_enc++) {
        /* for (auto& [lseed_enc, acr] : prefix_cat_anchors[seed_enc]) { */
        for (auto& lseed_enc_acr : prefix_cat_anchors[seed_enc]) {
          anchors.emplace(std::move(lseed_enc_acr));
        }
      }
      timer.end();
      /* mark valid_base_on_querys using checked_loci */ 
      log(ref.size(), "ref.size()");
      log(acc_lengthOfTemplateReads.back(), "acc_lengthOfTemplateReads.back()");
      valid_base_on_querys.resize(treads.size());
      for (int i = 0; i < treads.size(); i++) {
        valid_base_on_querys[i] = std::vector<bool>(treads[i].seq.size(), false);
        for (int j = 0; j < treads[i].seq.size() - params.kmer_size + 1; j++) {
          valid_base_on_querys[i][j] = checked_loci[acc_lengthOfTemplateReads[i] + j];
        }
      }
    };
    unordered_map_policy();

    auto kmer_size = this->params.kmer_size;

    /* auto save_anchor_offsets_table = [&anchors, &valid_base_on_querys, &kmer_size] () { */
    /*   Timer timer; */
    /*   timer.start("save_anchor_offsets_table"); */
    /*   /1* std::unordered_map<size_t, Anchor*>& anchors, *1/ */
    /*   /1* std::vector<std::vector<bool>>& valid_base_on_querys, *1/ */
    /*   auto fout = std::ofstream{"../midFiles/devnet/anchor_offsets_table_k" + std::to_string(kmer_size) + ".txt"}; */
    /*   for (const auto& anchor : anchors) { */
    /*     fout << anchor.first << '\n'; */
    /*     fout << anchor.second->pivot_anchors.size() << '\n'; */
    /*     for (const auto& pivot : anchor.second->pivot_anchors) { */
    /*       fout << pivot.index << '\t'; */
    /*       fout << pivot.read << '\t'; */
    /*     } */
    /*     fout << '\n'; */

    /*     fout << anchor.second->offsets.size() << '\n'; */
    /*     for (const auto& offset : anchor.second->offsets) { */
    /*       fout << offset << '\t'; */
    /*     } */
    /*     fout << '\n'; */
    /*   } */
    /*   fout.close(); */
    /*   timer.end(); */
    /*   timer.start("valid base build"); */
    /*   auto fout2 = std::ofstream{"../midFiles/devnet/valid_base_file_k" + std::to_string(kmer_size) + ".txt"}; */
    /*   for (const auto& read_vb : valid_base_on_querys) { */
    /*     fout2 << read_vb.size() << '\n'; */
    /*     for (const auto& validity : read_vb) { */
    /*       fout2 << validity << '\t'; */
    /*     } */
    /*     fout2 << '\n'; */
    /*   } */
    /*   timer.end(); */
    /* }; */
    /* save_anchor_offsets_table(); */

    /* auto load_anchor_offsets_table = [&anchors, &valid_base_on_querys, kmer_size] () { */
    /*   Timer timer; */
    /*   timer.start("load_anchor_offsets_table"); */
    /*   auto fin = std::ifstream{"../midFiles/devnet/anchor_offsets_table_k" + std::to_string(kmer_size) + ".txt"}; */
    /*   auto str = std::string(); */
    /*   while (std::getline(fin, str)) { */
    /*     auto seed = std::stoull(str); */
    /*     auto anchor_inst = std::make_unique<Anchor>(); */
    /*     std::getline(fin, str); */
    /*     for (uint32_t i = 0; i < std::stoul(str); i++) { */
    /*       std::getline(fin, str, '\t'); */
    /*       auto index = std::stoul(str); */
    /*       std::getline(fin, str, '\t'); */
    /*       auto read = std::stoul(str); */
    /*       anchor_inst->pivot_anchors.push_back(Pivot_anchor( index, read )); */
    /*     } */
    /*     std::getline(fin, str); */

    /*     std::getline(fin, str); */
    /*     for (uint32_t i = 0; i < std::stoul(str); i++) { */
    /*       std::getline(fin, str, '\t'); */
    /*       anchor_inst->offsets.push_back(std::stoul(str)); */
    /*     } */
    /*     std::getline(fin, str); */
    /*     anchors.insert({seed, std::move(anchor_inst)}); */
    /*   } */
    /*   timer.end(); */

    /*   timer.start("read_valid_base"); */
    /*   fin = std::ifstream{"../midFiles/devnet/valid_base_file_k" + std::to_string(kmer_size) + ".txt"}; */
    /*   auto inc = uint32_t(0); */
    /*   while (std::getline(fin, str)) { */
    /*     auto read_size = std::stoul(str); */
    /*     for (uint32_t i = 0; i < read_size; i++) { */
    /*       std::getline(fin, str, '\t'); */
    /*       valid_base_on_querys[inc][i] = str[0] - '0'; */
    /*     } */
    /*     std::getline(fin, str); */
    /*   } */
    /*   timer.end(); */
    /* }; */
    /* load_anchor_offsets_table(); */
  }

  void detection_self(
      /* std::vector<Anchor*> anchors, */ 
      std::unordered_map<size_t, std::unique_ptr<Anchor>>& anchors, 
      uint16_t templateBatchId, 
      std::vector<biovoltron::FastaRecord<true>>& reads,
      std::vector<std::vector<bool>>& valid_base_on_querys) {
    // this template batch number
    // treads[i] for i : [0, N]
    auto& kmer_size = params.kmer_size;
    auto ov_graph = std::vector<std::vector<size_t>>(reads.size());
    auto timer = Timer();
    InputParams& params = this->params;
    auto read_id_offset = params.batchSize_readsNum * templateBatchId;

    auto detection_simple_bin = [&] () {
      timer.start("detection");

      struct SeedChain {
        uint32_t id;
        std::pair<uint32_t, uint32_t> head_dot;
        /* Assume it's rare for a seed chain to be misleaded by inappropriate tail dot. */
        float acc_prob_score = 0;
        /* amount of continuous dots for the current group, -1 stands for no unscored point. */
        int cur_cont_dots_amt = 0;
        uint32_t total_cont_dots_amt = 0;
        uint32_t group_amount = 1;
        std::pair<uint32_t, uint32_t> tail_dot;
        SeedChain (uint32_t id_, const std::pair<uint32_t, uint32_t>& head_dot_)
          : id(id_), head_dot(head_dot_) {
            tail_dot = head_dot;
          }
      };

      /* local parameters */
      /* [max_shift / 10] == [long offset / short offset] */
      auto max_shift = 15;
      /* auto distance_to_cut_seed_chain = 500; */
      /* auto mut_cent = uint8_t(params.kmer_size / 2); */
      /* auto mut_max_shift = 12; */

      /* gapped seed related */
      auto get_cleaner = [&] (uint8_t k, uint8_t pos) {
        auto res = size_t(0);
        auto str = biovoltron::istring();
        for (int i = 0; i < k - pos; i++) {
          str.push_back(3);
        }
        for (int i = 0; i < pos; i++) {
          str.push_back(0);
        }
        res = biovoltron::Codec::hash(str);
        return res;
      };
      auto dist_gap = 3;
      /* auto get_gap_operator = [&get_cleaner, &kmer_size, &dist_gap] (const size_t& eseed) { */
      /*   auto cleaners = std::vector<size_t>{}; */
      /*   auto cent_pos = kmer_size / 2; */
      /*   for (int i = cent_pos - dist_gap + 1; i < cent_pos + dist_gap; i++) { */
      /*     cleaners.push_back(get_cleaner(kmer_size, i)); */
      /*   } */
      /*   auto res = std::vector<size_t>(); */
      /*   for (const auto& cleaner : cleaners) { */
      /*     auto left = eseed & (cleaner << 2); */
      /*     auto right = (eseed & ~cleaner) << 2; */
      /*     auto concat = left | right; */
      /*     res.emplace_back(std::move(concat)); */
      /*     for (int i = 1; i < 4; i++) { */
      /*       res.push_back(concat | i); */
      /*     } */
      /*   } */
      /*   return res; */
      /* }; */

      /* snpped seed related */
      auto i_sigma = std::vector<uint32_t>{0,1,2,3};
      auto mut_pos = std::vector<uint8_t>{7, 8, 9};
      auto get_mutation = [&i_sigma, &mut_pos] () {
        auto res = std::vector<size_t>();
        for (auto& p : mut_pos) {
          for (auto& i : i_sigma) {
            /* two bits as a unit for biovoltron istring residue */
            res.push_back(i << (2*p));
          }
        }
        return res;
      };
      auto mut_poses = get_mutation();


      auto score_single_hit = -std::log(std::pow(0.85, 2 * (params.kmer_size)));

      auto score_bins = [&anchors, &params, &max_shift, &score_single_hit] (
          const uint32_t& it_read,
          const uint32_t& it_base,
          const size_t& seed_enc, 
          std::vector<std::vector<SeedChain>>& seed_chains,
          std::vector<Dot_plot_meta>& dot_plots,
          bool forward, 
          bool mut_mode) {

        /* ckech if this seed exists in anchors table */
        if (!anchors.contains(seed_enc)) {
          return;
        }
        auto& anchor = anchors.at(seed_enc);
        auto& pivots = anchor->pivot_anchors;
        if (pivots.size() < 2) {
          return;
        }

        /* search the first anchor beyond the start of the current query read. */
        auto upper = std::upper_bound(
            pivots.begin(), pivots.end(), 
            it_read, [] (size_t read_id, const Pivot_anchor& item) { return read_id < item.read; });
        auto init_i_pv = size_t(std::distance(pivots.begin(), upper));

        /* iterate on pivot anchors w.r.t. it's corresponding template reads */
        for (size_t i_pv = init_i_pv; i_pv < pivots.size() - 1; i_pv++) {
          auto& cur_pivot = pivots[i_pv];
          auto& id_tread = cur_pivot.read;
          auto it_tread = id_tread - 1 - it_read;
          auto& dot_plot = dot_plots[it_tread];
          /* if (dot_plot.times_switch_chain > 5) { */
          /*   continue; */
          /* } */
          auto& y = it_base;
#define THRESHOLD_TIMES_SWITCH_CHAIN 30
          /* if (dot_plot.times_switch_chain > THRESHOLD_TIMES_SWITCH_CHAIN) { */
          /*   return; */
          /* } */

          /* iterate on dots in the same tread for this anchor */
          for (size_t it_index = pivots[i_pv].index; it_index < pivots[i_pv + 1].index; it_index++) {
            auto& x = anchor->offsets[it_index];
            auto append_to_some_chain = false;

            /* if no rectangle-scattered dots, start diagonal colinear clustering */
            for (auto& seed_chain : seed_chains[it_tread]) {
              if (seed_chain.tail_dot.first == y) {
                /* possible consistent y value when mutation seed is considered, which is repeatitive iteration on seed hits for this y-axis. */
                continue;
              }
              auto colinear_possible = forward ? (seed_chain.tail_dot.second < x) : (seed_chain.tail_dot.second > x);
              if (colinear_possible) {
                auto hori_diff = forward ? (x - seed_chain.tail_dot.second) : (seed_chain.tail_dot.second - x);
                auto vert_diff = y - seed_chain.tail_dot.first;
                if (hori_diff == 1 && vert_diff == 1) {
                  seed_chain.cur_cont_dots_amt += 1;
                  seed_chain.total_cont_dots_amt += 1;
                  seed_chain.tail_dot = {y, x};
                  append_to_some_chain = true;
                  if (seed_chain.id != dot_plot.id_last_appended_chain) {
#define MID_THRES_CHAIN_SCORE 15
                    if (dot_plot.last_chain_score < MID_THRES_CHAIN_SCORE && seed_chain.acc_prob_score < MID_THRES_CHAIN_SCORE) {
                      dot_plot.times_switch_chain += 1;
                    }
                    dot_plot.id_last_appended_chain = seed_chain.id;
                    /* dot_plot.times_switch_chain++; */
                  }
                  /* break; */
                  continue;
                }
                /* consider the following hori_diff as short_diff, vert_diff as long_diff */
                if (hori_diff > vert_diff) {
                  std::swap(hori_diff, vert_diff);
                }
                if (hori_diff * max_shift / 10 > vert_diff) {
                  seed_chain.acc_prob_score += -std::log(std::pow(0.85, 2 * (params.kmer_size + seed_chain.cur_cont_dots_amt)));
                  seed_chain.cur_cont_dots_amt = 0;
                  seed_chain.group_amount += 1;
                  seed_chain.tail_dot = {y, x};
                  append_to_some_chain = true;
                  if (seed_chain.id != dot_plot.id_last_appended_chain) {
                    if (dot_plot.last_chain_score < MID_THRES_CHAIN_SCORE && seed_chain.acc_prob_score < MID_THRES_CHAIN_SCORE) {
                      dot_plot.times_switch_chain += 1;
                    }
                    dot_plot.id_last_appended_chain = seed_chain.id;
                    /* dot_plot.times_switch_chain++; */
                  }
                  dot_plot.last_chain_score = seed_chain.acc_prob_score;
                  /* break; */
                  continue;
                }
              }
            }
            /* if this dot has no cluster to append, initialize a new one on it's own */
            if (!append_to_some_chain) {
              auto id_new_chain = seed_chains[it_tread].size();
              seed_chains[it_tread].emplace_back(SeedChain(id_new_chain, std::pair<uint32_t, uint32_t>(y, x)));
              dot_plot.id_last_appended_chain = id_new_chain;
              if (dot_plot.last_chain_score < MID_THRES_CHAIN_SCORE) {
                dot_plot.times_switch_chain += 1;
              }
              dot_plot.last_chain_score = score_single_hit;
            }
          }
        }
      };

      /* parameters for scoring mechanism */
      /* auto prob_seed_hit = std::pow(0.85f, 2 * params.kmer_size); */
      /* auto score_prob_sole_anchor = -std::log(std::pow(0.85, 2 * (params.kmer_size))); */
      /* auto sum_k = uint32_t(0); */
      auto error_rate = 0.15f;
      auto dec_expect = 0.0f;
      for (uint32_t i = 1; i < params.kmer_size; i++) {
        dec_expect += std::pow(1 - error_rate, 2 * i);
      }

      /* block list of seeds from file */
      /* via modifying valid_base_on_querys */
      /* of format begin lines of [seed, biovoltron::Codec::hash(seed)] (e.g.) [AAAT 3] */
      /* auto mark_blocked_seeds_from_file = [&anchors, &valid_base_on_querys, &reads] (const std::filesystem::path& file) { */
      /*   Timer t; */
      /*   t.start("mark_blocked_seeds_from_file"); */
      /*   // build intervals_template_reads */
      /*   auto acc_lengthOfTemplateReads = std::vector<uint32_t>{0}; */
      /*   for (const auto& tread : reads) { */
      /*     acc_lengthOfTemplateReads.push_back(acc_lengthOfTemplateReads.back() + tread.seq.length()); */
      /*   } */
      /*   auto fin = std::ifstream{file}; */
      /*   assert(fin.is_open()); */
      /*   auto seed = ""s; */
      /*   auto seed_enc = uint64_t(0); */
      /*   auto to_be_blocked_seeds_enc = std::vector<uint64_t>(); */
      /*   auto inc = 0; */
      /*   while (fin >> seed >> seed_enc) { */
      /*     if (inc++ < 5) { */
      /*       std::cout << seed << ',' << seed_enc << std::endl; */
      /*     } */
      /*     to_be_blocked_seeds_enc.emplace_back(std::move(seed_enc)); */
      /*   } */
      /*   log("to_be_blocked_seeds_enc.size()", to_be_blocked_seeds_enc.size()); */
      /*   for (const auto& seed_enc : to_be_blocked_seeds_enc) { */
      /*     auto it_seed = anchors.find(seed_enc); */
      /*     if (it_seed != anchors.end()) { */
      /*       anchors.erase(it_seed); */
      /*     } */
      /*   } */
      /*   t.end(); */
      /*   return; */
      /* }; */
      /* mark_blocked_seeds_from_file("/mnt/es/ness/johnson/thesis/overlaps/overlap_ver4/tests/iterKmer/output/noises_seeds.txt"); */


      #pragma omp parallel for schedule (dynamic, 1)
      for (uint32_t it_read = 0; it_read < reads.size(); it_read++) {
        // std::cout << it_read << std::endl;
        auto seed_chains_fw = std::vector<std::vector<SeedChain>>(reads.size() - it_read - 1);
        auto seed_chains_rc = std::vector<std::vector<SeedChain>>(reads.size() - it_read - 1);
        auto dot_plots_fw = std::vector<Dot_plot_meta>(reads.size() - it_read - 1);
        auto dot_plots_rc = std::vector<Dot_plot_meta>(reads.size() - it_read - 1);
        for (uint32_t it_base = 0; it_base < reads[it_read].seq.size() - kmer_size + 1; it_base++) {
          if (valid_base_on_querys[it_read][it_base]) {
            auto seed_fw = biovoltron::istring_view(reads[it_read].seq.c_str() + it_base, kmer_size);
            auto seed_enc_fw = biovoltron::Codec::hash(seed_fw);
            score_bins(it_read, it_base, seed_enc_fw, seed_chains_fw, dot_plots_fw, true, false);

            auto seed_rc = biovoltron::Codec::rev_comp(seed_fw);
            auto seed_enc_rc = biovoltron::Codec::hash(seed_rc);
            score_bins(it_read, it_base, seed_enc_rc, seed_chains_rc, dot_plots_rc, false, false);
            
            /* auto gapped_seeds_fw = get_gap_operator(seed_enc_fw); */
            /* for (const auto& s : gapped_seeds_fw) { */
            /*   score_bins(it_read, it_base, s, seed_chains_fw, dot_plots_fw, true, true); */
            /* } */
            /* auto gapped_seeds_rc = get_gap_operator(seed_enc_rc); */
            /* for (const auto& s : gapped_seeds_rc) { */
            /*   score_bins(it_read, it_base, s, seed_chains_rc, dot_plots_rc, false, true); */
            /* } */
            /* for (const auto& m : mut_poses) { */
            /*   if ((seed_enc_fw ^ m) != seed_enc_fw) { */
            /*     score_bins(it_read, it_base, seed_enc_fw ^ m, seed_chains_fw, dot_plots_fw, true, true); */
            /*   } */
            /*   if ((seed_enc_rc ^ m) != seed_enc_rc) { */
            /*     score_bins(it_read, it_base, seed_enc_rc ^ m, seed_chains_rc, dot_plots_rc, false, false); */
            /*   } */
            /* } */

          }
        } // end iter base
            
        /* treads for this query read has finished bin-level scoring. */
        /* Iterate all treads, select bin with maximum score and judge overlap with thrshold on dot count. */
        /* ovg_thisQuery declaration for this query read */
        for (size_t it_tread = 0; it_tread < seed_chains_fw.size(); it_tread++) {
#define LOW_SCORE_THRESHOLD 30
          for (auto& seed_chain : seed_chains_fw[it_tread]) {
#define DOT_COUNT_POLICY true
#if DOT_COUNT_POLICY
            seed_chain.acc_prob_score += -std::log(std::pow(0.85, 2 * (params.kmer_size + seed_chain.cur_cont_dots_amt)));
            if (seed_chain.acc_prob_score > params.thres_bin_score) {
              ov_graph[it_read].push_back(it_tread + 1 + it_read);
              continue;
            }
#else
            seed_chain.acc_prob_score += -std::log(std::pow(0.85, 2 * (params.kmer_size + seed_chain.cur_cont_dots_amt)));
            auto short_L = seed_chain.tail_dot.first - seed_chain.head_dot.first;
            auto long_L = seed_chain.tail_dot.second - seed_chain.head_dot.second;
            if (short_L > long_L) {
              std::swap(short_L, long_L);
            }
            if (short_L < 200) {
              continue;
            }
            auto L_prom = short_L * params.kmer_size;
            auto mu = L_prom * prob_seed_hit - seed_chain.group_amount * dec_expect;
            if (seed_chain.acc_prob_score > params.thres_bin_score) {
              ov_graph[it_read].push_back(it_tread + 1 + it_read);
              continue;
            }
#endif
          }
          for (auto& seed_chain : seed_chains_rc[it_tread]) {
#if DOT_COUNT_POLICY
            seed_chain.acc_prob_score += -std::log(std::pow(0.85, 2 * (params.kmer_size + seed_chain.cur_cont_dots_amt)));
            if (seed_chain.acc_prob_score > params.thres_bin_score) {
              ov_graph[it_read].push_back(it_tread + 1 + it_read);
              continue;
            }
#else
            seed_chain.acc_prob_score += -std::log(std::pow(0.85, 2 * (params.kmer_size + seed_chain.cur_cont_dots_amt)));
            auto short_L = seed_chain.tail_dot.first - seed_chain.head_dot.first
            auto long_L = seed_chain.head_dot.second - seed_chain.tail_dot.second;
            if (short_L > long_L) {
              std::swap(short_L, long_L);
            }
            if (short_L < 200) {
              continue;
            }
            auto L_prom = short_L * params.kmer_size;
            auto mu = L_prom * prob_seed_hit - seed_chain.group_amount * dec_expect;
            if (seed_chain.acc_prob_score > params.thres_bin_score) {
              ov_graph[it_read].push_back(it_tread + 1 + it_read);
              continue;
            }
#endif
          }
        }
        // critical section to insert into global ovg
      } // end iter read
      timer.end();
    };
    detection_simple_bin();

    timer.start("write ovlp");
    auto it_read = 0;
    for (const auto& r2s : ov_graph) {
      for (const auto& r2 : r2s) {
        fout << read_id_offset + it_read << ',';
        fout << read_id_offset + r2 << '\n';
      }
      it_read++;
    }
    fout.flush();
    timer.end();
  } // end detection_self
  
  /* void detection_QtoT( */
  /*     std::vector<Pivot_anchor>& anchorIntervalsOnReads, */
  /*     uint16_t queryBatchId, */ 
  /*     uint16_t templateBatchId, */ 
  /*     std::vector<biovoltron::FastaRecord<true>>& reads, */
  /*     std::vector<std::vector<bool>>& valid_base_on_querys) { */
  /*   // query batch number, template batch number */
  /*   // query batch : B[i], in which i in [0 .. N - 1], template batch B[j], j in [i + 1 .. N] */
  /* } */

  void batchelization() {
    auto fin = std::ifstream{params.reads_fasta};
    assert(fin);
    auto& batchSize_readsNum = params.batchSize_readsNum;
    auto processBatchPair = [] (int b1, int b2) {
      std::cout << "processBatchPair " << b1 << ", " << b2 << '\n';
    };
    auto templateBatchId = uint16_t(0);
    auto cursor_treadEnd = fin.tellg();
    /* auto showReads = [] (auto& reads) { */
    /*   for (const auto& read : reads) { */
    /*     std::cout << read << '\n'; */
    /*   } */
    /* }; */
    while (true) { // iter read template batch
      auto iter_trec = 1; // Record iter for template reads, start from read for this batch as 1.
      /* auto anchors = std::vector<Anchor*>{}; */
      auto anchors = std::unordered_map<size_t, std::unique_ptr<Anchor>>{};
      //        r1 -------------------- r2 --------- r3 -------------- r4 ------
      // ATCTG        * (hit)   *            *                              * *
      // anchors::pivot_anchors [ATCTG] = { {0th*, r1}, {2th*, r2}, {3th*, r4} }
      //
      fin = std::ifstream{params.reads_fasta};
      { // from file read in the first batch, after self detection, ref, treads can be ignored for the following qread detection.
        // On the other hand, anchors::pivot_anchors should be retained during this template reads as batch.
        fin.seekg(cursor_treadEnd); // This should move cursor to the last time the curosr moved to the end number of template reads.
        auto ref = biovoltron::istring{};
        // read in the first batch.
        auto treads = std::vector<biovoltron::FastaRecord<true>>();
        for (auto& r : std::ranges::istream_view<biovoltron::FastaRecord<true>>(fin)) {
          // std::cout << r.name << '\n';
          iter_trec++;
          auto hpc_seq = get_hpc(r.seq);
          ref += hpc_seq;
          /* ref += r.seq; */
          r.seq = hpc_seq;
          treads.push_back(r);
          if (iter_trec > batchSize_readsNum) { break; }
        }
        auto valid_base_on_querys = std::vector<std::vector<bool>>();
        cursor_treadEnd = fin.tellg(); // Save the cursor position for template read number.
        {
          // build fm-index for this template batch, if not exists; otherwise, load existing fm-index file.
          // after this scope, fmi should be freed. Keep Anchor::pivot_anchors only.
          auto timer = Timer();
          timer.start("get fmi");
          auto fmi = get_fmi(ref);
          timer.end();
          // [Anchor::pivot_anchors, seedToOffsets] = get_anchorIntervalsOnReads(fmi, treads);
          timer.start("get_anchorIntervalsOnReads");
          get_anchorIntervalsOnReads(fmi, treads, anchors, valid_base_on_querys, ref);
          timer.end();
        }
        // processBatchPair
        processBatchPair(templateBatchId, templateBatchId);
        // detection_itPivot(
        detection_self(
            anchors,
            templateBatchId,
            treads,
            valid_base_on_querys);
        // showReads(treads);
      }
      //   0 1 2 3
      // 0 x x x x
      // 1   x x x
      // 2     x x
      // 3       x
      return;
      // seedsOffsets for this template of reads should be calculated and stored in this stage, 
      // for the following query batches to use, the scope of this seeds ranges offsets should end before next template batch.
      // fmi also.
      auto queryBatchId = templateBatchId + 1;
      while (true) {
        auto iter_qrec = 1;
        auto qreads = std::vector<biovoltron::FastaRecord<true>>();
        for (const auto& r : std::ranges::istream_view<biovoltron::FastaRecord<true>>(fin)) {
          iter_qrec++;
          qreads.emplace_back(std::move(r));
          if (iter_qrec > batchSize_readsNum) { break; }
        }
        if (iter_qrec == 1) { break; }
        // processBatchPair
        processBatchPair(templateBatchId, queryBatchId);
        // showReads(qreads);
        queryBatchId += 1;
      }
      if (iter_trec == 1) { break; }
      templateBatchId += 1;
      std::cout << " ------------ batch ------------ \n";
    }
  }
};

