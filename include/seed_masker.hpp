#pragma once
#include <iostream>
#include <fstream>
#include <filesystem>
#include <unordered_set>
#include <vector>
#include <cassert>
#include <ranges>
#include <algorithm>

#include <utils.hpp>

class SeedMasker {
 private:
  std::filesystem::path base_seed_rule_file;
  std::vector<std::string> base_seeds;
  std::unordered_set<std::string> derived_seeds;
  std::vector<std::string> sigma = {"A", "C", "G", "T"};
 public:
  std::unordered_set<std::string>& get_derived_seeds() {
    return derived_seeds;
  }
  void set_base_seed_rule_file(const std::filesystem::path& file) {
    if (!std::filesystem::exists(file)) {
      std::cerr << "K-mer block rule1 file not exists, exit. " << std::endl;
      exit(0);
    }
    this->base_seed_rule_file = file;
  }
  void set_short_seeds() {
    auto fin = std::ifstream{base_seed_rule_file};
    assert(fin.is_open());
    // base_seeds in file are endline delimited.
    auto sd = std::string();
    while (std::getline(fin, sd)) {
      if (sd.size() < 7) {
        std::cout << "seed size larger than 7, generates too much blocked"
          << "seeds, program exits." << std::endl;
        exit(0);
      }
      base_seeds.push_back(sd);
    }
  }
  
  /**
   * Extend all seeds by one residue in HPC mode.
   *
   * @param v list of seeds
   * @return v with the original seeds substituted with all possible seeds with whose length increased by 1 w.r.t. the original seed, the tail residue should not be the same as the original tail residue.
   */
  std::vector<std::string> append_by_one (std::vector<std::string>& v) {
    if (v.size() == 0) {
      for (int i = 0; i < 4; i++) {
        v.push_back(sigma[i]);
      }
      return v;
    } 
    auto segments = std::vector<std::string>();
    int N = v.size();
    for (int i = 0; i < N; i++) {
      auto& seed = v[i];
      auto orig_tail_c = seed.back();
      seed.push_back('A');
      for (int j = 0; j < 4; j++) {
        auto c = sigma[j];
        if (c[0] != orig_tail_c) {
          seed.back() = c[0];
          segments.push_back(seed);
        }
      }
    }
    return segments;
  };

  /**
   * Generate combinations of segments as sequence over `sigma`.
   * The generated segment are used to concatenated to the blanking region of `base_seed`.
   *
   * @param k Length of the segment.
   * @return list of segment with possible combinations.
   */
  std::vector<std::string> gen_combs_k (int k) {
    auto res = std::vector<std::string>();
    for (int i = 0; i < k; i++) {
      res = append_by_one(res);
    }
    return res;
  };

  /**
   * Generate list of seeds of combinations on the blanking segment of the `base_seed`.
   * Prefix refers to the segment on left and Postfix for the right.
   * combinations[prefix] + `base_seed` + combinations[posifix] would build a set of size 
   * || comb.[prefix] ||  X  || comb.[postfix] ||, as which would be appended to.
   *
   * @param prefix_len Length of prefix
   * @param posifix_len Length of posifix
   * @param base_seed Seed to be append.
   * @param res Result set of seeds.
   */
  void derive_for_fixed_base_seed_pos (
      uint8_t prefix_len, 
      uint8_t postfix_len,
      const std::string& base_seed,
      std::unordered_set<std::string>& res) {
    auto prefix_seed_set = gen_combs_k(prefix_len);
    auto postfix_seed_set = gen_combs_k(postfix_len);
    if (prefix_seed_set.size() == 0) { prefix_seed_set.push_back({""}); }
    if (postfix_seed_set.size() == 0) { postfix_seed_set.push_back({""}); }
    for (const auto& pre_s : prefix_seed_set) {
      if (pre_s.back() == base_seed.front()) { continue; }
      for (const auto& post_s : postfix_seed_set) {
        if (base_seed.back() == post_s.front()) { continue; }
        res.insert(pre_s + base_seed + post_s);
      }
    }
  };

  /**
   * Slide the `base_seed` from left to right in the process of generate combinations
   * of prefix and postfix-composed seeds. Put the derived seeds into set `res`.
   *
   * @param seed_size Size of the composed seed.
   * @param base_seed Base seed to be appended.
   * @param res Result set of derived seeds.
   */
  void derive_all_comb_from_base_seed (
      uint8_t seed_size,
      const std::string& base_seed,
      std::unordered_set<std::string>& res) {
    for (int prefix_len = 0, postfix_len = seed_size - base_seed.size(); 
        postfix_len >= 0; 
        prefix_len++, postfix_len--) {
      derive_for_fixed_base_seed_pos(prefix_len, postfix_len, base_seed, res);
    }
  };

  /**
   * Derive seed of length k with base_seed whose length should be smaller than k.
   * The base_seed should slide from left-most to right-most position of k-length seed.
   * The class member `derived_seeds` should be appended from empty. 
   *
   * @params seed_size Length of the derived seed of length k in the previous desciption.
   */
  void derive_blocked_seeds(uint8_t seed_size) {
    for (const auto& base_seed : base_seeds) {
      derive_all_comb_from_base_seed(seed_size, base_seed, this->derived_seeds);
    }
  }

  void derive_blocked_seeds_rule2(uint8_t seed_size) {
    /* auto num_repeat_base = 7; */
    /* auto the_repeat_base = "G"; */
    /* append pattern (e.g. G.G.G.G.G.G) to base_seeds. After this, calling derive_blocked_seeds would append seeds generated from base_seeds to derived(blocked)_seeds_set. */
    auto append_block_with_rule2 = [this] (uint8_t num_repeat_base, std::string the_repeat_base) {
      auto deriving_seeds = std::vector<std::string>{the_repeat_base};
      auto append_by_one_specific = [] (std::vector<std::string>& v, char base) {
        for (auto& i : v) {
          i.push_back(base);
        }
      };
      for (int i = 0; i < num_repeat_base - 1; i++) {
        deriving_seeds = append_by_one(deriving_seeds);
        append_by_one_specific(deriving_seeds, the_repeat_base[0]);
      }
      std::ranges::move(deriving_seeds, std::back_inserter(base_seeds));
    };
    log(derived_seeds.size(), "rule1 blocked seeds size");

    append_block_with_rule2(7, "A");
    append_block_with_rule2(7, "C");
    append_block_with_rule2(7, "G");
    append_block_with_rule2(7, "T");
    derive_blocked_seeds(seed_size);
    log(derived_seeds.size(), "rule2 blocked seeds size");
  }
};
