#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>

#include <biovoltron/file_io/fasta.hpp>

class FileIO {
 public:
  static std::vector<biovoltron::FastaRecord<true>> get_fasta_records(const std::filesystem::path& file_reads) {
    auto res = std::vector<biovoltron::FastaRecord<true>>();
    auto fin = std::ifstream{file_reads};
    for (auto& r : std::ranges::istream_view<biovoltron::FastaRecord<true>>(fin)) {
      res.push_back(r);
    }
    return res;
  }
};
