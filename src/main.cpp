#include <utils.hpp>
#include <options.hpp>
#include <seed_masker.hpp>
#include <file_io.hpp>
#include <reads_manip.hpp>
#include <config.hpp>
#include <detector.hpp>
#include <anchor_map.hpp>

/* #include <old_detect.hpp> */

int main(int argc, char** argv) { 
  auto timer = Timer();
  timer.start("Total");

  std::cout << " ------------ PARAMS ------------\n";
  auto params = InputParams(); 
  auto parser = ParseInput(argc, argv);
  parser.parse(params);

  std::cout << " ------------ CREATE_SEED_MASKER ------------\n";
  SeedMasker seed_masker;
  /* seed_masker.set_base_seed_rule_file(config::project_source_dir / "data/kmer_block_rule/rule1.txt"); */
  /* seed_masker.set_short_seeds(); */
  /* seed_masker.derive_blocked_seeds(params.kmer_size); */
  /* seed_masker.derive_blocked_seeds_rule2(params.kmer_size); */
  auto& blocked_seeds = seed_masker.get_derived_seeds();
  log(blocked_seeds.size(), "blocked_seeds.size()");

  std::cout << " ------------ INPUT_READS ------------\n";
  auto reads = FileIO::get_fasta_records(params.reads_fasta);

  std::cout << " ------------ PROCESS_READS ------------\n";
  ReadsManip::homopolymer_compression(reads);

  std::cout << " ------------ CREATE_SEED_TO_OFFSETS_MAP ------------\n";
  AnchorsMap anchors_map;
  anchors_map.build(params, reads, blocked_seeds);
  auto& anchors = anchors_map.get_anchors();
  
  std::cout << " ------------ DETECTION ------------\n";
  Detector detector;
  detector.detection(params, anchors_map, reads);

  timer.end();
}
