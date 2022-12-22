#pragma once
#include <iostream>
#include <boost/program_options.hpp>
#include <filesystem>

#define VERBOSE false

class InputParams {
 public:
  std::filesystem::path reads_fasta = "/mnt/es/ness/johnson/thesis/data/Ecoli/reads/aligned_and_sampled/minimap2/DevNet-P6C4/rngAll/aln.lsmq10.fasta"; 
  std::filesystem::path nid_sampled_read_pairs;
  std::filesystem::path output_dir = "../output/devnet";
	std::filesystem::path fmi_dir = "../mid_files/devnet/";
	std::string fmi_src = "load";
  std::string anchors_offsets_source = "build";
  uint16_t kmer_size = 16;
  uint16_t dist_hori_filter = 800; 
	uint16_t bin_width = 400;
  uint16_t thres_one_bin_least_score = 75;
  uint16_t thres_bin_score = 75; 
  uint16_t numBin_sliding = 3;
	uint16_t thrsDcnt_toGraphBin = 200;
	uint32_t minNumSeedAnchors = 1;
	uint32_t maxNumSeedAnchors = 1000000;
	uint32_t batchSize_readsNum = 70000;
	int min_antiDiag_space = 30;
};

class ParseInput {
 public:
  int argc;
  char** argv;
  ParseInput(int argc_, char** argv_) 
  : argc(argc_), argv(argv_) {}

  void parse(InputParams& params) {
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("reads_fasta", po::value<std::filesystem::path>(), "nrd format template and query reads, both of which are the same for now.")
        ("fmi_dir", po::value<std::filesystem::path>(), "The directory for fm-index file to be saved or loaded.")
        ("fmi_src", po::value<std::string>(), "if this is experimental for accuracy test, one may want to skip the process of building fm-index built before, rather then building another one, while the opposite is for time performance estimation")
        ("kmer_size", po::value<uint16_t>(), "seed size")
        ("dist_hori_filter", po::value<uint16_t>(), "distance of horizontal anchor (anchors in one template search) for filtering.")
        ("bin_width", po::value<uint16_t>(), "bin_width")
        ("thres_one_bin_least_score", po::value<uint16_t>(), "threshold of count of continuous dots for deciding whether a seed should be considered.")
        ("maxNumSeedAnchors", po::value<uint32_t>(), "threshold maxNumSeedAnchors")
        ("nid_sampled_read_pairs", po::value<std::filesystem::path>(), "file containing numeric id of the sampled reads, should be endline separated.")
        ("output_dir", po::value<std::filesystem::path>(), "the directory of storing horizontally filtered hcrds.")
        ("thres_bin_score", po::value<uint16_t>(), "thres_bin_score")
        ("numBin_sliding", po::value<uint16_t>(), "numBin_sliding")
        ("thrsDcnt_toGraphBin", po::value<uint16_t>(), "thrsDcnt_toGraphBin")
        ("min_antiDiag_space", po::value<int>(), "min_antiDiag_space")
        ("batchSize_readsNum", po::value<uint32_t>(), "batchSize_readsNum")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(0);
    }

		auto set_vm = [&vm] <typename T>(std::string param_name, T& param, bool manditory) {
			if (vm.count(param_name)) {
				param = vm[param_name].as<T>();
				std::cout << param_name << " set to : " << param << "\n";
			} else {
				if (manditory) {
					std::cout << "loss manditory paramater " << param_name << ", exit(0)\n";
					exit(0);
				} else {
					std::cout << param_name << " (default): " << param << "\n";
				}
			}
		};

		set_vm("reads_fasta", params.reads_fasta, true);
		set_vm("fmi_dir", params.fmi_dir, true);
		set_vm("fmi_src", params.fmi_src, true);
		set_vm("kmer_size", params.kmer_size, true);
		set_vm("dist_hori_filter", params.dist_hori_filter, false);
		/* set_vm("bin_width", params.bin_width, false); */
		/* set_vm("thres_one_bin_least_score", params.thres_one_bin_least_score, false); */
		/* set_vm("maxNumSeedAnchors", params.maxNumSeedAnchors, true); */
		set_vm("thres_bin_score", params.thres_bin_score, true);
		/* set_vm("numBin_sliding", params.numBin_sliding, false); */
		set_vm("output_dir", params.output_dir, true);
		/* set_vm("thrsDcnt_toGraphBin", params.thrsDcnt_toGraphBin, true); */
		/* set_vm("min_antiDiag_space", params.min_antiDiag_space, true); */
		/* set_vm("batchSize_readsNum", params.batchSize_readsNum, true); */
  }
};

