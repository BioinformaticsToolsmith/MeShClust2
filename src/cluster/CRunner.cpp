/* -*- C++ -*-
 *
 * Runner.cpp
 *
 * Author: Benjamin T James
 */
#include <vector>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <cstdlib>
#include "../nonltr/ChromListMaker.h"
#include "../clutil/Datatype.h"
#include "../clutil/Loader.h"
#include "../clutil/Clock.h"
#include "CRunner.h"
#include "Trainer.h"
#include "ClusterFactory.h"
#include "bvec.h"
#include "../clutil/Progress.h"
#include "../predict/Predictor.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * Constructor for runner.
 * Gets the options
 * If --recover is passed, set the K using that
 * Otherwise, if the K wasn't set, find the K by iterating through all sequences
 */
Runner::Runner(int argc, char **argv)
{
	get_opts(argc, argv);
	if (pred64 != NULL) {
		k = pred64->get_k();
	} else if (k == -1) {
		k = find_k();
	}
	srand(10);
}


/*
 * Main entry point for MeShClust2.
 * Sets datatype and identity if --recover was used
 * If datatype wasn't set, run through and find the max histogram
 * Based on the max histogram, call do_run with the smallest possible
 * histogram that will fit all sequences
 */
int Runner::run()
{
	if (pred64 != NULL) {
		Datatype::set(pred64->get_datatype());
		similarity = pred64->get_id();
	} else if (Datatype::get() == "") {
		largest_count = 0;
		Progress progress(all_files.size(), "Reading in sequences");
#pragma omp parallel for
		for (auto i = 0; i < all_files.size(); i++) {
			auto f = all_files.at(i);
			ChromListMaker maker(f, is_single_file);
			auto chromList = maker.makeChromOneDigitDnaList();


//		cout << "Reading in sequences from " << f << "..." << endl;
			uint64_t local_largest_count = 0;
//#pragma omp parallel for reduction(max:local_largest_count)
			for (int i = 0; i < chromList->size(); i++) {
				std::vector<uint64_t> values;
				KmerHashTable<unsigned long, uint64_t> table(k, 1);
				ChromosomeOneDigit *chrom = dynamic_cast<ChromosomeOneDigit *>(chromList->at(i));
				fill_table<uint64_t>(table, chrom, values);
				uint64_t l_count = *std::max_element(values.begin(), values.end());
				if (l_count > local_largest_count) {
					local_largest_count = l_count;
					// #pragma omp critical
					// cout << "local largest count reset to " << local_largest_count << endl;
				}
			}

			#pragma omp critical
			{
				if (local_largest_count > largest_count) {
					largest_count = local_largest_count;
				// #pragma omp critical
				// cout << "largest count updated to " << largest_count << endl;
				}
				progress++;
			}
		}
		progress.end();
		cout << "Largest count: " << largest_count << endl;
	}

	if (Datatype::get() != "") {
		std::string type = Datatype::get();
		if (type == "uint8_t") {
			largest_count = std::numeric_limits<uint8_t>::max();
		} else if (type == "uint16_t") {
			largest_count = std::numeric_limits<uint16_t>::max();
		} else if (type == "uint32_t") {
			largest_count = std::numeric_limits<uint32_t>::max();
		} else if (type == "uint64_t") {
			largest_count = std::numeric_limits<uint64_t>::max();
		}
	}
	if (largest_count <= std::numeric_limits<uint8_t>::max()) {
		cout << "Using 8 bit histograms" << endl;
		Datatype::set("uint8_t");
		return do_run<uint8_t>();
	} else if (largest_count <= std::numeric_limits<uint16_t>::max()) {
		cout << "Using 16 bit histograms" << endl;
		Datatype::set("uint16_t");
		return do_run<uint16_t>();
	} else if (largest_count <= std::numeric_limits<uint32_t>::max()){
	       	cout << "Using 32 bit histograms" << endl;
		Datatype::set("uint32_t");
		return do_run<uint32_t>();
	} else if (largest_count <= std::numeric_limits<uint64_t>::max()) {
	       	cout << "Using 64 bit histograms" << endl;
		Datatype::set("uint64_t");
		return do_run<uint64_t>();
	} else {
		throw "Too big sequence";
	}
}

void usage(std::string progname)
{
	std::cout << "Usage: " << progname << " --id 0.x [OPTIONS] *.fasta" << std::endl << std::endl;
	#ifndef VERSION
        #define VERSION "(undefined)"
        #endif
        std::cout << "Version " << VERSION << " compiled on " << __DATE__ << " " << __TIME__;
        #ifdef _OPENMP
        std::cout << " with OpenMP " << _OPENMP;
        #else
        std::cout << " without OpenMP";
        #endif
	std::cout << std::endl << std::endl;

	std::string raw = R"(--id          The most important parameter, --id, controls the identity cutoff of the sequences.
              Needs to be between 0 and 1.
              If it is not specified, an identity of 0.9 is used.

--kmer        decides the size of the kmers. It is by default automatically decided by average sequence
              length, but if provided, MeShClust can speed up a little by not having to find the largest
              sequence length. Increasing kmer size can increase accuracy, but increases memory consumption.

--dump       Run until the classifier is trained, and then dump the weights to the file,
             default 'weights.txt'. Can be used with --recover to recover the weights
             instead of re-training.

--recover    Recover weights for the classifier trained by a previous run which used --dump to dump
             the weights.

--list       Instead of specifying files as extra arguments, provide a text file with
             a list of files. Can use pipes or process substitutions such as "--list <(ls *.fasta) "

--no-train-list    Same as --list, but these files are not passed to the classifier,
                   e.g. unassembled genomes

--mut-type   {single, both, nonsingle-typical, nonsingle-all, all-but-reversion, all-but-translocation}
             changes the mutation generation algorithm. By default, "both" is used, utilizing
             single point and block mutations. On higher identity data sets, "single", which includes only single point mutations,
             is preferable. The option "nonsingle-typical" uses only block mutations,
             disallowing single point mutations. Other options include "all", which includes single,
             block, and nontypical mutations translocation and reversion.

--feat       determines the combinations of features to be used. By default, "slow" allows 11
             combinations to be selected from. "fast" removes 2 slower features from "slow"
             which include logarithm based features.

--single-file  Using this option, (no value is needed), each file is treated as a single sequence.
               If multiple sequences in a file are encountered, they are joined with 50 Ns,
               and the k-mers are not counted in that region.
               However, to be most accurate, it is advised to not use these sequences in the
               training step (for mutations) and instead 1) train using un-joined sequences and
               use --dump to dump to a file, and 2) use --recover with --single-file for the
               file list.

--sample     selects the total number of sequences used for both training and testing.
             2000 is the default value. That is, --sample 2000 provides 2000 training
             pairs and 2000 testing pairs.

--num-templates   selects the number of "template" sequences from which to mutate.
             For example, if 300 (the default) templates are requested, and the number of
             "samples" is requested to be 2000 (the default), 300 sequences will be read in
             and mutated 2000/300 times each to create 2000 semi-synthetic pairs.

--min-feat   (default 4) sets the minimum feature pairs to be used. If set to 2, at least 2 feature pairs
             will be used. Recall that features include pairwise combinations of the "feat" option.

--max-feat   (default 4) sets the maximum feature pairs to be used. Diminishing returns appears quickly,
             so a very large maximum (>10) is not advised.

--min-id     (default 0.35) sets the lower bound for mutation identity scores to be calculated.
             Shouldn't need to be set normally, as lower identites take much longer,
             especially with single mutations only.

--datatype   (8,16,32,64) Decides the integer size of the histograms. If not provided,
             all sequences are read in and counted to ensure the largest k-mer does not
             overflow. If the provided k-mer is too small, it will overflow.

--threads    sets the number of threads to be used. By default OpenMP uses the number of available cores
             on your machine, but this parameter overwrites that.

--output     specifies the output file, in CD-HIT's CLSTR format, described below:
             A '>Cluster ' followed by an increasing index designates a cluster.
             Otherwise, the sequence is printed out.
             A '*' at the end of a sequence designates the center of the cluster.
             An example of a small data set:

             >Cluster 0
             0       993nt, >seq128 template_6... *
             >Cluster 1
             0       1043nt, >seq235 template_10...
             1       1000nt, >seq216 template_10... *
             2       1015nt, >seq237 template_10...

--delta      decides how many clusters are looked around in the final clustering stage.
             Increasing it creates more accuracy, but takes more time. Default value is 5.

--iterations specifies how many iterations in the final stage of merging are done until convergence.
             Default value is 15.



If the argument is not listed here, it is interpreted as an input (FASTA format) file.


If you find this tool helpful, please cite:

James, Benjamin T. et al. (2018), MeShClust2: Application of alignment-free identity scores in clustering long DNA sequences. bioRxiv, 451278.

)";

	std::cout << raw << endl;
}


void Runner::get_opts(int argc, char **argv)
{
	for (int i = 1; i < argc; i++) {
		string arg = argv[i];
		if (arg == "--id" && i + 1 < argc) {
			try {
				std::string opt = argv[i+1];
				similarity = std::stod(opt);
				if (similarity <= 0 || similarity >= 1) {
					throw std::invalid_argument("");
				}
			} catch(std::exception e) {
				cerr << "Similarity must be between 0 and 1" << endl;
				exit(EXIT_FAILURE);
			}
			i++;
		} else if (arg == "--single-file") {
			is_single_file = true;
		} else if ((arg == "--list" || arg == "-l") && i + 1 < argc) {
			std::ifstream in(argv[++i]);
			std::string line;
			while (getline(in, line)) {
				files.push_back(line);
			}
		} else if ((arg == "--no-train-list" || arg == "--notrain-list") && i + 1 < argc) {
			std::ifstream in(argv[++i]);
			std::string line;
			while (getline(in, line)) {
				notrain_files.push_back(line);
			}
		} else if (arg == "--dump") {
			if (i + 1 < argc && argv[i+1][0] != '-') {
				dump_str = argv[++i];
			}
			dump = true;
		} else if ((arg == "--datatype") && i + 1 < argc) {
			std::string val = argv[++i];
			if (val == "uint8_t" || val == "8" || val == "uint8") {
				Datatype::set("uint8_t");
			} else if (val == "uint16_t" || val == "16" || val == "uint16") {
				Datatype::set("uint16_t");
			} else if (val == "uint32_t" || val == "32" || val == "uint32") {
				Datatype::set("uint32_t");
			} else if (val == "uint64_t" || val == "64" || val == "uint64") {
				Datatype::set("uint64_t");
			} else {
				cerr << "Histogram data type must have a valid data type or size: one of 8, 16, 32, 64" << endl;
				exit(EXIT_FAILURE);
			}
	        } else if ((arg == "-r" || arg == "--recover") && i + 1 < argc) {
			dump_str = argv[++i];
			recover = true;
			pred64 = new Predictor<uint64_t>(dump_str);
			similarity = pred64->get_id();
			k = pred64->get_k();
		} else if (arg == "--min-id" && i + 1 < argc) {
			try {
				std::string opt = argv[i+1];
				min_id = std::stod(opt);
				if (min_id <= 0 || min_id >= 1) {
					throw std::invalid_argument("");
				}
			} catch(std::exception e) {
				cerr << "Similarity must be between 0 and 1" << endl;
				exit(EXIT_FAILURE);
			}
			i++;
		} else if ((arg == "-b" || arg == "--bias") && i + 1 < argc) {
			bias = std::stod(argv[++i]);
		} else if ((arg == "-k" || arg == "--kmer") && i + 1 < argc) {
			k = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (k <= 0) {
				fprintf(stderr, "K must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
			output = string(argv[i+1]);
			i++;
		} else if ((arg == "--num-templates") && i + 1 < argc) {
			n_templates = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (n_templates <= 0) {
				fprintf(stderr, "Number of templates must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else if ((arg == "-s" || arg == "--sample") && i + 1 < argc) {
			total_sample_size = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (total_sample_size <= 0) {
				fprintf(stderr, "Sample size must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else if ((arg == "--mut-type") && i + 1 < argc) {
			std::string opt = argv[i+1];
			if (opt == "all") {
				mut_type = HandleSeq::BOTH | HandleSeq::ATYPICAL;
			} else if (opt == "both") {
				mut_type = HandleSeq::BOTH;
			} else if (opt == "snp" || opt == "single") {
				mut_type = HandleSeq::SINGLE;
			} else if (opt == "nonsingle-typical") {
				mut_type = HandleSeq::NON_SINGLE;
			} else if (opt == "nonsingle-all") {
				mut_type = HandleSeq::NON_SINGLE | HandleSeq::ATYPICAL;
			} else if (opt == "all-but-reversion") {
				mut_type = HandleSeq::BOTH | HandleSeq::TRANSLOCATION;
			} else if (opt == "all-but-translocation") {
				mut_type = HandleSeq::BOTH | HandleSeq::REVERSION;
			} else {
				cerr << "Options for mutation type are \"single\", \"nonsingle-typical\", \"both\" (for single and nonsingle-typical), \"nonsingle-all\", and \"all\" (single, nonsingle, and atypical nonsingle)." << endl;
				exit(1);
			}
			i++;
		} else if ((arg == "--feat" || arg == "-f") && i + 1 < argc) {
			std::string opt = argv[i+1];
			if (opt == "fast") {
				feat_type = PRED_FEAT_FAST;
			} else if (opt == "slow") {
				feat_type = PRED_FEAT_FAST | PRED_FEAT_DIV;
			} else if (opt == "extraslow") {
				feat_type = PRED_FEAT_ALL;
			} else {
				cerr << "Options for feature sets are \"fast\", \"slow\", and \"extraslow\"." << endl;
				exit(1);
			}
			i++;
		} else if ((arg == "--min" || arg == "--min-feat") && i + 1 < argc) {
			try {
				std::string opt = argv[i+1];
				int xx = std::stoi(opt);
				if (xx <= 0) {
					throw std::invalid_argument("");
				}
				min_n_feat = xx;
			} catch (std::exception e) {
				cerr << "Minimum number of features must be greater than 0." << endl;
				exit(1);
			}

			i++;
		} else if ((arg == "--max" || arg == "--max-feat") && i + 1 < argc) {
			try {
				std::string opt = argv[i+1];
				int xx = std::stoi(opt);
				if (xx <= 0) {
					throw std::invalid_argument("");
				}
				max_n_feat = xx;
			} catch (std::exception e) {
				cerr << "Maximum number of features must be greater than 0." << endl;
				exit(1);
			}

			i++;
		} else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
			try {
				std::string opt = argv[i+1];
				int threads = std::stoi(opt);
				if (threads <= 0) {
					throw std::invalid_argument("");
				}
				#ifdef _OPENMP
				omp_set_num_threads(threads);
				#endif
			} catch (std::exception e) {
				cerr << "Number of threads must be greater than 0." << endl;
				exit(1);
			}

			i++;

		} else if ((arg == "-d" || arg == "--delta") && i + 1 < argc) {
			delta = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (delta <= 0) {
				fprintf(stderr, "Delta must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else if ((arg == "-i" || arg == "--iter" || arg == "--iterations") && i + 1 < argc) {
			iterations = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (iterations <= 0) {
				fprintf(stderr, "Iterations must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else {
			struct stat st;
			stat(argv[i], &st);
			if (S_ISREG(st.st_mode)) {
				files.push_back(argv[i]);
			} else {
				usage(*argv);
				exit(EXIT_FAILURE);
			}
		}
	}
	set<std::string> file_list(files.begin(), files.end());
	set<std::string> notrain_file_list;
	set<std::string> all_files_list(files.begin(), files.end());
	for (std::string val : notrain_files) {
		if (file_list.find(val) == file_list.end()) {
			notrain_file_list.insert(val);
			all_files_list.insert(val);
		}
	}
	files.assign(file_list.begin(), file_list.end());
	notrain_files.assign(notrain_file_list.begin(),
			     notrain_file_list.end());
	all_files.assign(all_files_list.begin(),
			 all_files_list.end());
	if (files.empty()) {
		usage(*argv);
		exit(EXIT_FAILURE);
	}
	if (min_n_feat > max_n_feat) {
		cerr << "Minimum number of features (" << min_n_feat << ") cannot be greater than maximum number of features (" << max_n_feat << ")" << endl;
		exit(1);
	}
}

int Runner::find_k()
{
	unsigned long long count = 0, length = 0, largest_count = 0;
#pragma omp parallel for
	for (size_t i = 0; i < all_files.size(); i++) {
//		cout << "Processing " << f << endl;
		ChromListMaker maker(all_files.at(i), is_single_file);
		auto chromList = maker.makeChromList();
		unsigned long long l = 0;
	        for (int i = 0; i < chromList->size(); i++) {
			Chromosome *chrom = dynamic_cast<Chromosome *>(chromList->at(i));
			auto sz = chrom->getEffectiveSize();
			l += sz;
		}
		l /= chromList->size();
#pragma omp atomic
		length += l;
	}
	length /= files.size();
	int newk = ceil(log(length) / log(4)) - 1;
	cout << "avg length: " << length << endl;
	cout << "Recommended K: " << newk << endl;
	return newk;
}

template<class T>
void get_points(std::vector<std::string> files, std::vector<Point<T>*> &points, int k, uintmax_t &_id, bool is_single_file, bool set_seq=true, int chunk_size=10000)
{
	if (files.empty()) {
		return;
	}
	auto sort_func = [](Point<T>*a, Point<T>*b) {
		return a->get_length() < b->get_length();
	};
	auto sort_hdr = [](Point<T>*a, Point<T>*b) {
		return a->get_header() < b->get_header();
	};
	int n_threads = omp_get_max_threads();
	std::ostringstream oss;
	oss << "Counting " << k << "-mers";
	Progress prog(files.size(), oss.str());
	#pragma omp parallel for
	for (size_t i = 0; i < files.size(); i++) {
		ChromListMaker maker(files.at(i), is_single_file);
		auto chromList = maker.makeChromOneDigitDnaList();
		for (auto elt : *chromList) {
			ChromosomeOneDigitDna* chrom = dynamic_cast<ChromosomeOneDigitDna*>(elt);
			Point<T>* pt = Loader<T>::get_point(chrom, _id, k);
#pragma omp critical
			points.push_back(pt);
		}
#pragma omp critical
		prog++;
	}
	prog.end();
	std::string warning = Loader<T>::get_warning();
	if (warning != "") {
		cout << warning << endl;
	}
	std::sort(points.begin(), points.end(), sort_hdr);
	std::sort(points.begin(), points.end(), sort_func);
	// for (auto seq : points) {
	// 	cout << "SEQ " << seq->get_header() << endl;
	// }

}

/*
 * Main launching point for the algorithm
 *
 * Reads in all points
 * Trains the model using a previous model if provided
 * Else trains the model from scratch
 * Initializes bvec search structure
 * Runs mean shift
 */
template<class T>
int Runner::do_run()
{
	using pvec = vector<Point<T> *>;
	using pmap = map<Point<T>*, pvec*>;
	uintmax_t s_id = 0;
	Predictor<T>::set_bias(bias);

	pvec points;
	get_points(files, points, k, s_id, is_single_file);
	Clock::stamp("read_in_points");
	Trainer<T> tr(points, total_sample_size, largest_count, similarity, n_templates, k);
	if (recover) {
		tr.train(dump_str);
	} else {
		// If we are working in low-identity space, get more room
		if (similarity < 0.6) {
			min_id = 0.2;
		}
		tr.train(min_n_feat, max_n_feat, feat_type, mut_type, min_id, dump ? dump_str : "");
	}
	get_points(notrain_files, points, k, s_id, is_single_file, false);
	vector<uint64_t> lengths;
	for (Point<T>* p : points) {
		if (!align) {
			p->set_data_str("");
		}
		lengths.push_back(p->get_length());
	}
	// Initializing BVec
	bvec<T> bv(lengths, 1000);
	lengths.clear();
	// Inserting points into BVec
	uint64_t idx = 0;
	for (Point<T>* p : points) {
		p->set_id(idx++);
		bv.insert(p);
	}
	bv.insert_finalize();
	ClusterFactory<T> factory(k);
	factory.MS(bv, bandwidth, similarity, tr, output, iterations, delta);
	return 0;
}


template<class T>
void Runner::print_output(const map<Point<T>*, vector<Point<T>*>*> &partition) const
{
	cout << "Printing output" << endl;
	std::ofstream ofs;
	ofs.open(output, std::ofstream::out);
	int counter = 0;
	for (auto const& kv : partition) {
		if (kv.second->size() == 0) {
			continue;
		}
		ofs << ">Cluster " << counter << endl;
		int pt = 0;
		for (auto p : *kv.second) {
			string s = p->get_header();
			ofs << pt << "\t"  << p->get_length() << "nt, " << s << "... " << endl;
//			string fa = am.get(p->get_id());
//			ofs << writefa(fa) << endl;
			pt++;
		}
		counter++;
	}
	ofs.close();
}
