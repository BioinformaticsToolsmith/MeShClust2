/* -*- C++ -*-
 *
 * Runner.h
 *
 * Author: Benjamin T James
 */
#ifndef CRUNNER_H
#define CRUNNER_H

#include <iostream>
#include <map>
#include "../clutil/Point.h"
#include "../predict/HandleSeq.h"
#include "../predict/Predictor.h"
using namespace std;

class Runner {
public:
	Runner(int argc, char** argv);
	~Runner() {};
	int run();
private:
	template<class T> int do_run();
	template<class T> void print_output(const map<Point<T>*, vector<Point<T>*>*> &m) const;
	int k = -1;
        int bandwidth;
	double similarity = 0.90;
	long largest_count = 0;
	int iterations = 15;
	int delta = 5;
	bool align = false;
	int total_sample_size = 2000;
	int n_templates = 300;
	int min_n_feat = 4;
	int max_n_feat = 4;
	bool is_single_file = false;
	double bias = 0;
	int mut_type = HandleSeq::BOTH;
	uint64_t feat_type = PRED_FEAT_FAST;
	double min_id = 0.35;
	std::vector<std::string> files, notrain_files, all_files;
	string output = "output.clstr";
	void get_opts(int argc, char** argv);
        int find_k();

	bool dump = false;
	bool recover = false;
	std::string dump_str = "weights.txt";
	Predictor<uint64_t> *pred64 = NULL;
};
#endif
