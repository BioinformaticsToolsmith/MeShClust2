/* -*- C++ -*-
 *
 * Runner.h
 *
 * Author: Benjamin T James
 */
#ifndef RUNNER_H
#define RUNNER_H

#include <iostream>
#include <map>
#include "Point.h"
#include "HandleSeq.h"
#include "Predictor.h"
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
	int sample_size = 0;
	int pivots = 40;
	int min_n_feat = 3;
	int max_n_feat = 5;
	int mut_type = HandleSeq::SINGLE;
	uint64_t feat_type = PRED_FEAT_FAST;
	double min_id = 0.35;
	std::vector<std::string> files;
	string output = "output.clstr";
	void get_opts(int argc, char** argv);
	pair<int,uint64_t> find_k();
};
#endif
