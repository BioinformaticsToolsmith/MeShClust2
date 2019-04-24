/* -*- C++ -*- */
/*
 * Trainer.h
 *
 * Author: Benjamin T James
 */

#ifndef TRAINER_H
#define TRAINER_H

#include "../clutil/Point.h"
#include "../predict/GLM.h"
#include "../predict/Feature.h"
#include "../predict/Predictor.h"
#include "bvec.h"
#include "Center.h"
#include <set>

template<class T>
class Trainer {
public:
	Trainer(std::vector<Point<T>*> v, size_t num_points, size_t largest_count, double cutoff_, size_t max_pts_from_one_, int ksize) : points(v), n_samples(num_points), cutoff(cutoff_), n_templates(max_pts_from_one_), k(ksize) {
		uintmax_t size = 1000 * 1000 * 10;
		feat = new Feature<T>(k);
	};
	~Trainer() { delete feat; }
	void train(std::string);
	void train(int min_n_feat, int max_n_feat, uint64_t feat_type, int mut_type, double min_id, std::string dump_str, double acc_cutoff=97.5);

	std::tuple<Point<T>*,double,size_t,size_t> get_close(Point<T>*, bvec_iterator<T> istart, bvec_iterator<T> iend,  bool& is_min) const;

	void filter(Point<T>*, vector<pair<Point<T>*,bool> >&) const;
	Point<T>* closest(Point<double>*, vector<pair<Point<T>*,bool> >&) const;
	long merge(vector<Center<T> > &centers, long current, long begin, long end) const;

	double get_id() const { return cutoff > 1 ? cutoff / 100.0 : cutoff; }
private:
	double classify(Point<T>*, Point<T>*) const;
	matrix::Matrix weights;
	Feature<T> *feat;
	std::vector<Point<T>*> points;
	size_t n_samples, n_templates;
	double cutoff;
	int k;
};
#endif
