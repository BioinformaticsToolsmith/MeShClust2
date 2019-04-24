/* -*- C++ -*- */
/*
 * Trainer.cpp
 *
 * Author: Benjamin T James
 */
#include "Trainer.h"
#include "../predict/HandleSeq.h"
#include "../clutil/Datatype.h"
#include "../clutil/Loader.h"
#include "ClusterFactory.h"
#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include "../predict/Predictor.h"
#include "../predict/GLM.h"
#include "../predict/Feature.h"
#include "../clutil/Progress.h"
#include "../clutil/Random.h"

template<class T>
std::tuple<Point<T>*,double,size_t,size_t> Trainer<T>::get_close(Point<T> *p, bvec_iterator<T> istart, bvec_iterator<T> iend, bool &is_min_r) const
{
	int ncols = weights.getNumRow();
#pragma omp declare reduction(pmax:std::tuple<Point<T>*,double,size_t,size_t>: \
			      omp_out = get<1>(omp_in) > get<1>(omp_out) ? omp_in : omp_out ) \
	initializer (omp_priv=std::make_tuple((Point<T>*)NULL,-1,0,0))

	std::tuple<Point<T>*,
		   double,
		   size_t,
		   size_t> result = std::tuple<Point<T>*, double, size_t, size_t>(NULL,
				     -1,
				     0,
				     0);
	bool has_found = false;
	bool is_min = true;
	uint64_t min_len = p->get_length() * cutoff;
	uint64_t max_len = p->get_length() / cutoff;
#pragma omp parallel for reduction(pmax:result), reduction(&&:is_min)
	for (bvec_iterator<T> i = istart; i < iend; ++i) {
		Point<T>* pt = (*i).first;

		uint64_t len = pt->get_length();
		if (len < min_len || len > max_len) {
			continue;
		}
		auto cache = feat->compute(*pt, *p);
		double dist = (*feat)(0, cache);
		double sum = classify(pt, p);
		double res = round(sum) > 0;
		// #pragma omp critical
		// cout << "Result: " << sum << " raw_sigmoid: " << matrix::GLM::logistic(sum) << " classify_sum: " << Predictor<T>::classify_sum(sum) << " final: " << res << endl;
// set second to true if result is not 1.0
// which means it will be removed
		result = (dist > std::get<1>(result)) ? std::make_tuple(pt, dist, i.r, i.c) : result;
		is_min = is_min && (res != 1.0);
//		has_found = has_found || (res != 1.0);
		if (res == 1.0) {
			*i = std::make_pair(pt, true);
//			(*i).second = true;
		}
	}

//	is_min = !has_found;
	is_min_r = is_min;
//	return get<0>(result);
	return result;

}

template<class T>
long Trainer<T>::merge(vector<Center<T> > &centers, long current, long begin, long last) const
{
#pragma omp declare reduction(ldpmax:std::pair<long,double>:			\
			      omp_out = omp_in.second > omp_out.second ? omp_in : omp_out ) \
	initializer (omp_priv=std::make_pair(0, std::numeric_limits<double>::min()))
	std::pair<long,double> best = std::make_pair(0, std::numeric_limits<double>::min());
	Point<T>* p = centers[current].getCenter();
	uint64_t cen_length = p->get_length();
	uint64_t min_length = cen_length * get_id();
	uint64_t max_length = cen_length / get_id();
#pragma omp parallel for reduction(ldpmax:best)
	for (long i = begin; i <= last; i++) {
		double sum = weights.get(0, 0);
		double dist = 0;

		Point<T>* cen = centers[i].getCenter();
		uint64_t cen_len = cen->get_length();
		bool length_pass = cen_len >= min_length && cen_len <= max_length;
		if (length_pass) {
			auto cache = feat->compute(*cen, *p);
			for (int col = 1; col < weights.getNumRow(); col++) {
				double d = (*feat)(col-1, cache);
				if (col == 1) {
					dist = d;
				}
				sum += weights.get(col, 0) * d;
			}
			double res = round(Predictor<T>::classify_sum(sum));

			if (res == 1) {
				best = best.second > dist ? best : std::make_pair(i, dist);
			}
		}
	}
	return best.first;
}

template<class T>
double Trainer<T>::classify(Point<T>*a, Point<T>*b) const
{
	double sum = weights.get(0, 0);
	auto cache = feat->compute(*a, *b);
	for (int col = 1; col < weights.getNumRow(); col++) {
		sum += weights.get(col, 0) * (*feat)(col-1, cache);
	}
	return Predictor<T>::classify_sum(sum);
}

template<class T>
void Trainer<T>::filter(Point<T> *p, vector<pair<Point<T> *, bool> > &vec) const
{
	uint64_t cen_length = p->get_length();
	uint64_t min_length = cen_length * get_id();
	uint64_t max_length = cen_length / get_id();
	for (auto& pt : vec) {
		uint64_t pt_len = pt.first->get_length();
		bool length_pass = pt_len >= min_length && pt_len <= max_length;
		pt.second = true;
		if (length_pass) {
			double sum = classify(p, pt.first);
			double res = round(sum);
			pt.second = (res == 0);
		}
	}
	vec.erase(std::remove_if(vec.begin(), vec.end(), [](pair<Point<T>*, bool> p) {
				return p.second;
			}), vec.end());
}

template<class T>
Point<T>* Trainer<T>::closest(Point<double> *p, vector<pair<Point<T> *, bool> > &vec) const
{
	Point<T>* best_pt = NULL;
	double best_dist = 0;
	for (auto& pt : vec) {
		double sum = weights.get(0, 0);
		double dist = pt.first->distance_d(*p);
		if (best_pt == NULL || dist < best_dist) {
			best_dist = dist;
			best_pt = pt.first;
		}
	}
	return best_pt;
}

template<class T>
void Trainer<T>::train(std::string dump_str)
{
	Predictor<T> pred(dump_str);
	delete feat;
	auto pr = pred.get_class();
	feat = pr.first;
	feat->set_save(false);
	matrix::GLM glm = pr.second;
	weights = glm.get_weights();
}

template<class T>
void Trainer<T>::train(int min_n_feat, int max_n_feat, uint64_t feat_type, int mut_type, double min_id, std::string dump_str, double acc_cutoff)
{
	std::cout << "Splitting data" << endl;
	uintmax_t _id = points.size();
	Predictor<T> pred(k, cutoff, PRED_MODE_CLASS, feat_type,
			  mut_type, min_n_feat, max_n_feat, min_id);
	pred.train(points, _id, n_samples, n_templates);
	delete feat;
	auto pr = pred.get_class();
	feat = pr.first;
	matrix::GLM glm = pr.second;
	weights = glm.get_weights();

	if (dump_str != "") {
		pred.save(dump_str, Datatype::get());
		exit(0);
	} else {
		pred.save("weights.txt", Datatype::get());
	}
}

template class Trainer<uint8_t>;
template class Trainer<uint16_t>;
template class Trainer<uint32_t>;
template class Trainer<uint64_t>;
template class Trainer<int>;
template class Trainer<double>;
