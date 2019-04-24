// -*- C++ -*-
/*
 * FeatureSelector.cpp
 *
 * Author: Benjamin T James
 */

#include "FeatureSelector.h"
template<class T>
std::pair<matrix::Matrix,matrix::Matrix> FeatureSelector<T>::generate_feat_mat(const vector<pra<T> > &data, Feature<T>& feat, double cutoff)
{
	bool classify = (cutoff > 0);
	int nrows = data.size();
	int ncols = feat.size()+1;
	matrix::Matrix feat_mat(nrows, ncols);
	matrix::Matrix labels(nrows, 1);
//	#pragma omp parallel for
	for (int row = 0; row < data.size(); row++) {
		auto kv = data.at(row);
		vector<double> cache;
 		// #pragma omp critical
		// {
			cache = feat.compute(*kv.first, *kv.second);
		// }
		feat_mat.set(row, 0, 1);
		if (classify) {
			labels.set(row, 0, kv.val >= cutoff ? 1 : -1);
		} else {
			labels.set(row, 0, kv.val);
			//	labels.set(row, 0, (kv.val - smin) / (smax - smin));
		}
		for (int col = 1; col < ncols; col++) {
			double val = feat(col-1, cache);
			feat_mat.set(row, col, val);
		}
	}
	return std::make_pair(feat_mat, labels);
}


template<class T>
std::pair<double, matrix::GLM> FeatureSelector<T>::regression_train(const vector<pra<T> > &data, Feature<T>& feat)
{
	auto pr = generate_feat_mat(data, feat, -1);
	matrix::GLM glm;
	glm.train(pr.first, pr.second);
	auto result1 = pr.first * glm.get_weights();
	auto diff1 = result1 - pr.second;
	double sum = 0;
	for (int i = 0; i < diff1.getNumRow(); i++) {
		sum += fabs(diff1.get(i, 0));
	}
	sum /= diff1.getNumRow();
	return {sum, glm};
}

template<class T>
std::pair<double, matrix::GLM> FeatureSelector<T>::class_train(const vector<pra<T> > &data, Feature<T>& feat, double cutoff)
{
	auto pr = generate_feat_mat(data, feat, cutoff);
	matrix::GLM glm;
	glm.train(pr.first, pr.second);
	matrix::Matrix p = glm.predict(pr.first);
	for (int row = 0; row < p.getNumRow(); row++) {
		if (p.get(row, 0) == 0) {
			p.set(row, 0, -1);
		}
	}
	auto tup = glm.accuracy(pr.second, p);
	double acc = get<0>(tup);
	double sens = get<1>(tup);
	double spec = get<2>(tup);
	return {acc, glm};
}

template<class T>
double FeatureSelector<T>::regression_test(const vector<pra<T> >& data, Feature<T>& feat, const matrix::GLM& glm)
{
	auto pr = generate_feat_mat(data, feat, -1);
	auto result1 = pr.first * glm.get_weights();
	auto diff1 = result1 - pr.second;
	double sum = 0;
	for (int i = 0; i < diff1.getNumRow(); i++) {
		sum += fabs(diff1.get(i, 0));
	}
	sum /= diff1.getNumRow();
	return sum;
}

template<class T>
tuple<double,double,double> FeatureSelector<T>::class_test(const vector<pra<T> >& data, Feature<T>& feat, const matrix::GLM& glm, double cutoff)
{
	auto pr = generate_feat_mat(data, feat, cutoff);
	matrix::Matrix p = glm.predict(pr.first);
	for (int row = 0; row < p.getNumRow(); row++) {
		if (p.get(row, 0) == 0) {
			p.set(row, 0, -1);
		}
	}
	auto tup = glm.accuracy(pr.second, p);
	return tup;

}

template class FeatureSelector<uint8_t>;
template class FeatureSelector<uint16_t>;
template class FeatureSelector<uint32_t>;
template class FeatureSelector<uint64_t>;
template class FeatureSelector<int>;
template class FeatureSelector<double>;
