// -*- C++ -*-
/*
 * FeatureSelector.h
 *
 * Author: Benjamin T James
 */

#ifndef FEATURE_SELECTOR_H
#define FEATURE_SELECTOR_H

#include "GLM.h"
#include "Feature.h"

template<class T>
class FeatureSelector {
public:
	virtual ~FeatureSelector() {};
	static std::pair<matrix::Matrix,matrix::Matrix> generate_feat_mat(const vector<pra<T> > &data, Feature<T>& feat, double cutoff);
	static std::pair<double, matrix::GLM> class_train(const vector<pra<T> > &data, Feature<T>& feat, double cutoff);
	static std::pair<double, matrix::GLM> regression_train(const vector<pra<T> > &data, Feature<T>& feat);
	static double regression_test(const vector<pra<T> >& data, Feature<T>& feat, const matrix::GLM& glm);
	static tuple<double,double,double> class_test(const vector<pra<T> >& data, Feature<T>& feat, const matrix::GLM& glm, double cutoff);

	virtual pair<Feature<T>*,matrix::GLM> train_regression(Feature<T>* tfeat, const vector<pra<T> > &training,const vector<pra<T> > &testing) = 0;
	virtual pair<Feature<T>*,matrix::GLM> train_class(Feature<T>* tfeat, const vector<pra<T> > &training,const vector<pra<T> > &testing, double id) = 0;
};
#endif
