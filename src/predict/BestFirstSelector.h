// -*- C++ -*-
/*
 * BestFirstSelector.h
 *
 * Author: Benjamin T James
 */

#ifndef BEST_FIRST_SELECTOR_H
#define BEST_FIRST_SELECTOR_H
#include "FeatureSelector.h"
#include <set>
template<class T>
class BestFirstSelector : public FeatureSelector<T> {
public:
	BestFirstSelector(vector<std::pair<uint64_t, Combo> > possible_feats_, int min_n_feat, int max_n_feat) : possible_feats(possible_feats_), min_num_feat(min_n_feat), max_num_feat(max_n_feat) {}
	~BestFirstSelector() {}

	pair<Feature<T>*,matrix::GLM> train_regression(Feature<T>* tfeat, const vector<pra<T> > &training,const vector<pra<T> > &testing);
	pair<Feature<T>*,matrix::GLM> train_class(Feature<T>* tfeat, const vector<pra<T> > &training,const vector<pra<T> > &testing, double id);

private:
	int max_num_feat, min_num_feat;
	vector<std::pair<uint64_t, Combo> > possible_feats;
};
#endif
