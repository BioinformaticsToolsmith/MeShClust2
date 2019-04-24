/* -*- C++ -*- */
/*
 * GreedySelector.cpp
 *
 * Author: Benjamin T James
 */
#include "GreedySelector.h"
#include "../clutil/Progress.h"

template<class T>
pair<Feature<T>*,matrix::GLM> GreedySelector<T>::train_regression(Feature<T>* feat, const vector<pra<T> > &training,const vector<pra<T> > &testing)
{
	auto c_size = feat->get_combos().size();
	for (int i = 0; i < c_size; i++) {
		feat->remove_feature();
	}
	vector<uintmax_t> used_list;
	double abs_best_regr = 1000000;
//	Progress prog(possible_feats.size() * max_num_feat, "Feature selection:");
	for (auto num_feat = 1; num_feat <= max_num_feat; num_feat++) {
		double best_regr_err = abs_best_regr;
		uintmax_t best_idx = -1, cur_idx = 1;
		auto best_regr_feat = possible_feats.front();
		for (uint64_t i = 0; i < possible_feats.size(); i++) {
			if (std::find(used_list.begin(), used_list.end(), i) != used_list.end()) {
				continue;
			}
			auto rfeat = possible_feats[i];
		        feat->add_feature(rfeat.first, rfeat.second);
			feat->normalize(training);
			feat->finalize();
			auto pr = FeatureSelector<T>::regression_train(training, *feat);
			auto name = feat->feat_names().back();
			double regr_mse = FeatureSelector<T>::regression_test(testing, *feat, pr.second);
			feat->remove_feature();
			//	prog++;
			//cout << "Feature: " << cur_idx++ << "/" << possible_feats.size() - used_list.size() << " " << num_feat << "/" << max_num_feat << " " << name << " err: " << regr_mse << endl;
			if (regr_mse < best_regr_err) {
				best_regr_err = regr_mse;
				best_regr_feat = rfeat;
				best_idx = i;
			}
		}
		if (best_regr_err < abs_best_regr) {
			feat->add_feature(best_regr_feat.first, best_regr_feat.second);
			feat->normalize(training);
			feat->finalize();
			abs_best_regr = best_regr_err;
			used_list.push_back(best_idx);
			//possible_feats.erase(std::remove(possible_feats.begin(), possible_feats.end(), best_regr_feat), possible_feats.end());
		}
	}
//	prog.end();

	Feature<T>* feat_r = new Feature<T>(*feat);
	feat_r->set_save(false);
	auto pr = FeatureSelector<T>::regression_train(training, *feat_r);
	matrix::GLM r_glm = pr.second;
	double tr_regr_mse = FeatureSelector<T>::regression_test(testing, *feat_r, r_glm); // "training"
	cout << "Training Mean Error: " << pr.first << endl;
	double regr_mse = FeatureSelector<T>::regression_test(testing, *feat_r, r_glm);//, "testing");
	cout << "Testing Mean Error: " << regr_mse << endl;
	cout << "Features: "<< endl;
	for (auto line : feat_r->feat_names()) {
		cout << "\t" << line << endl;
	}
	auto w = r_glm.get_weights();
	for (int r = 0; r < w.getNumRow(); r++) {
		cout << "weight: ";
		for (int c = 0; c < w.getNumCol(); c++) {
			cout << w.get(r, c) << " ";
		}
		cout << endl;
	}

}

template<class T>
std::pair<Feature<T>*,matrix::GLM> GreedySelector<T>::train_class(Feature<T>* feat, const vector<pra<T> > &training,const vector<pra<T> > &testing, double id)
{
	auto c_size = feat->get_combos().size();
	for (int i = 0; i < c_size; i++) {
		feat->remove_feature();
	}
	vector<uintmax_t> used_list;
	double abs_best_acc = 0;
//	cout << "possible feats at one step: " << possible_feats.size() << endl;
	Progress prog(possible_feats.size() * max_num_feat, "Feature selection:");

	std::ostringstream oss;
	for (auto num_feat = 1; num_feat <= max_num_feat; num_feat++) {
		double best_class_acc = abs_best_acc;
		uintmax_t best_idx = -1, cur_idx = 1;
		auto best_class_feat = possible_feats.front();
		for (uint64_t i = 0; i < possible_feats.size(); i++) {
			if (std::find(used_list.begin(), used_list.end(), i) != used_list.end()) {
				continue;
			}
			auto rfeat = possible_feats[i];
		        feat->add_feature(rfeat.first, rfeat.second);
			feat->normalize(training);
			feat->finalize();
			auto name = feat->feat_names().back();
			auto pr = FeatureSelector<T>::class_train(training, *feat, id);
			auto class_ac = FeatureSelector<T>::class_test(testing, *feat, pr.second, id);
			double class_accuracy = get<0>(class_ac);//sqrt(get<1>(class_ac) * get<2>(class_ac));
			feat->remove_feature();
			prog++;
//			cout << "Feature: " << cur_idx++ << "/" << possible_feats.size() - used_list.size() << " " << num_feat << "/" << max_num_feat << " " << name  << " acc: " << get<0>(class_ac) << " sens: " << get<1>(class_ac) << " spec: " << get<2>(class_ac) << endl;
			if (class_accuracy > best_class_acc) {
				best_class_acc = class_accuracy;
				best_class_feat = rfeat;
				best_idx = i;
			}
		}
		/* accept the feature if either 1. we don't have enough features
		 * or 2. it improves accuracy by over 0.5%
		 */
		if (best_class_acc > abs_best_acc || num_feat <= min_num_feat) {
			feat->add_feature(best_class_feat.first, best_class_feat.second);
			feat->normalize(training);
			feat->finalize();
			abs_best_acc = best_class_acc;
			used_list.push_back(best_idx);
			oss << "Feature added: " << best_class_feat.first << " " << (int)best_class_feat.second << endl;
			oss << "Accuracy: " << best_class_acc << endl;
			possible_feats.erase(std::remove(possible_feats.begin(), possible_feats.end(), best_class_feat), possible_feats.end());
		}
	}
	prog.end();
	cout << oss.str();
	Feature<T>* feat_c = new Feature<T>(*feat);
	feat_c->set_save(false);
	auto pr = FeatureSelector<T>::class_train(training, *feat_c, id);
	matrix::GLM c_glm = pr.second;
	auto train_results = FeatureSelector<T>::class_test(training, *feat_c, c_glm, id);//, "train");
	cout << "Training ACC: " << get<0>(train_results) << " " << get<1>(train_results) << " " << get<2>(train_results) << endl;
	auto test_results = FeatureSelector<T>::class_test(testing, *feat_c, c_glm, id);//, "test");
	double class_acc = get<0>(test_results);
	cout << "Testing ACC: " << class_acc << " " << get<1>(test_results) << " " << get<2>(test_results) << endl;

	cout << "Features: "<< endl;
	for (auto line : feat_c->feat_names()) {
		cout << "\t" << line << endl;
	}
	return std::make_pair(feat_c, c_glm);
}

template class GreedySelector<uint8_t>;
template class GreedySelector<uint16_t>;
template class GreedySelector<uint32_t>;
template class GreedySelector<uint64_t>;
template class GreedySelector<int>;
template class GreedySelector<double>;
