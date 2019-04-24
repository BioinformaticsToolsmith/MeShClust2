// -*- C++ -*-
/*
 * BestFirstSelector.cpp
 *
 * Author: Benjamin T James
 */

#include "BestFirstSelector.h"
#include "../clutil/Progress.h"
#include <set>
#include <queue>

template<class T>
pair<Feature<T>*,matrix::GLM> BestFirstSelector<T>::train_regression(Feature<T>* feat, const vector<pra<T> > &training,const vector<pra<T> > &testing)
{
	matrix::GLM glm;
	return {NULL,glm};
}

using FeatPair = std::pair<uint64_t, Combo>;

struct Compare {
	bool operator()(const std::pair<set<FeatPair>,double> &a, const std::pair<set<FeatPair>,double> &b) {
		return a.second < b.second;
	}
};

using pqueue = std::priority_queue<std::pair<set<FeatPair>,double>, std::vector<std::pair<set<FeatPair>,double> >, Compare>;

vector<set<FeatPair> > children_of(set<FeatPair> feat,
				       const vector<FeatPair>& all_feats,
				       const set<set<FeatPair> >& closed_list,

				       const set<set<FeatPair> > &open_list)
{
	vector<set<FeatPair> > out;
	for (auto fp : all_feats) {
		set<FeatPair> temp = feat;
		auto pos = std::find(temp.begin(), temp.end(), fp);
		if (pos == temp.end()) {
			temp.insert(fp);
		} else {
			temp.erase(pos);
		}
		auto pos_bad = std::find(closed_list.begin(), closed_list.end(), temp);
		if (!temp.empty() && pos_bad == closed_list.end()) {
			auto pos_good = std::find(open_list.begin(), open_list.end(), temp);
			if (pos_good == open_list.end()) {
				out.push_back(temp);
			}
		}
	}
        return out;
}
template<class T>
std::string feat_name(Feature<T>* feat)
{
	std::ostringstream oss;
	auto feat_names = feat->feat_names();
	for (int i = 0; i < feat_names.size(); i++) {
		oss << feat_names[i];
		if (i < feat_names.size() - 1) {
			oss << " + ";
		}
	}
	return oss.str();
}
template<class T>
std::string feature_name(const set<FeatPair>& feat_list, int k)
{
	Feature<T> feat(k);
	for (auto fpair : feat_list) {
		feat.add_feature(fpair.first, fpair.second);
	}
	return feat_name(&feat);
}

template<class T>
Feature<T>* load_feat(Feature<T>* old_feat, const set<FeatPair>& feat_list, const vector<pra<T> > &training)
{

	Feature<T>* feat = new Feature<T>(*old_feat);
	feat->set_save(true);
	auto c_size = feat->get_combos().size();
	for (int i = 0; i < c_size; i++) {
		feat->remove_feature();
	}
	set<uint64_t> single_feats;
	for (FeatPair fp : feat_list) {
		for (uint64_t i = 1; i <= fp.first; i *= 2) {
			if (i & fp.first) {
				single_feats.insert(i);
			}
		}
		feat->add_feature(fp.first, fp.second);
	}
	for (uint64_t i : single_feats) {
		auto minmax = old_feat->get_normal(i);
		// #pragma omp critical
		// {
		// 	cout << "Feature " << Feature<T>::log2(i) << " min: " << minmax.first << " max: " << minmax.second << endl;
		// }
		feat->set_normal(i, minmax.first, minmax.second);
	}

//	feat->normalize(training);
	feat->finalize();

	return feat;
}

template<class T>
void calculate_table(Feature<T>* feat, vector<FeatPair> possible_feats, const vector<pra<T> > &training,const vector<pra<T> > &testing)
{
	auto c_size = feat->get_combos().size();
	for (int i = 0; i < c_size; i++) {
		feat->remove_feature();
	}
	for (FeatPair fp : possible_feats) {
		feat->add_feature(fp.first, fp.second);
	}
	feat->set_save(true);
	feat->normalize(training);
	feat->finalize();
	for (auto pr : testing) {
		feat->compute(*pr.first, *pr.second);
	}
}
template<class T>
pair<string,double> feature_accuracy(const set<FeatPair>& feat_set, Feature<T>* old_feat, const vector<pra<T> > &training,const vector<pra<T> > &testing, double id)
{
	Feature<T>* feat = NULL;
	feat = load_feat(old_feat, feat_set, training);
//	cout << "Considering " << name << " ";
	auto pr = FeatureSelector<T>::class_train(training, *feat, id);
	auto class_ac = FeatureSelector<T>::class_test(testing, *feat, pr.second, id);
	double class_accuracy = get<0>(class_ac);
//	cout << "Accuracy: " << class_accuracy << endl;
	std::string name = feat_name(feat);
	delete feat;
	return make_pair(name, class_accuracy);
}

pair<int,int> feat_list_sizes(const set<set<FeatPair> >& flist)
{
	int minimum = flist.begin()->size();
	int maximum = 0;
	for (auto item : flist) {
		if (item.size() < minimum) {
			minimum = item.size();
		}
		if (item.size() > maximum) {
			maximum = item.size();
		}
	}
	return make_pair(minimum, maximum);
}
template<class T>
void evaluate(const vector<set<FeatPair> > &item_list, set<set<FeatPair> > &open_map, pqueue &open_heap, Feature<T>* feat, const vector<pra<T> > &training,const vector<pra<T> > &testing, double id, int max_num_feat)
{
	std::ostringstream oss;
// of size " << item_list[0].size();


//	Progress prog(item_list.size(), oss.str());
	auto minmax = feat_list_sizes(open_map);
	oss << "Evaluating features " << minmax.second << "/" << max_num_feat;
	Progress prog(item_list.size(), oss.str());
	#pragma omp parallel for
	for (int i = 0; i < item_list.size(); i++) {
		const set<FeatPair>& item = item_list[i];
		auto feat_acc = feature_accuracy(item, feat, training, testing, id);
		double acc = feat_acc.second;
		string name = feat_acc.first;
		#pragma omp critical
		{
			prog++;
//			cout << name << ": " << acc << endl;
			open_map.insert(item);
			open_heap.push(std::make_pair(item, acc));
		}
	}
	prog.end();
}

template<class T>
std::pair<Feature<T>*,matrix::GLM> BestFirstSelector<T>::train_class(Feature<T>* feat, const vector<pra<T> > &training,const vector<pra<T> > &testing, double id)
{
	set<FeatPair> feat_set, best_feat_set;
	set<set<FeatPair> > closed_list, open_list;
	pqueue open_heap;

	int last_best_changed = 0;
	double best_acc = -100;
	const double eps = 0;

	cout << "Calculating all features" << endl;
	calculate_table(feat, possible_feats, training, testing);
	// prime the open_map
	vector<set<FeatPair> > children = children_of(feat_set, possible_feats, closed_list, open_list);
	evaluate(children, open_list, open_heap, feat, training, testing, id, max_num_feat);
        for (int iteration = 0; !open_list.empty(); iteration++) {

		auto minmax = feat_list_sizes(open_list);

		// stopping criteria: if we have already met the maximum number of features
		// or if no features changed in the last 3 iterations of having a minimum number of features
		if (minmax.second > max_num_feat || (iteration - last_best_changed >= 3 && minmax.second > min_num_feat)) {
			break;
		}
		//cout << "Features: " << minmax.first << " to " << minmax.second << endl;

		// Peek at the maximum-accuracy feature
		auto ptr = open_heap.top();
		feat_set = ptr.first;
		double acc = ptr.second;

		// Remove the best item from the open list/heap and add to the closed list
		open_heap.pop();
		open_list.erase(feat_set);
		closed_list.insert(feat_set);

		if (acc - eps > best_acc && feat_set.size() >= min_num_feat && feat_set.size() <= max_num_feat) {
			//cout << "New Best feature: " << feature_name<T>(feat_set, feat->get_k()) << endl;
			best_feat_set = feat_set;
			best_acc = acc;
			last_best_changed = iteration;
		}

		vector<set<FeatPair> > children = children_of(feat_set, possible_feats, closed_list, open_list);
		evaluate(children, open_list, open_heap, feat, training, testing, id, max_num_feat);
	}// while (iteration++ - last_best_changed < 2);

	Feature<T>* feat_c = load_feat(feat, best_feat_set, training);
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


template class BestFirstSelector<uint8_t>;
template class BestFirstSelector<uint16_t>;
template class BestFirstSelector<uint32_t>;
template class BestFirstSelector<uint64_t>;
template class BestFirstSelector<int>;
template class BestFirstSelector<double>;
