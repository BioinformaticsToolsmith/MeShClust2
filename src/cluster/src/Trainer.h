/* -*- C++ -*- */
#ifndef TRAINER_H
#define TRAINER_H

#include "Point.h"
#include "GLM.h"
#include "Feature.h"
#include "bvec.h"
#include "Center.h"
#include "LogTable.h"
#include <set>

template<class T>
class Trainer {
public:
	Trainer(std::vector<Point<T>*> v, size_t num_points, size_t largest_count, double cutoff_, size_t max_pts_from_one_, double (&matrix)[4][4], double sig, double eps, int ksize) : points(v), n_points(num_points), cutoff(cutoff_), max_pts_from_one(max_pts_from_one_), k(ksize) {
		init(matrix, sig, eps);
		uintmax_t size = 1000 * 1000 * 10;
		feat = new Feature<T>(k);
	};
	~Trainer() { delete feat_mat; delete feat; }
	std::pair<std::map<std::pair<Point<T>*, Point<T>*>, double>,
		  std::map<std::pair<Point<T>*, Point<T>*>, double> > split_old();
        vector<std::pair<Point<T>*,Point<T>*> > split();
	double train_n(pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
		       Point<T>*> > > &data, int ncols);
	void train(int min_n_feat, int max_n_feat, uint64_t feat_type, int mut_type, double min_id, double acc_cutoff=97.5);
	void mutate_seqs(Point<T>* p, size_t num_seq, vector<pra<T> > &pos_buf, vector<pra<T> > &neg_buf, double id_begin, double id_end, uintmax_t& _id);
	std::tuple<Point<T>*,double,size_t,size_t> get_close(Point<T>*, bvec_iterator<T> istart, bvec_iterator<T> iend,  bool& is_min) const;
//	vector<pair<int, double> > get_close(Point<T>*, const vector<pair<Point<T>*,int> > &,  bool& is_min) const;
	std::pair<std::pair<vector<pair<Point<T>*,Point<T>*> >,
		    vector<pair<Point<T>*,Point<T>*> > >,
	  std::pair<vector<pair<Point<T>*,Point<T>*> >,
		    vector<pair<Point<T>*,Point<T>*> > > >
	new_get_labels(std::vector<Point<T>*> &points, size_t num_sample, double id, uintmax_t &_id);
	void filter(Point<T>*, vector<pair<Point<T>*,bool> >&) const;
	Point<T>* closest(Point<double>*, vector<pair<Point<T>*,bool> >&) const;
	long merge(vector<Center<T> > &centers, long current, long begin, long end) const;
//	Point<T>* merge(Point<T>*, vector<pair<Point<T>*,double> >&) const;
private:
	matrix::GLM glm;
	matrix::Matrix weights;
	double align(Point<T>* a, Point<T>* b) const;
	std::pair<matrix::Matrix,matrix::Matrix> generate_feat_mat(pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
					 Point<T>*> > > &data, int ncols);
	void init(double (&matrix)[4][4], double sig, double eps);
	pair<vector<pair<Point<T>*,
			 Point<T>*
			 > >,
	     vector<pair<Point<T>*,
			 Point<T>*> > > get_labels(vector<std::pair<Point<T>*,Point<T>*> >&, double cutoff) const;
	Feature<T> *feat;
	int mat[4][4];
	int sigma, epsilon;
	std::vector<Point<T>*> points;
	matrix::Matrix *feat_mat = NULL;
	size_t n_points, max_pts_from_one;
	double cutoff;
	int k;
};
#endif
