/* -*- C++ -*-
 *
 * Loader.cpp
 *
 * Author: Benjamin T James
 *
 * Class which can 'preload' chunks of sequences from a file list,
 * and then count the k-mers separately, which can be done in
 * multiple threads
 */
#include "Loader.h"
#include "ClusterFactory.h"
#include "DivergencePoint.h"
#include <omp.h>

template<class T>
bool Loader<T>::done() const
{
	return file_idx == files.size();
}

template<class T>
void Loader<T>::preload(int tid)
{
	if (file_idx == files.size()) {
		return;
	}
	for (uint64_t j = 0; j < chunk_size; j++) {
		auto chrom = next();
		if (chrom.first == "") {
			return;
		}
		cache_list.at(tid).emplace_back(chrom.first, chrom.second);
	}
}


template<class T>
Point<T>* Loader<T>::get_point(std::string header, const std::string &base, uintmax_t& id, int k)
{
	KmerHashTable<unsigned long, T> table(k, 1);
	KmerHashTable<unsigned long, uint64_t> table_k1(1, 0);
	std::vector<T> values;
	vector<uint64_t> values_k1;
	values.clear();
	ChromosomeOneDigit chrom;
	chrom.setHeader(header);
	chrom.appendToSequence(base);
	chrom.finalize();
	fill_table<T>(table, &chrom, values);
	fill_table<uint64_t>(table_k1, &chrom, values_k1);
//	int tmplate = get_template(chrom->getHeader(), templates);
	Point<T> *p = new DivergencePoint<T>(values, chrom.size());
//	cout << "mag: " << ((DivergencePoint<T>*)p)->getPseudoMagnitude() << std::endl;
	p->set_1mers(values_k1);
	p->set_header(header);
	p->set_length(chrom.getBase()->length());
	p->set_data_str(*chrom.getBase());
	DivergencePoint<T>* q = dynamic_cast<DivergencePoint<T>*>(p);
	const auto N = q->points.size();
	double aq = (double) q->getPseudoMagnitude() / N;
	double sq = 0;
	for (auto i = 0; i < N; i++) {
		double qdiff = q->points[i] - aq;
		sq += qdiff * qdiff;
	}
	sq = sqrt(sq / N);
	q->set_stddev(sq);
	p->set_id(id);
	#pragma omp atomic
	id++;
	return p;
}

template<class T>
std::vector<Point<T>*> Loader<T>::load_next(int tid)
{
	std::vector<Point<T>*> points;
	for (size_t i = 0; i < cache_list.at(tid).size(); i++) {
	        auto pr = cache_list.at(tid).at(i);
		Point<T>* p = get_point(pr.first, *pr.second, id_list.at(tid), k);
		points.push_back(p);
		delete pr.second;
	}
	cache_list.at(tid).clear();
	return points;
}

template<class T>
std::pair<std::string,std::string*> Loader<T>::next()
{
	auto n = maker->next();
	if (n.first != "") {
		return n;
	}
	delete maker;
	maker = NULL;
	file_idx++;
	if (file_idx >= files.size()) {
		return n;
	}
	maker = new SingleFileLoader(files.at(file_idx));
	return maker->next();
}

template class Loader<double>;
template class Loader<int>;
template class Loader<uint64_t>;
template class Loader<uint32_t>;
template class Loader<uint16_t>;
template class Loader<uint8_t>;
