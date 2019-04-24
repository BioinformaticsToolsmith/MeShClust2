/* -*- C++ -*-
 *
 * bvec.cpp
 *
 * Author: Benjamin T James
 */
#include "bvec.h"
#include <algorithm>
template<class T>
bvec<T>::bvec(vector<uint64_t>& lengths, uint64_t bin_size)
{
	uint64_t num_points = lengths.size();
	std::sort(std::begin(lengths), std::end(lengths));
	for (uint64_t i = 0; i < lengths.size(); i += bin_size) {
		begin_bounds.push_back(lengths[i]);
		// uint64_t last_index = std::min((uint64_t)lengths.size() - 1,
		// 			       i + bin_size - 1);
		//std::cout << "[" << i << " " << last_index << "]" << std::endl;
	}
	data.reserve(begin_bounds.size());
	for (uint64_t i = 0; i < begin_bounds.size(); i++) {
		data.push_back({});
	}
}

template<class T>
Point<T>* bvec<T>::pop()
{
	for (auto& bin : data) {
		if (!bin.empty()) {
			Point<T>* p = bin[0].first;
			bin.erase(std::begin(bin));
			return p;
		}
	}
	return NULL;
}

template<class T>
Point<T>* bvec<T>::peek() const
{
	for (auto& bin : data) {
		if (!bin.empty()) {
			Point<T>* p = bin[0].first;
			return p;
		}
	}
	return NULL;
}

template<class T>
bool bvec<T>::inner_index_of(uint64_t length, size_t &idx, size_t *pfront, size_t *pback) const
{

	if (data.at(idx).empty() || idx == data.size()) {
		if (pfront) {
			for (size_t i = 0; i < data.size(); i++) {
				if (!data.at(i).empty()) {
					idx = i;
					*pfront = 0;
					break;
				}
			}
		}
		if (pback) {
			for (int i = data.size()-1; i >= 0; i--) {
				if (!data.at(i).empty()) {
					idx = i;
					*pback = 0;
					break;
				}
			}
		}
		return true;
	}
	size_t front = 0, back = 0;
	size_t low = 0, high = data.at(idx).size() - 1;
	bool found = false;
	if (length < data[idx][low].first->get_length() && pfront != NULL) {
		*pfront = low;
	}
	if (length > data[idx][high].first->get_length() && pback != NULL) {
		*pback = high;
	}
	for (;low <= high;) {
		size_t mid = (low + high) / 2;
		uint64_t d = data[idx][mid].first->get_length();
		if (d == length) {
			front = mid;
			back = mid;
			found = true;
			break;
		} else if (length < d) {
			high = mid;
		} else if (length > d) {
			low = mid + 1;
		}
		if (low == high) {
			found = true;
			front = low;
			back = high;
			break;
		}
	}
	if (pfront) {
		for (long i = front; i >= 0
			     && data[idx][i].first->get_length() == length; i--) {
			front = i;
		}
		*pfront = front;
	}
	if (pback) {
		for (long i = back; i < data[idx].size()
			     && data[idx][i].first->get_length() == length; i++) {
			back = i;
		}
		*pback = back;
	}
	return true;
}

template<class T>
bool bvec<T>::index_of(uint64_t point, size_t* pfront, size_t* pback) const
{
	size_t low = begin_bounds.size()-1, high = 0;
	size_t prev = 0;
	size_t prev_index = 0;

	for (size_t i = 1; i < begin_bounds.size(); i++) {
		prev_index = i - 1;
		prev = begin_bounds[i-1];
		if (point >= prev && point < begin_bounds[i]) {
			low = std::min(low, prev_index);
			high = std::max(high, prev_index);
		}
	}
	if (point >= begin_bounds[begin_bounds.size()-1]) {
		high = std::max(high, begin_bounds.size()-1);
	}
	if (pfront) {
		*pfront = low;
	}
	if (pback) {
		*pback = high;
	}
	return true;
}

template<class T>
void bvec<T>::insert(Point<T> *p)
{
	uint64_t len = p->get_length();
	size_t front = 0, back = 0;
	bool good = index_of(len, &front, &back);
	if (!good || front > back) {
		std::cerr << "error: list is not sorted" << std::endl;
	}
	std::vector<size_t> min_sizes;
	size_t minimum = std::numeric_limits<size_t>::max();
	for (size_t i = front; i <= back; i++) {
		size_t sz = data[i].size();
		if (sz < minimum) {
			minimum = sz;
			min_sizes.clear();
			min_sizes.push_back(i);
		} else if (sz == minimum) {
			min_sizes.push_back(i);
		}
	}
	if (min_sizes.empty()) {
		std::cerr << "error: no bins to insert into, item not inserted" << std::endl;
	}
	auto mid_min = min_sizes[min_sizes.size() / 2];
	data.at(mid_min).push_back(std::make_pair(p, false));
	if (begin_bounds.at(mid_min) > len) {
		cerr << "Begin Insertion of " << len << " should not be in bin " << begin_bounds.at(mid_min) << endl;
		throw std::exception();
	}
	if (mid_min < begin_bounds.size() - 1 && begin_bounds.at(mid_min+1) < len) {
		cerr << "End Insertion of " << len << " should not be in bin " << begin_bounds.at(mid_min+1) << endl;
		throw std::exception();
	}

}

template<class T>
size_t bvec<T>::size() const
{
	size_t num_bins = data.size();
	size_t total_size = 0;
	for (size_t i = 0; i < num_bins; i++) {
		total_size += data[i].size();
	}
	return total_size;
}

template<class T>
size_t bvec<T>::report() const
{
	cout << "BVec: ";
	size_t num_bins = data.size();
	cout << "num_bins=" << num_bins << endl;
	size_t total_size = 0;
	for (size_t i = 0; i < num_bins; i++) {
		uint64_t next_bound = std::numeric_limits<uint64_t>::max();
		if (i + 1 < num_bins) {
			next_bound = begin_bounds[i+1];
		}
		cout << "Bin " << i << ": [" << begin_bounds[i] << " " << next_bound << "] size=" << data[i].size() << endl;
		total_size += data[i].size();
	}
	cout << "total_size=" << total_size << endl;
	return total_size;
}
template<class T>
void bvec<T>::insert_finalize()
{
	auto sorter = [](const std::pair<Point<T>*,bool> a, const std::pair<Point<T>*,bool> b) {
		return a.first->get_length() < b.first->get_length();
	};
	for (size_t i = 0; i < data.size(); i++) {
		std::sort(std::begin(data[i]), std::end(data[i]), sorter);
		data[i].shrink_to_fit();
/*		if (data[i][0].first->get_length() < begin_bounds[i]) {
			cerr << "Length " << data[i][0].first->get_length() << " should not be in bin " << begin_bounds[i] << endl;
			throw std::exception();
		}
		if (i < data.size()-1 && data[i][data[i].size()-1].first->get_length() > begin_bounds[i+1]) {
			cerr << "Length " << data[i][0].first->get_length() << " should not be in bin " << begin_bounds[i] << " to " << begin_bounds[i+1] << endl;
			throw std::exception();
		}*/
	}
}

template<class T>
bool bvec<T>::empty() const
{
	bool is_empty = true;
	for (auto bin : data) {
		if (!bin.empty()) {
			is_empty = false;
			break;
		}
	}
	return is_empty;
}


template<class T>
uint64_t bvec<T>::absolute_idx(bvec_idx_t idx) const
{
	uint64_t ptr = 0;
	for (int i = 0; i < idx.first; i++) {
		ptr += data[i].size();
	}
	ptr += idx.second;
	return ptr;
}

template<class T>
std::pair<bvec_idx_t, bvec_idx_t>
bvec<T>::get_range(uint64_t begin_len, uint64_t end_len) const
{
	/* perform binary search to find bin */
	bvec_idx_t front, back;
	front.first = 0;
	front.second = 0;
	back.first = data.size()-1;
        back.second = data[back.first].size() - 1;

	// Determination of the outer indices
	if (!index_of(begin_len, &front.first, NULL)) {
		throw 100;
	}
	if (!index_of(end_len, NULL, &back.first)) {
		throw 100;
	}
	// if (begin_len < begin_bounds.at(front.first)) {
	// 	cerr << "Low index is not accurate" << endl;
	// 	throw std::exception();
	// }
	// if (front.first > back.first) {
	// 	cerr << "Front index is greater than back index" << endl;
	// 	throw std::exception();
	// }
	// Determination of the inner indices
	if (!inner_index_of(begin_len, front.first, &front.second, NULL)) {
		throw 100;
	}
	if (!inner_index_of(end_len, back.first, NULL, &back.second)) {
		throw 100;
	}
	// if (back.first != data.size()) { // ++ to make it an end iterator
	// 	if (back.second != data[back.first].size()) {
	// 		back.second++;
	// 	} else {
	// 		back.first++;
	// 		back.second = 0;
	// 	}
	// } else {
	// 	throw 101;
	// }

	if (back.first == (uint64_t)-1 || back.second == (uint64_t)-1) {
		back.is_empty = true;
	}
// 	for (uint64_t i = front.first; i <= back.first; i++) {
// 		uint64_t j = 0;
// 		uint64_t end = data.at(i).size();
// 		if (i == front.first) {
// 			j = front.second;
// 		}
// 		if (i == back.first) {
// 			end = min(back.second, end);
// 		}
// 		for (; j < end; j++) {
// 			uint64_t len = data.at(i).at(j).first->get_length();
// 			if (len < begin_len || len > end_len) {
// //				cerr << "Warning: Length in BVec " << len << " is not in [" << begin_len << ", " << end_len << "]. The classifier will not select these points." << endl;
// 				// if (i == front.first) {
// 				// 	cerr << "Front Bounds of selected bin: " << j << " -> " << data.at(i).at(j).first->get_length();
// 				// }
// 				// if (i == back.first) {
// 				// 	cerr << "End Bounds of selected bin: " << end-1 << " -> " << data.at(i).at(end-1).first->get_length();
// 				// }
// 			}
// 		}
// 	}
	return std::make_pair(front, back);
}

template<class T>
void bvec<T>::erase(size_t r, size_t c)
{
	data.at(r).erase(data.at(r).begin() + c);
}

/*
 * TODO: change available to Center class so no intermediate copying is done
 */
template<class T>
void bvec<T>::remove_available(bvec_idx_t begin, bvec_idx_t end, std::vector<Point<T>*> &available)
{
	size_t a = begin.first;
	size_t b = end.first;
	if (begin.is_empty || end.is_empty) {
		return;
	}
	int num = 0, new_num = 0;
	auto func = [](const bv_data_type<T> d) { return d.second; };
	auto inserter = [&](const std::pair<Point<T>*,bool> p) {
		if (p.second) {
#pragma omp critical
			available.push_back(p.first);
		}
	};
	// uint64_t begin_len = data.at(a).at(begin.second).first->get_length();
	// uint64_t end_len;
	// if (end.second == data.at(b).size()) {
	// 	end_len = data.at(b).at(end.second-1).first->get_length();
	// } else {
	// 	end_len = data.at(b).at(end.second).first->get_length();
	// }
	// cout << "Boundary: " << begin_len << " -> " << end_len;
	#pragma omp parallel for
	for (size_t i = a; i <= b; i++) {
		/* move marked points to end of vector, then copy, then erase */
		//const auto last = std::remove_if(std::begin(data[i]), std::end(data[i]), func);
		for (int j = 0; j < data[i].size(); j++) {
			auto kv = data[i][j];
			if (kv.second) {
#pragma omp critical
				{
					//cout << "Added " << kv.first->get_header() << " length: " << kv.first->get_length() << endl;
					// if (kv.first->get_length() > end_len || kv.first->get_length() < begin_len) {
					// 	cerr << "Error in bvec" << endl;
					// }
					available.push_back(kv.first);
				}
			}
		}
		data[i].erase(std::remove_if(std::begin(data[i]), std::end(data[i]), func), std::end(data[i]));
	}
}


template<class T>
bvec_iterator<T> bvec<T>::iter(bvec_idx_t idx)
{
	return bvec_iterator<T>(idx.first, idx.second, &data);
}


template class bvec<uint8_t>;
template class bvec<uint16_t>;
template class bvec<uint32_t>;
template class bvec<uint64_t>;
template class bvec<int>;
template class bvec<double>;
