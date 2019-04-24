/*
 * ChromListMaker.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Hani Zakaira Girgis
 */

#include "ChromListMaker.h"

namespace nonltr {

ChromListMaker::ChromListMaker(string seqFileIn, bool is_oneseq_) {
	seqFile = seqFileIn;
	is_oneseq = is_oneseq_;
	chromList = new vector<Chromosome *>();
}

ChromListMaker::~ChromListMaker() {
	Util::deleteInVector(chromList);
	delete chromList;
}


std::istream& safe_getline(std::istream& is, std::string& t)
{
	t.clear();
	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();
	for(;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if (sb->sgetc() == '\n') {
				sb->sbumpc();
			}
			return is;
		case std::streambuf::traits_type::eof():
			if (t.empty()) {
				is.setstate(std::ios::eofbit);
			}
			return is;
		default:
			t += (char)c;
		}
	}
}

const vector<Chromosome *> * ChromListMaker::makeChromList() {
	ifstream in(seqFile.c_str());
	bool isFirst = true;
	Chromosome * chrom;
	vector<uint64_t> size_list = getSize();
	uint64_t cur_seq = 0;
	if (is_oneseq) {
		uint64_t sum = 0;
		for (uint64_t len : size_list) {
			sum += len + 50;
		}
		size_list.clear();
		size_list.push_back(sum);
	}
	while (in.good()) {
		string line;
		safe_getline(in, line);
		if (line[0] == '>') {
			if (!isFirst) {
				if (is_oneseq) {
					std::string interseq(50, 'N');
					//	chrom->insert(interseq);
					chrom->appendToSequence(interseq);
				} else {
					chrom->finalize();
					chromList->push_back(chrom);
					chrom = new Chromosome(size_list.at(cur_seq++));
					chrom->setHeader(line);
				}
			} else {
				isFirst = false;
				chrom = new Chromosome(size_list.at(cur_seq++));
				chrom->setHeader(line);
			}
		} else if (line[0] == ' ' || line[0] == '\t') {
		} else {
			//	chrom->insert(line);
			chrom->appendToSequence(line);
		}
	}
	chrom->finalize();
	chromList->push_back(chrom);
	in.close();

	return chromList;
}

const vector<uint64_t> ChromListMaker::getSize() {
	ifstream in(seqFile.c_str());
	vector<uint64_t> size_list;
	uint64_t current_size = 0;
	while (in.good()) {
		string line;
		safe_getline(in, line);
		if (line[0] == '>') {
			if (current_size > 0) {
				size_list.push_back(current_size);
			}
			current_size = 0;
		} else if (line[0] == ' ' || line[0] == '\t') {
		} else {
			current_size += line.length();
		}
	}
	size_list.push_back(current_size);
	return size_list;
}
const vector<Chromosome *> * ChromListMaker::makeChromOneDigitDnaList() {
	ifstream in(seqFile.c_str());
	bool isFirst = true;
	ChromosomeOneDigitDna * chrom;
	vector<uint64_t> size_list = getSize();
	uint64_t cur_seq = 0;
	if (is_oneseq) {
		uint64_t sum = 0;
		for (uint64_t len : size_list) {
			sum += len + 50;
		}
		if (sum > 0) {
			sum -= 50;
		}
		size_list.clear();
		size_list.push_back(sum);
	}
	while (in.good()) {
		string line;
		safe_getline(in, line);
		if (line[0] == '>') {
			if (!isFirst) {
				if (is_oneseq) {
					std::string interseq(50, 'N');
					chrom->insert(interseq);
				} else {
					chrom->finalize();
					chromList->push_back(chrom);
					chrom = new ChromosomeOneDigitDna(size_list.at(cur_seq++));
					chrom->setHeader(line);
				}
			} else {
				isFirst = false;
				chrom = new ChromosomeOneDigitDna(size_list.at(cur_seq++));
				chrom->setHeader(line);

			}
		} else if (line[0] == ' ' || line[0] == '\t') {
		} else {
			chrom->insert(line);
//			chrom->appendToSequence(line);
		}
	}
	chrom->finalize();
	chromList->push_back(chrom);
	in.close();

	return chromList;
}

const vector<Chromosome *> * ChromListMaker::makeChromOneDigitProteinList() {
	ifstream in(seqFile.c_str());
	bool isFirst = true;
	ChromosomeOneDigitProtein * chrom;

	while (in.good()) {
		string line;
		safe_getline(in, line);
		if (line[0] == '>') {
			if (!isFirst) {
				chrom->finalize();
				chromList->push_back(chrom);
			} else {
				isFirst = false;
			}

			chrom = new ChromosomeOneDigitProtein();
			chrom->setHeader(line);
		} else {
			chrom->appendToSequence(line);
		}
	}

	chrom->finalize();
	chromList->push_back(chrom);
	in.close();

	return chromList;
}

}
/* namespace nonltr */
