/* -*- C++ -*-
 *
 * Progress.h
 *
 * Author: Benjamin T James
 */
#include <iostream>
#ifndef PROGRESS_H
#define PROGRESS_H

class Progress {
public:
	static void set_quiet(bool is_quiet_=true);
	static bool is_quiet;
	Progress(long num, std::string prefix_);
	~Progress() { end(); }
	void end();
	void operator++();
	void operator++(int);
	void operator+=(size_t);
private:
	void print();
	long pmax;
	long pcur;
	bool ended;
	std::string prefix;
	int barWidth;

	std::string last;
};
#endif
