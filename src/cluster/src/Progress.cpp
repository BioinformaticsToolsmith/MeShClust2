#include "Progress.h"
#include <iostream>
#include <sstream>

Progress::Progress(long num, std::string prefix_)
{
	pmax = num;
	ended = 0;
	pcur = 0;
	prefix = prefix_;
	last = "";
	barWidth = 70 - (prefix.size()+1);
	print();
}

void Progress::print()
{
	std::ostringstream oss;
	double prog = (double)pcur / pmax;
	oss << prefix << " [";
	int pos = barWidth * prog;
	for (int i = 0; i < barWidth; i++) {
		if (i < pos) {
			oss << "=";
		} else if (i == pos) {
			oss << ">";
		} else {
			oss << " ";
		}
	}
	oss << "] " << int(prog * 100.0) << " %\r";
	if (oss.str() != last) {
		last = oss.str();
		std::cout << last;
		std::cout.flush();
	}
}

void Progress::end()
{
	if (!ended) {
		pcur = pmax;
		print();
		std::cout << std::endl;
	}
	ended = true;
}

void Progress::operator++()
{
	pcur++;
	print();
}
void Progress::operator++(int)
{
	print();
	pcur++;
}


void Progress::operator+=(size_t num)
{
	pcur += num;
	print();
}
