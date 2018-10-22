#ifndef RANDOM_H // -*- C++ -*-
#define RANDOM_H
#include <random>

class Random {
	std::mt19937 rng;
public:
	Random() : rng(std::random_device()()) {}

	template<class T>
	T randMod(T max) {
		std::uniform_int_distribution<T> distribution(0, max-1);
		return distribution(rng);
	}

	double random() {
		std::uniform_real_distribution<double> distribution(0.0, 1.0);
		return distribution(rng);
	}
};

#endif
