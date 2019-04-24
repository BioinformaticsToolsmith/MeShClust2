/* -*- C++ -*- */
/*
 * Clock.cpp
 *
 * Author: Benjamin T James
 */

#include "Clock.h"
#include <chrono>
#include <ctime>

static const auto _begin = std::chrono::system_clock::now();

void Clock::stamp(std::string desc)
{
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> diff = end - _begin;
	std::cout << "timestamp " << desc << " " << diff.count() << std::endl;
}
