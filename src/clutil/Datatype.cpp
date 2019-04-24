/* -*- C++ -*- */
/*
 * Datatype.cpp
 *
 * Author: Benjamin T James
 */

#include "Datatype.h"
std::string _dt_datatype = "";

std::string Datatype::get()
{
	return _dt_datatype;
}

void Datatype::set(std::string s)
{
	_dt_datatype = s;
}
