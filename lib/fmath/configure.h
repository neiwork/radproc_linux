#pragma once

#include "RungeKutta.h"
#include <boost/property_tree/ptree.hpp>


void fmath_configure(boost::property_tree::ptree& cfg)
{
	DefOpt_RungeKutta.samples_x = cfg.get<int>("math.runge-kutta-2.samples.x", DefOpt_RungeKutta.samples_x);
	DefOpt_RungeKutta.samples_y = cfg.get<int>("math.runge-kutta-2.samples.y", DefOpt_RungeKutta.samples_y);

	DefOpt_RungeKuttaSimple.samples_x = cfg.get<int>("math.runge-kutta-1.samples.x", DefOpt_RungeKuttaSimple.samples_x);
}
