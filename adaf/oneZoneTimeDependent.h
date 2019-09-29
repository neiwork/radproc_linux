#pragma once
#include "State.h"
#include <fstream>

void oneZoneDist(Particle& p, State& st);
void flareEmission(State& st, Particle& p, ofstream& file);