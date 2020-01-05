#pragma once
#include "State.h"
#include <fstream>

void oneZoneDist(Particle& p, State& st);
void multiZoneInjection(Particle& p, State& st);
void multiZoneDist(Particle& p, State& st, double dt);
void flareEmission(State& st, Particle& p, ofstream& file);
void flareEmission2(State& st,Particle& p,ofstream& file1,ofstream& file2,ofstream& file3,ofstream& file4);