#pragma once
#include "State.h"

void thermalRadiation(State&, const std::string&);
void thermalTargetField(Particle&, const std::string&);
void writeBlob(State& st, Vector energies, Matrix lumOut,double tAccBlob);