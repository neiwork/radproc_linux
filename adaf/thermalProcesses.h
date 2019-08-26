#pragma once
#include "State.h"

void thermalProcesses(State&, const std::string&);
void writeBlob(State& st, Vector energies, Matrix lumOut,double tAccBlob);