#pragma once

#include "State.h"

/* todas las perdidas estan en s⁻1*/
void timescales(Particle& p, State& st, const std::string& filename);
void radiativeLossesNeutron(Particle& n, State& st, const std::string& filename);



