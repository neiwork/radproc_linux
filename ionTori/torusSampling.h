#pragma once

#include <time.h>
#include <stdlib.h>
#include "State.h"
#include "torusParameters.h"
#include <fparticle/Particle.h>


//void torusSampling(Particle& p, double **prob);
void torusSampling(State& st, Matrix& prob);