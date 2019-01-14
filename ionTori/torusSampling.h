#pragma once

#include <cstddef>
#include <time.h>
#include <stdlib.h>
#include "State.h"
#include "torusFunctions.h"
#include <fparticle/Particle.h>

void torusSampling(State& st, Matrix& prob, Vector& escape);