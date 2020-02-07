#pragma once

#include "State.h"

void distributionFast(Particle& p, State& st);
void distributionDetailed(Particle& p, State& st);
void distributionOneZone(Particle& p, State& st);
void distributionMultiZone(Particle& p, State& st);
void distributionNew(Particle& p, State& st);
void distributionNewE(Particle& p, State& st);
void distributionGAMERA(Particle& p, State& st);
void distributionFokkerPlanck(Particle& p, State& st);
void distributionFokkerPlanckGamma(Particle& p, State& st);
void distributionFokkerPlanckMultiZone(Particle& p, State& st);
void distributionFokkerPlanckTimeDependent(Particle& p, State& st);
void distributionFokkerPlanckMultiZoneTimeDependent(Particle& p, State& st);
