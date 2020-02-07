#pragma once

/// Header Messages
#include <time.h>
#include <ostream>

enum MessageType {msgStart, msgEnd, msgError};

enum ProgramModule {Module_Main,
					Module_state,
					Module_targetField,
					Module_radLosses,
					Module_electronInjection,
					Module_electronDistribution,
					Module_protonInjection,
					Module_protonDistribution,
					Module_comptonScattMatrix,
					Module_photonInjection,
					Module_neutronInjection,
					Module_neutronDistribution,
					Module_thermalLuminosities,
					Module_thermalCompton,
					Module_luminosities,
					Module_secondariesLuminosities,
					Module_secondariesTimescales,
					Module_pionInjection,
					Module_pionDistribution,
					Module_muonInjection,
					Module_muonDistribution,
					Module_neutrinoInjection,
					Module_pairInjection,
					Module_pairDistribution,
					Module_Message}; 


extern void show_message(MessageType message, ProgramModule module);
//extern void show_error(ErrorType error, ProgramModule module);

extern void show_time(void);
extern double diff_seconds(clock_t timeStart, clock_t timeEnd);
extern void show_module(ProgramModule module);
extern void timestamp_stream(std::ostream& out);

