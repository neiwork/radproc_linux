/// Messages
#define MAX_BUFFER_SIZE 255
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "messages.h"

void show_message(MessageType message, ProgramModule module)
{
    show_time();
    static clock_t startTime[255];

    switch(message)
    {
        case msgStart:
            startTime[module] = clock();
            printf("Starting ");
            show_module(module);
            break;

        case msgEnd:
            printf("Finished ");
            show_module(module);
            printf("(%gs) ", diff_seconds(startTime[module], clock()) );
            break;

        case msgError:
            printf("ERROR in: ");
            show_module(module);
            break;
        
//        default:
  //          show_error(Invalid_MessageType, Module_Message);
    //        break;
    }
    printf("\n");    
    return;
}


// Imprime la hora %H:%M:%S en pantalla.
void timestamp_stream(std::ostream& out)
{

	time_t rawTime;
	struct tm *timeInfo;
	char buffer[MAX_BUFFER_SIZE];

	//version windows
	//struct tm timeInfo;
	//time(&rawTime);
	//localtime_s(&timeInfo, &rawTime);
	//strftime(buffer, MAX_BUFFER_SIZE, "%H:%M:%S ", &timeInfo);

	time(&rawTime);     
	timeInfo = localtime(&rawTime);

	strftime(buffer, MAX_BUFFER_SIZE, "%H:%M:%S ", timeInfo);
	//printf("%s", buffer);

	out << buffer << std::endl;
}

// Imprime la hora %H:%M:%S en pantalla.
void show_time(void)
{
	
    time_t rawTime;
	struct tm *timeInfo;
    char buffer[MAX_BUFFER_SIZE];

    time(&rawTime);     
    timeInfo = localtime(&rawTime);
    
    strftime(buffer, MAX_BUFFER_SIZE, "%H:%M:%S ", timeInfo);
    printf("%s", buffer);
    
    return;
}

void show_module(ProgramModule module)
{

    switch(module)
    {
        case Module_Main:
            printf("main ");
            break;

		case Module_state:
            printf("Model structure and parameters");
            break;
        
		case Module_targetField:
            printf("Target fields");
            break;
			
		case Module_radLosses:
            printf("Non-thermal radiative losses");
            break;
	
		case Module_electronInjection:
            printf("Injection of electrons");
            break;
			
		case Module_comptonScattMatrix:
			printf("Compton Scattering Probability Matrix");
			break;
      
        case Module_Message:
            printf("message ");
            break;

		case Module_electronDistribution:
            printf("Electron Distribution ");
            break;
		
		case Module_protonDistribution:
            printf("Proton Distribution ");
            break;

		case Module_thermalLuminosities:
			printf("Thermal SEDs");
			break;
			
		case Module_thermalCompton:
			printf("Thermal Compton");
			break;
			
		case Module_neutronInjection:
			printf("Neutron Injection");
			break;
			
		case Module_luminosities:
			printf("Non-thermal luminosities");
			break;
  /*      default:
            printf("\n");
            show_error(Invalid_Module, Module_Message);
            break;*/ 
	}

//    return;
}

double diff_seconds(clock_t timeStart, clock_t timeEnd)
{
    return ((double )(timeEnd - timeStart))/CLOCKS_PER_SEC;
}

/* void show_error(ErrorType error, ProgramModule module)
{
    printf("\t ");
    switch(error)
    {
        case Invalid_Energy:
            show_message(msgError, module);
            printf("Invalid energy.");
            break;
        case Invalid_Particle:
            show_message(msgError, module);
            printf("Invalid particle type.");
            break;
        case Invalid_Pointer:
            show_message(msgError, module);
            printf("Invalid pointer.");
            break;
        case Invalid_MessageType:
            show_message(msgError, module);
            printf("Invalid message type.");
            break;
        case Invalid_Module:
            show_message(msgError, Module_Message);
            printf("Invalid Module.");
            break;
        case Invalid_Process:
            show_message(msgError, module);
            printf("Invalid process.");
            break;
        default:
            show_message(msgError, Module_Message);
            printf("Unknown error message.");
            break;
    }
    printf("\n");
    system("pause");
    exit(1);
}*/ 