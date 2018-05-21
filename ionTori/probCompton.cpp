#include "probCompton.h"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

void probCompton(double **prob, int size)
{
    ifstream file("prob.txt");
    
    for (int i=0;i<size;i++) {
        for (int j=0;j<size;j++) {
            file >> prob[i][j];
        }
    }
}