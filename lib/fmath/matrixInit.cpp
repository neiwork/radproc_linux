#include "matrixInit.h"

void matrixInit(Matrix& m, int height, int width, double initValue ) {

		Vector row(width,initValue);
		Matrix aux(height,row);
		m = aux;

}

using namespace std;
#include <fstream>

void matrixRead(Matrix& m, int rows, int columns) {
    
    ifstream f("probMatrix2.txt");
    
    double aux;
    
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
            f >> m[i][j];
}