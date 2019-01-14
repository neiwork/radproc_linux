#include "matrixInit.h"

void matrixInit(Matrix& m, size_t height, size_t width, double initValue )
{
		Vector row(width,initValue);
		Matrix aux(height,row);
		m = aux;
}

void matrixInitCopy(Matrix& m, size_t height, size_t width, Matrix copy)
{
		Vector row(width,0.0);
		Matrix aux(height,row);
		m = aux;
		
		for (std::size_t i=0;i<height;i++) {
			for (std::size_t j=0;j<width;j++) {
				m[i][j] = copy[i][j];
			}
		}
}

void matrixInitSum3(Matrix& m, size_t height, size_t width, Matrix sum1, Matrix sum2, Matrix sum3)
{
		Vector row(width,0.0);
		Matrix aux(height,row);
		m = aux;
		
		for (std::size_t i=0;i<height;i++) {
			for (std::size_t j=0;j<width;j++) {
					m[i][j] = sum1[i][j]+sum2[i][j]+sum3[i][j];
			}
		}
}

void matrixInit4(Matrix& m1, Matrix& m2, Matrix& m3, Matrix& m4,
				size_t height, size_t width, double initValue)
{
		Vector row(width,initValue);
		Matrix aux(height,row);
		m1 = aux;  m2 = aux;  m3 = aux;  m4 = aux;
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