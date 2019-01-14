#pragma once

#include <cstddef>
#include "mathematics.h"

void matrixInit(Matrix& m, size_t width, size_t height, double initValue);
void matrixInit4(Matrix& m1, Matrix& m2, Matrix& m3, Matrix& m4,
				size_t width, size_t height, double initValue);
void matrixInitCopy(Matrix& m, size_t width, size_t height, Matrix copy);
void matrixInitSum3(Matrix& m, size_t width, size_t height, Matrix sum1, Matrix sum2, Matrix sum3);
void matrixRead(Matrix& m, int rows, int columns);
