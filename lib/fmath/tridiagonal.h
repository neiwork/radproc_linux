#include <iostream>
#include <fmath/mathematics.h>

using namespace std;

//===========================================================================
void TriDiagSys(Vector a, Vector b, Vector c, Vector& d, int n)
//---------------------------------------------------------------------------
// Solves a system with tridiagonal matrix by LU factorization (diag(L) = 1).
// a - lower codiagonal (i=2,n)
// b - main diagonal (i=1,n)
// c - upper codiagonal (i=1,n-1)
// d - constant terms (i=1,n); solution on exit
// n - order of system.
//---------------------------------------------------------------------------
{
	int i;
	if (b[1] == 0e0) {
		cout << "TriDiagSys: singular matrix !" << endl; 
		return;
	}
	for (i=2; i<=n; i++) {				// factorization
		a[i] /= b[i-1];
		b[i] -= a[i]*c[i-1];
		if (b[i] == 0e0) {
			cout << "TriDiagSys: singular matrix !" << endl;
			return;
		}
		d[i] -= a[i]*d[i-1];
	}
	d[n] /= b[n];						// backward substitution
	for (i=n-1;i>=1;i--)
		d[i] = (d[i] - c[i]*d[i+1])/b[i];
}