#include "tridiagonal.h"
#include "matrixInit.h"

using namespace std;

//===========================================================================
void TriDiagSys2(Vector a, Vector b, Vector c, Vector& d, int n)
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

//===========================================================================
void TriDiagSys(Vector a, Vector b, Vector c, Vector& d, int n)
//---------------------------------------------------------------------------
// Solves a system with tridiagonal matrix by LU factorization (diag(L) = 1).
// a - lower codiagonal (i=1,n-1)
// b - main diagonal (i=0,n-1)
// c - upper codiagonal (i=0,n-2)
// d - constant terms (i=0,n-1); solution on exit
// n - order of system.
//---------------------------------------------------------------------------
{
	int i;
	if (b[0] == 0e0) {
		cout << "TriDiagSys: singular matrix !" << endl; 
		return;
	}
	for (i=1;i<n;i++) {				// factorization
		a[i] /= b[i-1];
		b[i] -= a[i]*c[i-1];
		if (b[i] == 0e0) {
			cout << "TriDiagSys: singular matrix !" << endl;
			return;
		}
		d[i] -= a[i]*d[i-1];
	}
	d[n-1] /= b[n-1];						// backward substitution
	for (i=n-2;i>=0;i--)
		d[i] = (d[i] - c[i]*d[i+1])/b[i];
}

//===========================================================================
void LUFactor(Matrix& a, vector<int>& ipivot, int n, double &det)
//---------------------------------------------------------------------------
// Performs LU factorization of (n x n) matrix a (diag(L) = 1). On exit,
// replaces upper triangle and diagonal with U, and lower triangle, with L.
// Uses partial pivoting on columns.
// a - coefficient matrix (n x n); LU decomposition on exit
// ipivot - array of pivot row indexes (output)
// det - determinant of coefficient matrix (output).
//---------------------------------------------------------------------------
{
	double amax, sum, t;
	int i, imax, j, k;
	det = 1e0;
	for (j=0;j<n;j++) {					// loop over columns
		for (i=0;i<=(j-1);i++) {			// elements of matrix U
			sum = a[i][j];
			for (k=0;k<=(i-1);k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		amax = 0e0;
		for (i=j;i<n;i++) {					// elements of matrix L
			sum = a[i][j];					// undivided by pivot
			for (k=0;k<=(j-1);k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;					// determine pivot
	
			if (amax < fabs(a[i][j])) {amax = fabs(a[i][j]); imax = i;}
		}
		if (amax == 0e0)
			{ printf("LUFactor: singular matrix !\n"); det = 0e0; return; }
		ipivot[j] = imax;					// store pivot row index
											// interchange rows imax and j
											// to put pivot on diagonal
		if (imax != j) {
			det = -det;
			for (k=0; k<n; k++)
				{ t = a[imax][k]; a[imax][k] = a[j][k]; a[j][k] = t; }
		}
		det *= a[j][j];							// multiply determinant with pivot
		t = 1e0/a[j][j];
		for (i=j+1; i<n; i++) a[i][j] *= t;	// divide elements of L by pivot
	}
}

//===========================================================================
void LUSystem(Matrix& a, vector<int>& ipivot, Vector& b, int n)
//---------------------------------------------------------------------------
// Solves linear system a x = b of order n by LU factorization.
// a - LU decomposition of coefficient matrix (returned by LUFactor)
// ipivot - array of pivot row indexes (input)
// b - vector of constant terms (input); solution x (on exit)
//---------------------------------------------------------------------------
{
	double sum;
	int i,j;
	for (i=0;i<n;i++) {
		sum = b[ipivot[i]];
		b[ipivot[i]] = b[i];
		for (j=0;j<=(i-1);j++) sum -= a[i][j]*b[j];
		b[i] = sum;
	} // solves Ly = b
	
	for (i=n-1;i>=0;i--) {
		sum = b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i] = sum/a[i][i];
	}
}

//===========================================================================
void MatInv(Matrix& a, int n, double &det)
//---------------------------------------------------------------------------
// Calculates inverse of real matrix by LU factorization.
// a - (n x n) matrix (input); a^(-1) (output)
// det - determinant of coefficient matrix (output).
// Calls: LUFactor, LUSystem.
//---------------------------------------------------------------------------
{
	Matrix ainv;
	Vector b(n,0.0);
	int i, j;
	vector<int> ipivot(n,0.0);
	matrixInit(ainv,n,n,0.0);
	LUFactor(a,ipivot,n,det);			// LU factorization of a
	
	if (det == 0e0)						// singular matrix
		{ printf("MatInv: singular matrix !\n"); return; }
	
	for (j=0;j<n;j++) {
		for (i=0;i<n;i++) b[i] = 0e0;
		b[j] = 1e0;
		LUSystem(a,ipivot,b,n);			// loop over columns of unit matrix
										// column j
										// solve systemSystems of Linear Equations
		for (i=0;i<n;i++) ainv[i][j] = b[i];			// column j of inverse
	}
	for (j=0;j<n;j++)
		for (i=0;i<n;i++) a[i][j] = ainv[i][j];		// copy inverse in a
}


void MatMul(Matrix a, Matrix& b, size_t m)
{
	Matrix c;
	matrixInitCopy(c,m,m,b);
	matrixInit(b,m,m,0.0);
	for (size_t k1=0;k1<m;k1++) {
		for (size_t k2=0;k2<m;k2++) {
			for (size_t j=0;j<m;j++)
					b[k1][k2] += a[k1][j] * c[j][k2];
		}
	}
}

//===========================================================================
void TriBlockDiagSys(Vector a, Vector ba, Vector bb, Vector bc, Vector c, Vector& d, size_t J, size_t M)
//---------------------------------------------------------------------------
// Solves a system with block tridiagonal matrix by LU factorization (diag(L) = 1).
// a - lower codiagonal (i=1,n-1)
// b - main diagonal (i=0,n-1)
// c - upper codiagonal (i=0,n-2)
// d - constant terms (i=0,n-1); solution on exit
// n - order of system.
//---------------------------------------------------------------------------
{
	Matrix alpha, beta;
	Vector x(J*M,0.0);
	matrixInit(beta,M,M,0.0);
	for (size_t row=0;row<M;row++) {
		for (size_t col=0;col<M;col++) {
			if (row == col) beta[row][col] = bb[row];
			else if (row == col+1) beta[row][col] = ba[row];
			else if (row == col-1) beta[row][col] = bc[row];
		}
	}
	
	double det = 0.0;
	for (size_t j=1;j<J;j++) {
		
		matrixInit(alpha,M,M,0.0);
		
		MatInv(beta,M,det);
		for (size_t col=0;col<M;col++)
			for (size_t row=0;row<M;row++)
				alpha[row][col] = beta[row][col] * a[j*M+col];
		for (size_t col=0;col<M;col++)
			for (size_t row=0;row<M;row++) {
				if (row == col) beta[row][col] = bb[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
				else if (row == col+1) beta[row][col] = ba[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
				else if (row == col-1) beta[row][col] = bc[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
				else beta[row][col] = - alpha[row][col] * c[(j-1)*M+col];
			}
		for (size_t k=0;k<M;k++)
			for (size_t kk=0;kk<M;kk++)
				d[j*M+k] -= alpha[k][kk]*d[(j-1)*M+kk];
	}
	
	MatInv(beta,M,det);
	for (size_t k=0;k<M;k++)
		for (size_t kk=0;kk<M;kk++)
			x[(J-1)*M+k] += beta[k][kk] * d[(J-1)*M+kk];
	
	for (size_t U=J-1;U>1;U--) {
		matrixInit(beta,M,M,0.0);
		for (size_t row=0;row<M;row++) {
			for (size_t col=0;col<M;col++) {
				if (row == col) beta[row][col] = bb[row];
				else if (row == col+1) beta[row][col] = ba[row];
				else if (row == col-1) beta[row][col] = bc[row];
			}
		}
		for (size_t j=1;j<U;j++) {
			
			matrixInit(alpha,M,M,0.0);
			
			MatInv(beta,M,det);
			for (size_t col=0;col<M;col++)
				for (size_t row=0;row<M;row++)
					alpha[row][col] = beta[row][col] * a[j*M+col];
			for (size_t col=0;col<M;col++)
				for (size_t row=0;row<M;row++) {
					if (row == col) beta[row][col] = bb[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
					else if (row == col+1) beta[row][col] = ba[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
					else if (row == col-1) beta[row][col] = bc[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
					else beta[row][col] = - alpha[row][col] * c[(j-1)*M+col];
				}
			for (size_t k=0;k<M;k++)
				for (size_t kk=0;kk<M;kk++)
					d[j*M+k] -= alpha[k][kk]*d[(j-1)*M+kk];
		}
		MatInv(beta,M,det);
		for (size_t k=0;k<M;k++)
			for (size_t kk=0;kk<M;kk++)
				x[(U-1)*M+k] += beta[k][kk] * ( d[(U-1)*M+kk] - c[(U-1)*M+kk] * x[U*M+kk] );
	}
	for (size_t i=0;i<J*M;i++) d[i] = x[i];
}

void GaussSeidel(Matrix& a, Vector& b, Vector& x, int n, int init, double &err)
//---------------------------------------------------------------------------
// Solves linear system a x = b by the Gauss-Seidel method.
// To ensure convergence, the system is left-multiplied with a^T.
// a - coefficient matrix (n x n)
// b - vector of constant terms
// x - initial approximation of solution (input); solution (output)
// n - order of system
// err - maximum relative error of the solution components
// init - initialization option: 0 - refines initial approximation
//                               1 - initializes solution
//---------------------------------------------------------------------------
{
	const double eps = 1e-2;			// relative precision criterion
	const int itmax = 1000;				// max no. of iterations
	double del, f;
	Vector t(n,0.0);
	Matrix s;
	matrixInit(s,n,n,0.0);
	int i, j, k;
	
	for (i=0;i<n;i++) {					// matrices of normal system
		for (j=0;j<=i;j++) {			// by multiplication with aT
			s[i][j] = 0e0;				// store result in s and t
			for (k=0;k<n;k++) s[i][j] += a[k][i]*a[k][j];
			s[j][i] = s[i][j];
		}
		t[i] = 0e0;
		for (j=0;j<n;j++) t[i] += a[j][i]*b[j];
	}
	
	for (i=0;i<n;i++) {				// matrices s and t of reduced system
		f = -1e0/s[i][i]; t[i] /= s[i][i];
		for (j=0;j<n;j++) s[i][j] *= f;
	}
	
	if (init) for (i=0;i<n;i++) x[i] = t[i];			// initialize solution
	
	for (k=1; k<=itmax; k++) {							// loop of iterations
		err = 0e0;
		for (i=0;i<n;i++) {
			del = t[i];				// correction
			for (j=0;j<n;j++) del += s[i][j]*x[j];
			x[i] += del;								// new approximation to solution
			if (x[i]) del /= x[i];						// relative error
			if (fabs(del) > err) err = log10(fabs(del));		// maximum error
		}
		if (err <= eps) break;							// check convergence
	}
	
	if (k > itmax) printf("GaussSeidel: max. no. of iterations exceeded !\n");
}