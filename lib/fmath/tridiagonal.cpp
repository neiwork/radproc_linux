#include "tridiagonal.h"
#include "matrixInit.h"
#include <math.h>

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
	for (j=0;j<n;j++) {						// loop over columns
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
			for (k=0;k<n;k++)
				{ t = a[imax][k]; a[imax][k] = a[j][k]; a[j][k] = t; }
		}
		det *= a[j][j];							// multiply determinant with pivot
		t = 1e0/a[j][j];
		for (i=j+1;i<n;i++) a[i][j] *= t;	// divide elements of L by pivot
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

void MatMul(Matrix a, Matrix& b, size_t m)
{
	Matrix c;
	matrixInitCopy(c,m,m,b);
	for (size_t i=0;i<m;i++) {
		for (size_t j=0;j<m;j++) {
			double sum = 0.0;
			for (size_t k=0;k<m;k++)
				sum += a[i][k] * c[k][j];
			b[i][j] = sum;
		}
	}
}

#define SQR(a) ((a)*(a))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))


double pythag(double a, double b)
// Computes (a 2 + b 2 ) 1/2 without destructive underflow or overflow.
{
	double absa,absb;
	
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa*(sqrt(1.0+SQR(absb/absa)));
	else return (absb == 0.0 ? 0.0 : absb*(sqrt(1.0+SQR(absa/absb))));
}

void svdcmp(Matrix& a, size_t m, size_t n, Vector& w, Matrix& v)
// Given a matrix a[0..m-1][0..n-1] , this routine computes its singular value decomposition,
// A = U·W·Vt. The matrix U replaces a on output. The diagonal matrix of singular values W is
// output as a vector w[0..n-1]. The matrix V (not the transpose Vt) is output as v[0..n-1][0..n-1].
{
	double pythag( double a, double b);
	int flag, i, its, j, jj, k, l, nm;
	double anorm, c, f, g, h, s, scale, x, y, z;
	Vector rv1(n,0.0);
	
	g = scale = anorm = 0.0;							// Householder reduction to bidiagonal form.

	for (i=0;i<n;i++) {
		
		l = i+1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s),f);
				h = f*g - s;
				a[i][i] = f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f = s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		
		w[i] = scale *g;
		g = s = scale = 0.0;
		if (i < m && i != n-1) {
			for (k=l;k<n;k++) scale += (fabs(a[i][k]));
			if (scale) {
				for (k=l;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s),f);
				h = f*g - s;
				a[i][l] = f-g;
				for (k=l;k<n;k++) rv1[k] = a[i][k]/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	
	for (i=n-1;i>=0;i--) {							// Accumulation of right-hand transformations.
		if (i < n-1) {
			if (g) {
				for (j=l;j<n;j++)					// Double division to avoid possible underflow.
					v[j][i] = (a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}

	for (i=IMIN(m,n)-1;i>=0;i--) {					// Accumulation of left-hand transformations.
		l = i+1;
		g =w[i];
		for (j=l;j<n;j++) a[i][j] = 0.0;
		if (g) {
			g = 1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f = (s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i] = 0.0;
		++a[i][i];
	}

	for (k=n-1;k>=0;k--) {							// Diagonalization of the bidiagonal form: Loop over
		for (its=1;its<=30;its++) {					// singular values, and over allowed iterations.
			flag=1;
			for (l=k;l>=0;l--) {					// Test for splitting.
				nm = l-1;							// Note that rv1[0] is always zero.
				if ((fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((fabs(w[nm])+anorm) == anorm) break;
			}
			
			if (flag) {
				c = 0.0;							// Cancellation of rv1[l], if l > 1.
				s = 1.0;
				for (i=l;i<=k;i++) {
					f = s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g = w[i];
					h = pythag(f,g);
					w[i] = h;
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y*c+z*s;
						a[j][i] = z*c-y*s;
					}
				}
			}
			z = w[k];
			if (l == k) {									// Convergence.
				if (z < 0.0) {								// Singular value is made nonnegative.
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) printf("no convergence in 30 svdcmp iterations");
			x = w[l];										// Shift from bottom 2-by-2 minor.
			nm = k-1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c = s = 1.0;									// Next QR transformation:
			for (j=l;j<=nm;j++) {
				i = j+1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x*c+z*s;
					v[jj][i] = z*c-x*s;
				}
				z = pythag(f,h);
				w[j] = z;								// Rotation can be arbitrary if z = 0.
				if (z) {
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				
				for (jj=0;jj<m;jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y*c+z*s;
					a[jj][i] = z*c-y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
}

void Transpose(Matrix& a, size_t M) {
	
	Matrix aaux;
	matrixInitCopy(aaux,M,M,a);
	for (size_t i=0;i<M;i++)
		for (size_t j=0;j<M;j++)
			a[i][j] = aaux[j][i];
}

//===========================================================================
void MatInv(Matrix& a, int n, double &det, int& err)
//---------------------------------------------------------------------------
// Calculates inverse of real matrix by LU factorization.
// a - (n x n) matrix (input); a^(-1) (output)
// det - determinant of coefficient matrix (output).
// Calls: LUFactor, LUSystem.
//---------------------------------------------------------------------------
{
	err = 0;
	Matrix ainv;
	Vector b(n,0.0);
	int i, j;
	vector<int> ipivot(n,0.0);
	matrixInit(ainv,n,n,0.0);
	LUFactor(a,ipivot,n,det);			// LU factorization of a
	
	if (det == 0e0) {						// singular matrix
		printf("MatInv: singular matrix !\n");
		/*rr = 1;
		Matrix v;			matrixInit(v,n,n,0.0);
		Vector w(n,0.0);
		Matrix aaux;		matrixInitCopy(aaux,n,n,a);
		svdcmp(aaux,n,n,w,v);
		double max, min;
		max = min = w[0];
		for (size_t i=1;i<n;i++) {
			if (w[i] > max) max = w[i];
			if (w[i] < min) min = w[i];
		}
		double cond = max/min;
		cout << "inverse of the condition number = " << 1.0/cond << endl;
		if (fabs(1.0/cond) > 1e-12) {
			Transpose(aaux,n);
			for (size_t row=0;row<n;row++)
				for (size_t col=0;col<n;col++)
					a[row][col] = (1.0/w[row]) * aaux[row][col];
			MatMul(v,a,n);
		}*/
		return;
	}
	
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
	Vector y(J*M,0.0);
	
	matrixInit(beta,M,M,0.0);
	for (size_t row=0;row<M;row++) beta[row][row] = 1.0;
	for (size_t m=0;m<M;m++) y[m] = d[m];
	
	double det = 1.0;
	int err = 0;
	for (size_t j=1;j<J;j++) {
		
		matrixInit(alpha,M,M,0.0);
		if (j>1) {
			MatInv(beta,M,det,err);
		}
		
		for (size_t col=0;col<M;col++)
			for (size_t row=0;row<M;row++)
				alpha[row][col] = beta[row][col] * a[j*M+col];
		
		matrixInit(beta,M,M,0.0);
		for (size_t col=0;col<M;col++) {
			double sum = 0.0;
			for (size_t row=0;row<M;row++) {
				if (row == col)
					beta[row][col] = bb[j*M+row] - alpha[row][row] * c[(j-1)*M+col];
				else if (row == col+1)
					beta[row][col] = ba[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
				else if (row == col-1)
					beta[row][col] = bc[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
				else
					beta[row][col] = - alpha[row][col] * c[(j-1)*M+col];
					
				sum += ( alpha[col][row] * y[(j-1)*M+row] );
			}
			y[j*M+col] = d[j*M+col] - sum;
		}
	}
	
	err = 0;
	MatInv(beta,M,det,err);
	
	for (size_t m=0;m<M;m++) {
		double sum = 0.0;
		for (size_t mm=0;mm<M;mm++)
			sum += beta[m][mm] * y[(J-1)*M+mm];
		x[(J-1)*M+m] = sum;
	}
	
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
			err = 0;
			MatInv(beta,M,det,err);
			
			for (size_t col=0;col<M;col++)
				for (size_t row=0;row<M;row++)
					alpha[row][col] = beta[row][col] * a[j*M+col];
			
			for (size_t col=0;col<M;col++) {
				double sum = 0.0;
				for (size_t row=0;row<M;row++) {
					if (row == col)
						beta[row][col] = bb[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
					else if (row == col+1)
						beta[row][col] = ba[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
					else if (row == col-1)
						beta[row][col] = bc[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
					else
						beta[row][col] = - alpha[row][col] * c[(j-1)*M+col];
					
					sum += ( alpha[col][row] * y[(j-1)*M+row] );
				}
				y[j*M+col] = d[j*M+col] - sum;
			}
		}
		
		err = 0;
		MatInv(beta,M,det,err);
		
		for (size_t m=0;m<M;m++) {
			double sum = 0.0;
			for (size_t mm=0;mm<M;mm++) {
				double sum_2 = c[(U-1)*M+mm] * x[U*M+mm];
				sum += ( beta[m][mm] * (y[(U-1)*M+mm] - sum_2) );
			}
			x[(U-1)*M+m] = sum;
		}
	}
	for (size_t i=0;i<J*M;i++) d[i] = x[i];
}


//===========================================================================
void TriBlockDiagSys2(Vector a, Vector ba, Vector bb, Vector bc, Vector c, Vector& d, size_t J, size_t M)
//---------------------------------------------------------------------------
// Solves a system with block tridiagonal matrix by LU factorization (diag(L) = 1).
// a - lower codiagonal (i=1,n-1)
// b - main diagonal (i=0,n-1)
// c - upper codiagonal (i=0,n-2)
// d - constant terms (i=0,n-1); solution on exit
// n - order of system.
//---------------------------------------------------------------------------
{
	Vector betaVec(J*M*M,0.0);
	Matrix alpha, beta;
	matrixInit(alpha,M,M,0.0);
	matrixInit(beta,M,M,0.0);
	Vector x(J*M,0.0);
	Vector y(J*M,0.0);
	
	for (size_t m=0;m<M;m++) {
		betaVec[m*M+m] = 1.0;
		y[m] = d[m];
		beta[m][m] = 1.0;
	}
	
	double det = 1.0;
	int err = 0;
	for (size_t j=1;j<J;j++) {
		
		if (j>1) {
			MatInv(beta,M,det,err);
		}
		
		for (size_t col=0;col<M;col++)
			for (size_t row=0;row<M;row++)
				alpha[row][col] = beta[row][col] * a[j*M+col];
		
		matrixInit(beta,M,M,0.0);
		for (size_t col=0;col<M;col++) {
			double sum = 0.0;
			for (size_t row=0;row<M;row++) {
				if (row == col)
					beta[row][col] = bb[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
				else if (row == col+1)
					beta[row][col] = ba[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
				else if (row == col-1)
					beta[row][col] = bc[j*M+row] - alpha[row][col] * c[(j-1)*M+col];
				else
					beta[row][col] = - alpha[row][col] * c[(j-1)*M+col];
				
				betaVec[j*M*M+row*M+col] = beta[row][col];
				sum += ( alpha[col][row] * y[(j-1)*M+row] );
			}
			y[j*M+col] = d[j*M+col] - sum;
		}
	}
	
	for (size_t row=0;row<M;row++)
		for (size_t col=0;col<M;col++)
			beta[row][col] = betaVec[(J-1)*M*M+row*M+col];
	MatInv(beta,M,det,err);
	
	for (size_t m=0;m<M;m++) {
		double sum = 0.0;
		for (size_t mm=0;mm<M;mm++) {
			sum += ( beta[m][mm] * y[(J-1)*M+mm] );
		}
		x[(J-1)*M+m] = sum;
	}
	
	for (size_t j=J-1;j>=1;j--) {
		for (size_t row=0;row<M;row++)
			for (size_t col=0;col<M;col++)
				beta[row][col] = betaVec[(j*M+row)*M+col];
		MatInv(beta,M,det,err);
		
		for (size_t m=0;m<M;m++) {
			double sum = 0.0;
			for (size_t mm=0;mm<M;mm++) {
				double sum_2 = c[(j-1)*M+mm] * x[j*M+mm];
				sum += ( beta[m][mm] * (y[(j-1)*M+mm] - sum_2) );
			}
			x[(j-1)*M+m] = sum;
		}
	}
	for (size_t i=0;i<J*M;i++) d[i] = x[i];
}

void Jacobi(Vector a, Vector ba, Vector bb, Vector bc, Vector c, Vector& d, size_t J, size_t M)
{
	Matrix LplusU;
	matrixInit(LplusU,(J+1)*(M+1),(J+1)*(M+1),0.0);
	
	for (size_t j1=0;j1<J+1;j1++) {
		for (size_t j2=0;j2<J+1;j2++) {
			for (size_t m1=0;m1<M+1;m1++) {
				for (size_t m2=0;m2<M+1;m2++) {
					size_t i1 = j1*(M+1)+m1;
					size_t i2 = j2*(M+1)+m2;
					if (j1 == j2) {					// B matrix
						if (m2 == m1-1) LplusU[i1][i2] = ba[j1*(M+1)+m1];
						else if (m2 == m1+1) LplusU[i1][i2] = bc[j1*(M+1)+m1];
						else LplusU[i1][i2] = 0.0;
					} else if (j2 == j1-1) {		// A matrix
						if (m2 == m1) LplusU[i1][i2] = a[j1*(M+1)+m1];
						else LplusU[i1][i2] = 0.0;
					} else if (j2 == j1+1) {		// C matrix
						if (m2 == m1) LplusU[i1][i2] = c[j1*(M+1)+m1];
						else LplusU[i1][i2] = 0.0;
					} else
						LplusU[i1][i2] = 0.0;
				}
			}
		}
	}
	
	for (int it=1;it<10;it++) {
		for (size_t j=0;j<J;j++) {
			for (size_t m=0;m<M;m++) {
				size_t l = j*M+m;
				double sum1 = 0.0;
				for (size_t k=0;k<M*J;k++) {
					sum1 += LplusU[l][k] * d[k];
				}
				d[l] = (1.0/bb[l]) * (d[l]-sum1);
			}
		}
	}
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
	const int itmax = 100;				// max no. of iterations
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