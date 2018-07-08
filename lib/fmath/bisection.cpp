#include "bisection.h"
#define NTRY 100
#define FACTOR 1.6

int zbrac2(fun1 func,double& x1,double& x2)
    /* Given a function func and an initial guessed range x1 to x2, the routine expands the range
    geometrically until a root is bracketed by the returned values x1 and x2 (in which case zbrac
    returns 1) or until the range becomes unacceptably large (in which case zbrac returns 0). */
{
    void nrerror(char error_text[]);
    int j;
    float f1,f2;
    if (x1 == x2) nrerror("Bad initial range in zbrac");
    f1=func(x1);
    f2=func(x2);
    for (j=1;j<=NTRY;j++) {
        if (f1*f2 < 0.0) return 1;
        if (fabs(f1) < fabs(f2))
            f1=func(x1 += FACTOR*(x1-x2));
        else
            f2=func(x2 += FACTOR*(x2-x1));
    }
    return 0;
}
#define JMAX 40

double bisection(fun1 func,double x1,double x2,double xacc)
    /* Using bisection, find the root of a function func known to lie between x1 and x2. The root,
    returned as rtbis, will be refined until its accuracy is ±xacc. */
{
    void nrerror(char error_text[]);
    int j;
    double dx,f,fmid,xmid,rtb;
    
    zbrac2(func,x1,x2);
    f=func(x1);
    fmid=func(x2);
    if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);            // Orient the search so that f>0
    for (j=1;j<=JMAX;j++) {                                                   // lies at x+dx.
        fmid=func(xmid=rtb+(dx *= 0.5));                             // Bisection loop.
        if (fmid <= 0.0) rtb=xmid;
        if (fabs(dx) < xacc || fmid == 0.0) return rtb;
    }
    nrerror("Too many bisections in rtbis");
    return 0.0;                                                                              // Never get here.
}