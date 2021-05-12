/*
 * 2D interpolation (slow)
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */


#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "interpolation2d.hpp"
#include <iostream>



double DipoleInterpolator2D::Evaluate(double x, double y)
{
	double n;
    
    if (x > maxx or y>maxy or x<minx or y<miny)
        return 0;
    
	int s = gsl_spline2d_eval_e(gslinterp,x,y, xacc, yacc, &n);
	if (s==GSL_EDOM) return 0;
	if (s) return 0;
    
    if (n<0 or n>1)
    {
        std::cerr << "What, N(b=" << x <<", r=" << std::exp(y) << " = " << n << std::endl;
    }
    return n; //gsl_spline2d_eval(gslinterp,x,y, xacc, yacc);


}


// Initialize, data format: [xind][yind]
// Note that we assume that xgrid and ygrid contain ONLY the
// distinct grid points
// z[x_i,y_j] = zgrid[j*xgridsize+i]
DipoleInterpolator2D::DipoleInterpolator2D(std::vector<double> xgrid,
                std::vector<double> ygrid, std::vector<double> zgrid)
{

    
    double* xdata = xgrid.data();
    double* ydata= ygrid.data();
    double* zdata = zgrid.data();
    
    const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    
    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();
    
    gslinterp = gsl_spline2d_alloc(T, xgrid.size(), ygrid.size());
    gsl_spline2d_init(gslinterp, xdata, ydata, zdata, xgrid.size(), ygrid.size());
    
    minx=xgrid[0]; maxx=xgrid[xgrid.size()-1];
    miny=ygrid[0]; maxy=ygrid[ygrid.size()-1];
}



DipoleInterpolator2D::~DipoleInterpolator2D()
{
    gsl_spline2d_free(gslinterp);
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);
}
