/*
 * Interpoate 2D data array
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _INTERPOLATION2D_H
#define _INTERPOLATION2D_H



#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <vector>


/**
 * 2D interpolator
 * Interpolates given data using spline (goes trough every data point)
 * Uses GSL
 *
 * Xdata and ydata pointers are saved, but not used for interpolation purposes
 * If user frees the allocated memory, one should be sure that these pointers
 * are not asked from this class!
 *
 * Uses Interpolator class to perform 1D interpolations
 */

class DipoleInterpolator2D
{
    public:
        // data format [xind][yind]
        // Grid is rectangular, that is, grid(x,y) = (grid[x], grid[y])
        /**
         * Initialize 2D interpolator
         * zgrid[i] is the point at xgrid[i],ygrid[i]
         */
        DipoleInterpolator2D(std::vector<double> xgrid, std::vector<double> ygrid,
                       std::vector<double> zgrid);
                       
        ~DipoleInterpolator2D();
        
        double Evaluate(double x, double y);
    
        double MinX() { return minx; };
        double MinY() { return miny; };
        double MaxX() { return maxx; };
        double MaxY() { return maxy; };

    private:
        gsl_spline2d *gslinterp;
        gsl_interp_accel *xacc;
        gsl_interp_accel *yacc;
        double minx,maxx,miny,maxy;


};




#endif
