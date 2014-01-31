/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _SOLVER_H
#define _SOLVER_H

#include "amplitude.hpp"
#include <tools/interpolation.hpp>

enum RunningCoupling
{
    CONSTANT,
    PARENT,
    KW,
    BALITSKY,
    MS,  // Motyka & Staśto, 0901.4949: kinematical constraint, bessel kernel
    JIMWLK_SQRTALPHA
};

class Solver
{
    public:
        Solver(AmplitudeR* N_);
        void Solve(REAL maxy);

        REAL RapidityDerivative(REAL y, REAL lnr01, REAL lnb01, REAL thetab,
            const REAL* data, Interpolator* interp);

        REAL Kernel(REAL r01, REAL r02, REAL r12, REAL alphas_r01=0,
            REAL alphas_r02=0, REAL alphas_r12=0,
            REAL y=0, REAL theta2=0, REAL b01=0, REAL thetab=0 );
        REAL InterpolateN(REAL lnr, REAL lnb, REAL thetab, const REAL* data);

        void SetRunningCoupling(RunningCoupling rc_);
        RunningCoupling GetRunningCoupling();

        void SetDeltaY(REAL dy);

        void SetBfkl(bool bfkl_);
        bool GetBfkl();

    private:
        AmplitudeR* N;

        RunningCoupling rc;
        REAL deltay;
        bool bfkl;

        

};

const int INTERPOLATION_POINTS = 6;

const double ALPHABAR_s = 0.2;

#endif
