/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2015
 */

#include "solver.hpp"
#include <tools/interpolation.hpp>
#include "interpolation2d.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h> // odeiv2 Requires GSL 1.15
#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include "vector.hpp"


/// Accuracy parameters

int THETAINTPOINTS = 200;
double THETAINTACCURACY = 0.00000001;
int RINTPOINTS = 200;
REAL RINTACCURACY = 0.00000001;
double DESOLVEACCURACY = 0.00000001; // orig 0.01, relative accuracy of de solver
double DESOLVERABSACCURACY = 0;

const int PHI_RB_AVERAGES=5;

using Amplitude::SQR;

/*
const int THETAINTPOINTS = 10;
const double THETAINTACCURACY = 0.005;
const int RINTPOINTS = 40;
const REAL RINTACCURACY = 0.005;
const double DESOLVEACCURACY = 0.005; // orig 0.01, relative accuracy of de solver
const double DESOLVERABSACCURACY = 0.00000001;
*/

Solver::Solver(AmplitudeR* N_, bool fast_solver)
{
    N=N_;
    deltay=0.2;
    bfkl=false;
    dndy=false;

    fast=fast_solver;
    if (fast)
    {
        THETAINTPOINTS = 8;
        THETAINTACCURACY = 0.05;
        RINTPOINTS = 8;
        RINTACCURACY = 0.05;
        DESOLVERABSACCURACY = 0;
        DESOLVEACCURACY = 0.45;

    }
}

/*
 * Solve the BK
 */

struct EvolutionHelperR
{
    AmplitudeR* N;
    Solver* S;
};

int EvolveR(REAL y, const REAL amplitude[], REAL dydt[], void *params);

void Solver::Solve(REAL maxy)
{
    /*
     * In order to be able to use Runge Kutta we create one large array
     * which is basically just vector that we evolve
     * Array size is RPoints()*BPoints()*ThetaPoints(), and indexes are
     * ind = rind*BPoitns()*ThetaPoints()+bind*ThetaPoints()+Theta
     *
     * Values are ln of amplitude
     */

    int vecsize = N->BPoints()*N->RPoints()*N->ThetaPoints();
    REAL *vec = new REAL[2*vecsize]; // First 1/2 is N, second 1/2 is v2
    for (int rind=0; rind<N->RPoints(); rind++)
    {
        for (int bind=0; bind<N->BPoints(); bind++)
        {
            for (int thetaind=0; thetaind < N->ThetaPoints(); thetaind++)
            {
                 vec[rind*N->BPoints()*N->ThetaPoints() + bind*N->ThetaPoints()+thetaind]
                    = N->Ntable(0, rind, bind, thetaind);
                
                vec[vecsize+rind*N->BPoints()*N->ThetaPoints() + bind*N->ThetaPoints()+thetaind]
                   = N->V2table(0, rind, bind, thetaind);
                
                if (thetaind > 0)
                {
                    std::cerr << "Note, this BBK code evolves N and v2, so it can't support thetaind grid!" << endl;
                    exit(1);
                }
            }
        }
    }

    REAL y=0; REAL step = deltay;  // We have always solved up to y

    // Intialize GSL
    EvolutionHelperR help; help.N=N; help.S=this;
    gsl_odeiv_system sys = {EvolveR, NULL, 2*vecsize, &help};
        
    gsl_odeiv_step_type * T;
    if (fast==false)
        T = (gsl_odeiv_step_type *)gsl_odeiv_step_rkf45; //2; //f45;
    else
        T = (gsl_odeiv_step_type *)gsl_odeiv_step_rk2;
    

    gsl_odeiv_step * s    = gsl_odeiv_step_alloc (T, 2*vecsize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (DESOLVERABSACCURACY, DESOLVEACCURACY);    //abserr relerr
    gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (2*vecsize);
    REAL h = step;  // Initial ODE solver step size
    
    do
    {
        REAL nexty = y+step;
        while (y<nexty)
        {
            int status = gsl_odeiv_evolve_apply(e, c, s, &sys,
                &y, nexty, &h, vec);
            if (status != GSL_SUCCESS) {
                cerr << "Error in gsl_odeiv_evolve_apply at " << LINEINFO
                << ": " << gsl_strerror(status) << " (" << status << ")"
                << " y=" << y << ", h=" << h << endl;
            }
            cout << "Evolved up to " << y << "/" << nexty << ", h=" << h << endl;
        }

        
        int yind = N->AddRapidity(nexty);
        for (int rind=0; rind<N->RPoints(); rind++)
        {
            for (int bind=0; bind<N->BPoints(); bind++)
            {
                for (int thetaind=0; thetaind < N->ThetaPoints(); thetaind++)
                {
                    int tmpind = rind*N->BPoints()*N->ThetaPoints()
                            + bind*N->ThetaPoints()+thetaind;
                    N->AddDataPoint(yind, rind, bind, thetaind, vec[tmpind]);
                    N->AddV2DataPoint(yind, rind, bind, thetaind, vec[vecsize+tmpind]);
                }
            }
        }
    } while (y <= maxy);

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    

    delete[] vec;


}

int EvolveR(REAL y, const REAL amplitude[], REAL dydt[], void *params)
{
    //cout << "Evolving with y=" <<y << endl;
    EvolutionHelperR* par = (EvolutionHelperR*)params;
    int vecsize = par->N->BPoints()*par->N->RPoints()*par->N->ThetaPoints();

    // Intialize 2d r,b interpolator
    std::vector<double> ndata;
    std::vector<double> v2data;
    for (int ri=0; ri < par->N->RPoints(); ri++)
    {
        for (int bi=0; bi < par->N->BPoints(); bi++)
        {
            int tmpind = ri*par->N->BPoints()+ bi;
            ndata.push_back(amplitude[tmpind]);
            v2data.push_back(amplitude[vecsize+tmpind]);
        }
    }
    // interp(b, ln r)
    // Todo: local for each thread, or global here?
    DipoleInterpolator2D dipinterp( par->N->BVals(), par->N->LogRVals(),ndata);
    DipoleInterpolator2D v2interp( par->N->BVals(), par->N->LogRVals(),v2data);
    
    
    int thetaind=0;
    #pragma omp parallel for collapse(2) //schedule(dynamic,2) // firstprivate(interp)
    for (int rind=0; rind < par->N->RPoints(); rind++)
    {
        //for (int thetaind=0; thetaind < par->N->ThetaPoints(); thetaind++)
        //{
            for (int bind=0; bind < par->N->BPoints(); bind++)
            {
                
                
                int tmpind = rind*par->N->BPoints()*par->N->ThetaPoints()
                    + bind*par->N->ThetaPoints()+thetaind;

                REAL tmplnr = par->N->LogRVal(rind);
                REAL tmpb = par->N->BVal(bind);
                REAL tmptheta = 0; // not used par->N->ThetaVal(thetaind);
                
               // if ( std::abs(std::exp(tmplnr)-1.9)>0.1 or std::abs(tmpb-2.2)>0.1) continue;


                dydt[tmpind] = par->S->RapidityDerivative(y, tmplnr, tmpb, tmptheta,
                    amplitude, &dipinterp, &v2interp, false);
                
                dydt[vecsize+tmpind] = par->S->RapidityDerivative(y, tmplnr, tmpb, tmptheta,
                    amplitude, &dipinterp, &v2interp, true);
                
#pragma omp parallel critical
                {
                    if (par->S->PrintDnDy() )
                    {
                        cout << std::exp(tmplnr) << " " << tmpb << " " << dydt[tmpind] << " " << dydt[vecsize+tmpind] << " " << amplitude[tmpind] << endl;
                    }
                }
                
                
                
                //cout << "dydt at r=" << std::exp(tmplnr)<< ", b=" << tmpb  << " dN0/dy " << dydt[tmpind] << " dv2/dy " << dydt[vecsize+tmpind] << " amp " << amplitude[tmpind]<< " v2 " << amplitude[vecsize+tmpind] << endl;
                //dipinterp.Evaluate(tmpb, tmplnr) << endl;

                
            }
        
       // }
    }
    
    if (par->S->PrintDnDy() )
    {
        exit(1);
    }
    
    
    return GSL_SUCCESS;
    
}

/*
 * Calculate \partial_Y N
 * Data contais the previously calculated data points at rapidity y, format
 * ind = rind*BPoitns()*ThetaPoints()+bind*ThetaPoints()+Theta
 */
struct Inthelper_rthetaint
{
    AmplitudeR* N;
    Solver* Solv;
    REAL y;
    REAL b01;   // impact parameter of the parent dipole
    REAL lnr01;   // parent dipole
    REAL lnr02;   // integrated over
    REAL n01;
    REAL n02;
    REAL thetab; // angle between r01 and b01
    REAL alphas_r01;
    REAL alphas_r02;
    const REAL* data;
    bool v2_comp; // True if we evolve modulation
    
    DipoleInterpolator2D *dipoleinterp;
    DipoleInterpolator2D *v2interp;
};
REAL Inthelperf_rint(REAL lnr, void* p);
REAL Inthelperf_thetaint(REAL theta, void* p);

REAL Solver::RapidityDerivative(REAL y,
            REAL lnr01, REAL b01, REAL thetab, const REAL* data,
            DipoleInterpolator2D *interp, DipoleInterpolator2D *v2interp, bool v2comp)
{

    //if (lnr01 <= N->LogRVal(0)) lnr01*=0.999;
    //else if (lnr01 >= N->LogRVal(N->RPoints()-1)) lnr01*=0.999;
    
    // Integrate first over r, then over \theta
    Inthelper_rthetaint helper;
    helper.N=N; helper.Solv=this;
    helper.y=y; helper.b01=b01; helper.lnr01=lnr01;
    helper.thetab = thetab; helper.data=data;
    helper.dipoleinterp = interp;
    helper.v2interp = v2interp;
    helper.v2_comp=v2comp;
    

    if (rc!=CONSTANT)
        helper.alphas_r01 = N->Alpha_s_ic(std::exp(2.0*lnr01));
    else
        helper.alphas_r01=0;

    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_rint;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);

    REAL minlnr = std::log( 0.5*N->RVal(0) );
    REAL maxlnr = std::log( 2.0*N->RVal(N->RPoints()-1) );

    int status; REAL result, abserr;
    status=gsl_integration_qag(&fun, minlnr,
            maxlnr, 0, RINTACCURACY, RINTPOINTS,
            GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        #pragma omp critical
        cerr << "RInt failed at " << LINEINFO << ": lnr=" << lnr01 <<", r= " << std::exp(lnr01) <<", thetab="
        << thetab <<", b01=" << b01 <<", result " << result << ", relerr "
        << std::abs(abserr/result) << endl;
    }

    return result;
}

REAL Inthelperf_rint(REAL lnr, void* p)
{
    Inthelper_rthetaint* par = (Inthelper_rthetaint*)p;


    par->lnr02=lnr;
    gsl_function fun;
    fun.function = Inthelperf_thetaint;
    fun.params = par;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

    REAL mintheta = 0.0000001;
    REAL maxtheta = 2.0*M_PI-0.0001;

    if (par->Solv->GetRunningCoupling()!=CONSTANT)
         par->alphas_r02 = par->N->Alpha_s_ic(std::exp(2.0*lnr));
    else
        par->alphas_r02=0;


    int status; REAL result, abserr;
    status = gsl_integration_qag(&fun, mintheta,
            maxtheta, 0, THETAINTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);
    if (status and std::abs(result)>1e-3) 
    {
        if (! (std::exp(par->lnr01) > 0.1 and std::exp(lnr) < 1e-5) )
            #pragma omp critical
            cerr << "Thetaint failed at " << LINEINFO <<": r=" << std::exp(lnr) 
                << " n01 " << par->n01 << " r01 " << std::exp(par->lnr01) <<
                " result " << result << " relerr " << std::abs(abserr/result) << endl;
    }

    result *= std::exp(2.0*lnr);        // Change of variables r->ln r

    return result;

}

REAL Inthelperf_thetaint(REAL theta, void* p)
{
    Inthelper_rthetaint* par = (Inthelper_rthetaint*)p;
    
    Vec b(par->b01, 0);
    double rlen = std::exp(par->lnr01);
    Vec x2(std::exp(par->lnr02)*cos(theta), std::exp(par->lnr02)*sin(theta));
    
    double result=0;
    

    double sum1=0;
    double sum2=0;
    double sum3=0;
    double sum4=0;
    for (int avg=0; avg < PHI_RB_AVERAGES; avg++)
    {
        // Todo: integral would be more accurate...
        double phirb = 2.0*M_PI*avg/(PHI_RB_AVERAGES-1);
        Vec r(rlen*std::cos(phirb), rlen*std::sin(phirb));
        
        Vec x0 = (r+b*2.0)*0.5;
        Vec x1 = (b*2.0-r)*0.5;
        
        
        
        Vec b1 = (x0+x2)*0.5;
        Vec r1 = x0-x2;
        double cos_b1_r1 = b1*r1/(b1.Len()*r1.Len());
        double cos2phi_b1_r1 = 2.0*pow(cos_b1_r1,2)-1.; // cos(2x) = 2cos^2(x)-1
        
        sum1 = par->dipoleinterp->Evaluate(b1.Len(), std::log(r1.Len()) )
            + 2.0*par->v2interp->Evaluate(b1.Len(), std::log(r1.Len())) * cos2phi_b1_r1;
        
        
        
        Vec b2 =(x1+x2)*0.5;
        Vec r2 = x1-x2;
        double cos_b2_r2 = b2*r2/(b2.Len()*r2.Len());
        double cos2phi_b2_r2 = 2.0*pow(cos_b2_r2,2)-1.; // cos(2x) = 2cos^2(x)-1
        
        sum2 = par->dipoleinterp->Evaluate(b2.Len(),std::log(r2.Len()) )
        + 2.0*par->v2interp->Evaluate(b2.Len(), std::log(r2.Len())) *  cos2phi_b2_r2;
        
        sum3 = par->dipoleinterp->Evaluate( par->b01, par->lnr01)
         + 2.0*par->v2interp->Evaluate(par->b01,par->lnr01) * std::cos(2.0 * phirb);
        // Note: v2 part above cancels only in <N0> component, not when computing v2 evolution
       
        
        sum4 =(par->dipoleinterp->Evaluate(b1.Len(),std::log(r1.Len()))
               + 2.0*par->v2interp->Evaluate(b1.Len(), std::log(r1.Len())) * cos2phi_b1_r1)
               *( par->dipoleinterp->Evaluate( b2.Len(), std::log(r2.Len()))
                 +2.0*par->v2interp->Evaluate(b2.Len(), std::log(r2.Len())) * cos2phi_b2_r2 );
        
        if (isnan(sum1) or isnan(sum2) or isnan(sum3) or isnan(sum4))
        {
            cout << "NAN! "<< endl;
            cout <<cos2phi_b1_r1 << " " <<cos2phi_b2_r2 << endl;
            cout << r1 << endl << r2 << endl << b1 << endl << b2 << endl; exit(1);
        }
        
        if (par->v2_comp)
        {
            sum1 *= std::cos(2.0*phirb);
            sum2 *= std::cos(2.0*phirb);
            sum3 *= std::cos(2.0*phirb);
            sum4 *= std::cos(2.0*phirb);
        }
        
        // Kernels
        double Xlen = (x0-x2).Len();
        double Ylen = (x1-x2).Len();
        result += (sum1+sum2-sum3-sum4)*par->Solv->Kernel(rlen, Xlen, Ylen, par->alphas_r01, par->alphas_r02, par->N->Alpha_s_ic(Ylen*Ylen),par->y, theta);
        
    }
    result /= PHI_RB_AVERAGES;
    return result;


    REAL r01 = std::exp(par->lnr01);
    REAL r02 = std::exp(par->lnr02);
    REAL r12sqr = SQR(r01)+SQR(r02)-2.0*r01*r02*std::cos(theta);

    REAL n12=0;
        
    REAL alphas_r12=0;
    if (par->Solv->GetRunningCoupling()!=CONSTANT
        and par->Solv->GetRunningCoupling() != PARENT)
                alphas_r12 = par->N->Alpha_s_ic(r12sqr);

    
    result *= par->Solv->Kernel(r01, r02, std::sqrt(r12sqr), par->alphas_r01,
            par->alphas_r02, alphas_r12, par->y, theta);

        return result;


    return 0;
}

/*
 * BK kernel
 * r01: parent dipole size, b01: parent dipole impact parameter,
 * thetab: angle between vecs r01 and b01, r02 integration variable, pos. of
 * point 2, theta2 is angle between r01 and r02, also integrated over
 * Rapidity y may be used with modified kernel
 * By default y=b01=tetab=theta=0
 *
 * Running coupling is also handled here
 */
REAL Solver::Kernel(REAL r01, REAL r02, REAL r12, REAL alphas_r01,
        REAL alphas_r02, REAL alphas_r12, REAL y,
        REAL theta2, REAL b01, REAL thetab)
{
    if (r12<N->MinR()/10 or r02 < N->MinR()/10)
        return 0;
    REAL result=0;

    // Ref for different prescriptions: 0704.0612
    // Convention: r01=r, r02 = r1, r12=r2
    REAL Rsqr, z, zsqrt;
    REAL costheta2; REAL r02dotr12;

    switch(rc)
    {
        case CONSTANT:
            result = ALPHABAR_s/(2.0*M_PI)
                    * SQR(r01) / ( 1e-20+ SQR(r12) * SQR(r02) );
			if (isnan(result) or isinf(result)) return 0;
			else return result;
            break;
        case PARENT:
            result = alphas_r01*Nc/(2.0*SQR(M_PI))
                    * SQR(r01) / ( SQR(r12) * SQR(r02) );
            break;
        case BALITSKY:
            result = Nc/(2.0*SQR(M_PI))*alphas_r01
            * (
            SQR(r01) / ( SQR(r12) * SQR(r02) )
            + 1.0/SQR(r02)*(alphas_r02/alphas_r12 - 1.0)
            + 1.0/SQR(r12)*(alphas_r12/alphas_r02 - 1.0)
            );
			if (isnan(result) or isinf(result)) return 0;
			else return result;
			return result;
            break;
        case JIMWLK_SQRTALPHA:  // Similar than BALITSKY, slightly simpler (TL&HM JIMWLK paper)
            result = Nc/(2.0*SQR(M_PI))
            * (
            alphas_r02/SQR(r02) + alphas_r12/SQR(r12)
				+ std::sqrt(alphas_r02*alphas_r12)* ( SQR(r01/(r02*r12)) - 1.0/SQR(r02) - 1.0/SQR(r12)) 	// = -2 (x-z).(y-z)/(x-z)^2(y-z)^2
            );
            break;
        case KW:
            if (std::abs(SQR(r02)-SQR(r12)) < SQR(N->MinR())/1000) return 0;

            costheta2 = std::cos(theta2);
            r02dotr12 = SQR(r02) - r02*r01*costheta2;

            /* calculate dot product using geometry, slower than vector
             * calculateion used above
             if (theta2>M_PI) theta2=2.0*M_PI-theta2;

            theta012 = std::acos(
                -0.9999*( SQR(r02) - SQR(r01) - SQR(r12) ) / (2.0*r01*r12) );

            cosr1r2 = -std::cos(theta2 + theta012);
            r02dotr12 = r02*r12*cosr1r2;   // r02 \cdot r12
            */
            

            Rsqr = r02*r12*std::pow((r12/r02),
                (SQR(r02)+SQR(r12))/(SQR(r02)-SQR(r12))
                    - 2.0*SQR(r02)*SQR(r12)/(r02dotr12*(SQR(r02)-SQR(r12)) )
                );
            /*if (std::abs(r02-r12) < 1e-3 and std::abs(r02-r12)>1e-5)
                cout << "r1 " << r02 << " r2 " << r12 << " theta " << theta2 << " Rsqr " <<
                Rsqr << " exp " << r02*r12*std::exp(-1.0+1.0/costheta2) << endl;*/
            if (Rsqr<1e-50) return 0;

            result = Nc/(2.0*SQR(M_PI)) * (
                alphas_r02 / SQR(r02)
                - 2.0*alphas_r02*alphas_r12 / N->Alpha_s_ic(Rsqr)
                 * r02dotr12 / (SQR(r02)*SQR(r12) )
                + alphas_r12 / SQR(r12)
                );
            break;
        case MS:
            // Motyka & Staśto, 0901.4949: kinematical constraint, bessel kernel
            ///FIXME: 0.01 factor????
            cerr << "Motyka&Stasto kin. constraint kernel is most likely wrong..."
				<<" Check before you use this! " << LINEINFO << endl;
				
            z = 0.01*std::exp(-y); zsqrt = std::exp(-0.5*y);

            // r02 dot r12
            costheta2 = std::cos(theta2);
            r02dotr12 = SQR(r02) - r02*r01*costheta2;

            REAL besselr02r01 = gsl_sf_bessel_K1(r02/r01*zsqrt);
            REAL besselr12r01 = gsl_sf_bessel_K1(r12/r01*zsqrt);

            result = z/SQR(r01) * (
                  SQR( besselr02r01 )
                + SQR( besselr12r01 )
                -  2.0*besselr02r01*besselr12r01
                  * r02dotr12 / (r02*r12)
                );

            // Parent dipole RC
            result *= Nc/(2.0*SQR(M_PI))*alphas_r01;
            break;
            
    }
    return result;
        
}

/*
 * Interpolate data array
 */
REAL Solver::InterpolateN(REAL lnr, REAL lnb, REAL thetab, const REAL* data)
{
    cerr << "Solver::InterpolateN called, why?" << endl;
}


void Solver::SetRunningCoupling(RunningCoupling rc_)
{
    rc=rc_;
}

RunningCoupling Solver::GetRunningCoupling()
{
    return rc;
}

void Solver::SetDeltaY(REAL dy)
{
    deltay = dy;
}

void Solver::SetBfkl(bool bfkl_)
{
    bfkl=bfkl_;
}

bool Solver::GetBfkl()
{
    return bfkl;
}
