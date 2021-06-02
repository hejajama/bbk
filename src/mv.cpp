#include "mv.hpp"
#include <sstream>
#include <string>
#include <cmath>
#include <tools/config.hpp>
#include "config.hpp"

using Amplitude::SQR;

double MV::DipoleAmplitude(double r, double b)
{
    // Berger, Stasto
    return 1.0 - std::exp(-10*r*r*std::exp(-b*b/2));
    
    
    return 1.0-std::exp(-10*r*r*std::exp(-b*b/2.));
	
	const double e = 2.7182818;
    double Bp=4.;
    return 1. - std::exp(-std::pow(SQR(r)*qs0sqr, anomalous_dimension)/4.0*std::exp(-b*b/Bp)
            * std::log( 1.0/(r*lambdaqcd) + ec*e) );
    
}

double MV::DipoleAmplitudeV2(double r, double b)
{
    return 0; // Berger, Stasto
    double v2 = r*r*b*b/40.;
    if (v2/DipoleAmplitude(r,b) > 0.5) return 0.5*DipoleAmplitude(r, b);
    else return v2;
}

void MV::SetQsqr(double qsqr)
{
	qs0sqr=qsqr;
}

void MV::SetAnomalousDimension(double gamma_)
{
	anomalous_dimension=gamma_;
}

void MV::SetLambdaQcd(double lambda)
{
	lambdaqcd=lambda;
}

void MV::SetE(double ec_)
{
	ec=ec_;
}

double MV::GetE()
{
    return ec;
}

std::string MV::GetString()
{
	std::stringstream ss;
	ss << "MV model, Q_s0^2 = " << qs0sqr << " GeV^2, \\gamma = " << anomalous_dimension
		<< ", coefficient of E inside Log is " << ec 
		<< ", x0=" << x0 << ", \\Lambda_QCD = " << lambdaqcd << " GeV";
	return ss.str();
}

/*
 * Set some reasonable parameters
 * That is, MV1 model for nucleus
 */
MV::MV()
{
	qs0sqr = 0.72;
	x0=0.007;
	ec=1.0;
	lambdaqcd=LAMBDAQCD;
	anomalous_dimension=1;
}
