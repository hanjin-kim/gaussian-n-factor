#include <calibrator/processes/generalhullwhiteprocess.hpp>

#include <boost/bind.hpp>

#include <ql/compounding.hpp>

namespace HJCALIBRATOR
{
	GeneralizedHullWhiteProcess::GeneralizedHullWhiteProcess( const Handle<YieldTermStructure>& h,
															  const IntegrableParameter& a,
															  const Parameter& sigma )
		: termStructure_( h ), a_( a ), sigma_( sigma )
		, r0_( h->forwardRate( 0.0, 0.0, Continuous, NoFrequency ) )
		, integrator_( SimpsonIntegral( 0.00001, 1000 ) )
		, Vrintegrand_( boost::bind( &GeneralizedHullWhiteProcess::VrIntegrand, this, _1 ) )
		, OneOverEintegrand_( boost::bind( &GeneralizedHullWhiteProcess::E, this, _1, -1. ) )
	{}
}