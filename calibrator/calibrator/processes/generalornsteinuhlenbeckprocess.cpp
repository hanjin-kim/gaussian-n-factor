#include <calibrator/processes/generalornsteinuhlenbeckprocess.hpp>

#include <boost/bind.hpp>

namespace HJCALIBRATOR
{
	GeneralizedOrnsteinUhlenbeckProcess::GeneralizedOrnsteinUhlenbeckProcess( const IntegrableParameter& a, // must be an IntegrableParameter
																			  const Parameter& sigma,
																			  const Real x0 )
		: x0_( x0 ), a_( a ), sigma_( sigma )
		, integrator_( GaussKronrodAdaptive( 0.00001, 1000 ) )
		, Vrintegrand_( boost::bind( &GeneralizedOrnsteinUhlenbeckProcess::VrIntegrand, this, _1 ) )
	{}
}