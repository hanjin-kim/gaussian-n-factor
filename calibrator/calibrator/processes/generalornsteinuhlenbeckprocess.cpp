#include <calibrator/processes/generalornsteinuhlenbeckprocess.hpp>

#include <boost/bind.hpp>

namespace HJCALIBRATOR
{
	GeneralizedOrnsteinUhlenbeckProcess::GeneralizedOrnsteinUhlenbeckProcess( const IntegrableParameter& a, // must be an IntegrableParameter
																			  const Parameter& sigma,
																			  const Real x0 )
		: x0_( x0 ), a_( a ), sigma_( sigma )
		, integrator_( SimpsonIntegral( 10, 10000 ) )
		//, Vrintegrand_( boost::bind( &GeneralizedOrnsteinUhlenbeckProcess::VrIntegrand, this, _1 ) )
		, int_wrkspcs_( gsl_integration_workspace_alloc(1000) )
	{
		boost::function<Real( Real )> tmpptr = boost::bind(&GeneralizedOrnsteinUhlenbeckProcess::VrIntegrand, this, _1);
		integrand_Vr_ = convertToGslFunction( tmpptr );
	}

	GeneralizedOrnsteinUhlenbeckProcess::~GeneralizedOrnsteinUhlenbeckProcess()
	{
		gsl_integration_workspace_free( int_wrkspcs_ );
	}
}