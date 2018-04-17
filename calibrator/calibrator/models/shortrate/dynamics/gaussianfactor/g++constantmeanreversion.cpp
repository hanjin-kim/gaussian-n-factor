#include <calibrator/models/shortrate/dynamics/gaussianfactor/g++constantmeanreversion.hpp>

namespace HJCALIBRATOR
{
	ParamVector&& RealVectorToParamVector( const RealVector & a, Constraint constraint )
	{
		ParamVector a_;

		for ( Real ai : a )
		{
			a_.push_back( ConstantParameter( ai, constraint ) );
		}

		return std::move( a_ );
	}

	Real GPPConstantMeanReversion::E( const Parameter & a, Time s, Time t ) const
	{
		return exp( a( 0.0 ) * (t - s) );
	}

	Real GPPConstantMeanReversion::B( const Parameter & a, Time s, Time t ) const
	{
		Real aval = a( 0.0 );
		return (1 - exp( -aval * (t - s) )) / aval;
	}

	Real GPPConstantMeanReversion::integralVariance( Size i, Size j, Time s, Time t ) const
	{
		Real a_i = a( i, 0.0 );
		Real a_j = a( j, 0.0 );

		const Parameter& sigma_i = sigma( i );
		const Parameter& sigma_j = sigma( j );

		auto integrand = [&, a_i, a_j, sigma_i, sigma_j, t]( Time u )
		{
			return sigma_i( u ) * sigma_j( u ) * (1 - exp( -a_i * (t - u) )) * (1 - exp( -a_j * (t - u) )) / a_i / a_j;
		};

		return integrator_( integrand, s, t );
	}

	Real GPPConstantMeanReversion::variance( Size i, Size j, Time s, Time t ) const
	{
		Real a_i = a( i, 0.0 );
		Real a_j = a( j, 0.0 );

		const Parameter& sigma_i = sigma(i);
		const Parameter& sigma_j = sigma(j);

		auto integrand = [&, a_i, a_j, sigma_i, sigma_j, t]( Time u )
		{
			return sigma_i( u ) * sigma_j( u ) * exp(  - ( a_i + a_j ) * (t - u) );
		};

		return integrator_( integrand, s, t );
	}
}
