#ifndef CALIBRATOR_MODELS_PARAMETERS_SHORTRATE_GAUSSIANFACTORFITTINGPARAMETER_HPP
#define CALIBRATOR_MODELS_PARAMETERS_SHORTRATE_GAUSSIANFACTORFITTINGPARAMETER_HPP

#include <gsl/gsl_integration.h>

#include <ql/math/matrix.hpp>
#include <ql/math/integrals/gaussianorthogonalpolynomial.hpp>
#include <ql/math/randomnumbers/sobolrsg.hpp>
#include <ql/math/integrals/kronrodintegral.hpp>


#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/models/parameter.hpp>

#include <calibrator/global.hpp>

namespace HJCALIBRATOR
{
	using ParamVector = std::vector<Parameter>;
	using ParamMatrix = std::vector<std::vector<Parameter>>;

	class GaussianFactorDynamics
	{
	public:
		GaussianFactorDynamics( const Handle<YieldTermStructure>& termStructure,
								const ParamVector& a,
								const ParamVector& sigma,
								const Matrix& rho );
	public:
		Size dof() { return a_.size(); }

		virtual Real E( Size i, Time s, Time t ) const { return E( a_[i], s, t ); }
		virtual Real B( Size i, Time s, Time t ) const { return B( a_[i], s, t ); }

		virtual Real phi( Time t ) const;
		virtual Real variance( Time s, Time t ) const;
		virtual Real integralVariance( Time s, Time t ) const;

	private:
		Real E( const Parameter& a, Time s, Time t ) const;
		Real B( const Parameter& a, Time s, Time t ) const;

		Real integralVariance( Size i, Size j, Time s, Time t ) const;
		Real variance( Size i, Size j, Time s, Time t ) const;
		Real phi( Size i, Size j, Time t ) const;

		Handle<YieldTermStructure> termStructure_;
		ParamVector a_;
		ParamVector sigma_;
		ParamMatrix rho_;

		GaussKronrodAdaptive integrator_;
	};

	Real GaussianFactorDynamics::phi( Time t ) const
	{
		Rate forwardRate = termStructure_->forwardRate( t, t, Continuous, NoFrequency );
		
		Real intsum = 0;
		
		for ( Size i = 0; i < a_.size(); i++ )
		{
			for ( Size j = i; j < a_.size(); j++ )
			{
				Real integ = phi( i, j, t );

				if(j == i) intsum += integ;
				else intsum += 0.5 * rho_[i][j]( 0.0 ) * integ;
			}
		}

		return forwardRate + intsum;
	}

	Real GaussianFactorDynamics::variance( Time s, Time t ) const
	{
		Real intsum = 0;

		for ( Size i = 0; i < a_.size(); i++ )
		{
			for ( Size j = i; j < a_.size(); j++ )
			{
				Real integ = variance( i, j, s, t );

				if ( j == i ) intsum += 2*integ;
				else intsum += rho_[i][j]( 0.0 ) * integ;
			}
		}

		return intsum;
	}

	Real GaussianFactorDynamics::integralVariance( Time s, Time t ) const
	{
		Real intsum = 0;

		for ( Size i = 0; i < a_.size(); i++ )
		{
			for ( Size j = i; j < a_.size(); j++ )
			{
				Real integ = integralVariance( i, j, s, t );

				if ( j == i ) intsum += 2 * integ;
				else intsum += rho_[i][j]( 0.0 ) * integ;
			}
		}

		return intsum;
	}

	Real GaussianFactorDynamics::E( const Parameter& a, Time t, Time T ) const
	{
		auto lambda = [&, a]( Time u )
		{
			return a( u );
		};

		return integrator_( lambda, t, T );
	}

	Real GaussianFactorDynamics::B( const Parameter& a, Time t, Time T ) const
	{
		auto lambda = [&, a, T]( Time u )
		{
			return 1 / E( a, t, u );
		};

		return integrator_( lambda, t, T );
	}

	Real GaussianFactorDynamics::integralVariance( Size i, Size j, Time s, Time t ) const
	{
		// todo
		return 0;
	}

	Real GaussianFactorDynamics::variance( Size i, Size j, Time s, Time t ) const
	{
		const Parameter& a_i = a_[i];
		const Parameter& a_j = a_[j];

		const Parameter& sigma_i = sigma_[i];
		const Parameter& sigma_j = sigma_[j];

		auto integrand = [&, a_i, a_j, sigma_i, sigma_j, t]( Time u )
		{
			return sigma_i( u ) * sigma_j( u ) / E( a_i, u, t ) / E( a_j, u, t );
		};

		return integrator_( integrand, s, t );
	}

	Real GaussianFactorDynamics::phi( Size i, Size j, Time t ) const
	{
		const Parameter& a_i = a_[i];
		const Parameter& a_j = a_[j];

		const Parameter& sigma_i = sigma_[i];
		const Parameter& sigma_j = sigma_[j];

		auto integrand = [&, a_i, a_j, sigma_i, sigma_j, t]( Time u, Time s )
		{	
			return sigma_i( u ) * sigma_j( u )
				* (1 / E( a_i, u, s ) + 1 / E( a_j, u, s ))
				/ E( a_i, u, t ) / E( a_i, u, t );
		};

		// monte carlo integration
		SobolRsg sobol( 3 );

		Size N = 0; Size n = 0;
		for ( Size i = 0; i < 10000; i++ )
		{
			auto seq = sobol.nextSequence();

			// value[0] : u - [0,t], value[1] : s - [u,t]
			Real u = t * seq.value[0];
			Real s = t * seq.value[1];
			Real val = t * seq.value[2]; // problematic..

			if ( u < s )
			{
				continue;
			}
			else
			{
				N++;
				if ( integrand( u, s ) < val ) n++;
			}
		}

		return (double( n ) / double( N )) * t * t;
	}
}

#endif // !CALIBRATOR_MODELS_PARAMETERS_SHORTRATE_GAUSSIANFACTORFITTINGPARAMETER_HPP
