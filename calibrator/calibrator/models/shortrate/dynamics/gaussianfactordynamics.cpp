#include <calibrator/models/shortrate/dynamics/gaussianfactordynamics.hpp>

namespace HJCALIBRATOR
{
	GaussianFactorDynamics::GaussianFactorDynamics( const Handle<YieldTermStructure>& termStructure,
													const ParamVector& a,
													const ParamVector& sigma,
													const Matrix& rho )
		: termStructure_( termStructure )
		, a_( a ), sigma_( sigma )
		, integrator_( GaussKronrodAdaptive( 0.01, 10000 ) )
	{
		QL_REQUIRE( a.size() == sigma.size(),
					"The number of a and sigma does not coinciede." );

		QL_REQUIRE( a.size() == rho.columns(),
					"The dimension of the correlation matrix is not equal to the dimension of the mean reversion parameter" );

		QL_REQUIRE( rho.rows() == rho.columns(),
					"The correlation matrix provided is not square matrix." );

		setupCorrelMatrix( rho );
	}

	void GaussianFactorDynamics::setupCorrelMatrix( const Matrix& rho )
	{
		for ( int i = 0; i < rho.rows(); i++ )
		{
			for ( int j = i; j < rho.columns(); j++ )
			{
				rho_[std::minmax(i,j)] = ConstantParameter( rho[i][j], BoundaryConstraint( -1, 1 ) );
			}
		}
	}

	void GaussianFactorDynamics::a( const Parameter& a, Size i )
	{
		QL_ENSURE( i < a_.size()+1,
				   "Memory for the " << i << "-th mean reversion parameter is not allocated" );
		
		a_[i] = a;
	}

	void GaussianFactorDynamics::sigma( const Parameter& sigma, Size i )
	{
		QL_ENSURE( i < sigma_.size() + 1,
				   "Memory for the " << i << "-th volatility parameter is not allocated" );
		
		sigma_[i] = sigma;
	}

	void GaussianFactorDynamics::rho( const Parameter& rho, Size i, Size j )
	{
		rho_[std::minmax( i, j )] = rho;
	}

	Parameter GaussianFactorDynamics::rho( Size i, Size j ) const
	{
		auto it = rho_.find( std::minmax( i, j ) );

		QL_ENSURE( it != rho_.end(),
				   "(" << i << "," << " j)-th correlation factor not found." );

		return it->second;
	}

	Real GaussianFactorDynamics::A( Time t, Time T ) const
	{
		Real discount_t = termStructure_->discount( t );
		Real discount_T = termStructure_->discount( T );
		Real exponent = integralVariance( 0, T ) - integralVariance( t, T ) - integralVariance( 0, t );

		return discount_T / discount_t * exp( -0.5 * exponent );
	}

	Real GaussianFactorDynamics::phi( Time t ) const
	{
		Rate forwardRate = termStructure_->forwardRate( t, t, Continuous, NoFrequency );

		Real intsum = 0;

		for ( Size i = 0; i < a_.size(); i++ )
		{
			for ( Size j = 0; j <= i; j++ )
			{
				Real integ = phi( i, j, t );

				if ( j == i ) intsum += 0.5 * integ;
				else intsum += rho(i,j)(0.0) * integ;
			}
		}

		return forwardRate + intsum;
	}

	Real GaussianFactorDynamics::variance( Time s, Time t ) const
	{
		Real intsum = 0;

		for ( Size i = 0; i < a_.size(); i++ )
		{
			for ( Size j = 0; j <= i; j++ )
			{
				Real integ = variance( i, j, s, t );

				if ( j == i ) intsum += integ;
				else intsum += 2 * rho(i,j)(0.0) * integ;
			}
		}

		return intsum;
	}

	Real GaussianFactorDynamics::integralVariance( Time s, Time t ) const
	{
		Real intsum = 0;

		for ( Size i = 0; i < a_.size(); i++ )
		{
			for ( Size j = 0; j <= i; j++ )
			{
				Real integ = integralVariance( i, j, s, t );

				if ( j == i ) intsum += integ;
				else intsum += 2 * rho(i,j)(0.0) * integ;
			}
		}

		return intsum;
	}

	Real GaussianFactorDynamics::E( const Parameter& a, Time t, Time T ) const
	{
		auto lambda = [&a]( Time u )
		{
			return a( u );
		};

		Real val = integrator_( lambda, t, T );
		return exp( val );
	}

	Real GaussianFactorDynamics::B( const Parameter& a, Time t, Time T ) const
	{
		auto lambda = [&, t, T]( Time u )
		{
			return 1 / E( a, t, u );
		};

		return integrator_( lambda, t, T );
	}

	Real GaussianFactorDynamics::meanTforward( Size i, Time T, Time s, Time t ) const
	{
		Real intsum = 0;
		for ( Size j = 0; j < dimension(); j++ )
		{
			intsum += rho(i,j)(0.0) * meanTforward( i, j, T, s, t );
		}

		return - intsum;
	}

	Real GaussianFactorDynamics::meanTforward( Size i, Size j, Time T, Time s, Time t ) const
	{
		const Parameter& a_i = a_[i];
		const Parameter& a_j = a_[j];

		const Parameter& sigma_i = sigma_[i];
		const Parameter& sigma_j = sigma_[j];

		auto integrand = [&, T]( Time u )
		{
			return sigma_i( u ) * sigma_j( u ) * B( a_j, u, T ) / E(a_i, u, t);
		};

		return integrator_( integrand, s, t );
	}

	Real GaussianFactorDynamics::integralVariance( Size i, Size j, Time s, Time t ) const
	{
		const Parameter& a_i = a_[i];
		const Parameter& a_j = a_[j];

		const Parameter& sigma_i = sigma_[i];
		const Parameter& sigma_j = sigma_[j];

		auto integrand = [&, t]( Time u )
		{
			return sigma_i( u ) * sigma_j( u ) * B( a_i, a_j, u, t );
		};

		return integrator_( integrand, s, t );
	}

	Real GaussianFactorDynamics::variance( Size i, Size j, Time s, Time t ) const
	{
		const Parameter& a_i = a_[i];
		const Parameter& a_j = a_[j];

		const Parameter& sigma_i = sigma_[i];
		const Parameter& sigma_j = sigma_[j];

		auto integrand = [&, t]( Time u )
		{
			return sigma_i( u ) * sigma_j( u ) / E( a_i, a_j, u, t );
		};

		return integrator_( integrand, s, t );
	}

	// must be corrected I guess..
	Real GaussianFactorDynamics::phi( Size i, Size j, Time t ) const
	{
		const Parameter& a_i = a_[i];
		const Parameter& a_j = a_[j];

		const Parameter& sigma_i = sigma_[i];
		const Parameter& sigma_j = sigma_[j];

		auto integrand = [&, t]( Time u, Time s )
		{
			return sigma_i( u ) * sigma_j( u )
				* (1 / E( a_i, u, s ) + 1 / E( a_j, u, s ))
				/ E( a_i, a_j, u, t );
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