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
		for ( Size i = 0; i < rho.rows(); i++ )
		{
			for ( Size j = i; j < rho.columns(); j++ )
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

	Real GaussianFactorDynamics::E( Size i, Time t, Time T ) const
	{
		Parameter a = a_[i];
		auto lambda = [&a]( Time u )
		{
			return a( u );
		};

		Real val = integrator_( lambda, t, T );
		return exp( val );
	}

	Real GaussianFactorDynamics::B( Size i, Time t, Time T ) const
	{
		auto lambda = [&, t, T]( Time u )
		{
			return 1 / E( i, t, u );
		};

		return integrator_( lambda, t, T );
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
			return sigma_i( u ) * sigma_j( u ) * B( j, u, T ) / E(i, u, t);
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
			return sigma_i( u ) * sigma_j( u ) * B( i, j, u, t );
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
			return sigma_i( u ) * sigma_j( u ) / E( i, j, u, t );
		};

		return integrator_( integrand, s, t );
	}

	Real GaussianFactorDynamics::phi( Size i, Size j, Time t ) const
	{
		const Parameter& a_i = a_[i];
		const Parameter& a_j = a_[j];

		const Parameter& sigma_i = sigma_[i];
		const Parameter& sigma_j = sigma_[j];

		auto integrand = [&, t]( Time u )
		{
			return sigma_i( u ) * sigma_j( u )
				* (B( j, u, t ) / E( i, u, t ) + B( i, u, t ) / E( j, u, t ));
		};

		return integrator_( integrand, 0, t );
	}
}