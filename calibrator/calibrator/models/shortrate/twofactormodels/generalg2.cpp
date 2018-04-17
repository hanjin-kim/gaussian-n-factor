#include <ql/pricingengines/blackformula.hpp>

#include <calibrator/models/shortrate/twofactormodels/generalg2.hpp>

namespace HJCALIBRATOR
{
	GeneralizedG2::GeneralizedG2( std::shared_ptr<Gaussian2FactorDynamics> dynamics )
		: TwoFactorModel( 5 )
		, AffineModel()
		, TermStructureConsistentModel( dynamics->termStructure() )
		, a_( arguments_[0] ), sigma_( arguments_[1] )
		, b_( arguments_[2] ), eta_( arguments_[3] )
		, rho_( arguments_[4] )
		, integrator_( GaussKronrodAdaptive( 0.01, 10000 ) )
	{
		a_ = dynamics->a(0);
		b_ = dynamics->a(1);
		sigma_ = dynamics->sigma(0);
		eta_ = dynamics->sigma(1);
		rho_ = dynamics->rho( 0, 1 );

		generateArguments();

		registerWith( dynamics->termStructure() );
	}

	void GeneralizedG2::generateArguments() 
	{
		dynamics_->a( a_, 0 );
		dynamics_->sigma( sigma_, 0 );
		dynamics_->a( b_, 1 );
		dynamics_->sigma( eta_, 1 );
		dynamics_->rho( rho_, 1, 0 );
	}

	Real GeneralizedG2::A( Time t, Time T ) const
	{
		return dynamics_->A( t, T );
	}

	Real GeneralizedG2::discountBond( Time now,
									  Time maturity,
									  Array factors ) const
	{
		Real x = factors[0];
		Real y = factors[1];

		Real Bx = dynamics_->B( 0, now, maturity );
		Real By = dynamics_->B( 1, now, maturity );

		return A( now, maturity ) * exp( - Bx * x - By * y );
	}


	Real GeneralizedG2::discountBondOption( Option::Type type,
											Real strike,
											Time maturity,
											Time bondMaturity ) const
	{
		auto dynamics = dynamics_;

		Real Bx = dynamics->B( 0, maturity, bondMaturity );
		Real By = dynamics->B( 1, maturity, bondMaturity );
		Real rho = dynamics->rho( 0, 1 )(0.0);

		Real Vpratio = 0;

		auto integrand00 = [maturity, dynamics]( Time u )
		{
			Real sigma = dynamics->sigma( 0 )(u);
			Real E = dynamics->E( 0, u, maturity );
			
			return sigma * sigma / E / E;
		};

		Vpratio += Bx * By * integrator_( integrand00, 0, maturity );

		auto integrand11 = [maturity, dynamics]( Time u )
		{
			Real sigma = dynamics->sigma( 1 )(u);
			Real E = dynamics->E( 1, u, maturity );

			return sigma * sigma / E / E;
		};

		Vpratio += Bx * By * integrator_( integrand11, 0, maturity );

		auto integrand01 = [maturity, dynamics]( Time u )
		{
			return dynamics->sigma( 0 )(u) * dynamics->sigma( 1 )(u) 
				/ dynamics->E( 0, u, maturity ) / dynamics->E( 1, u, maturity );
		};

		Vpratio += 2 * Bx * Bx * rho * integrator_( integrand01, 0, maturity );

		Real stdDev = sqrt( std::max( Vpratio, 0.0 ) );

		Real f = termStructure()->discount( bondMaturity );
		Real k = termStructure()->discount( maturity )*strike;

		return blackFormula( type, k, f, stdDev );
	}


	Real GeneralizedG2::discountBondOption( Option::Type type,
											Real strike,
											Time maturity,
											Time bondStart,
											Time bondMaturity ) const
	{
		return discountBondOption( type, strike, bondStart, bondMaturity );
	}
}