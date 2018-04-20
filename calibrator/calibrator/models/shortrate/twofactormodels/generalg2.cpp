#include <boost/bind.hpp>

#include <ql/time/date.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/pricingengines/blackformula.hpp>
#include <ql/math/solvers1d/brent.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/math/integrals/segmentintegral.hpp>

#include <calibrator/models/shortrate/twofactormodels/generalg2.hpp>

namespace HJCALIBRATOR
{
	GeneralizedG2::GeneralizedG2( boost::shared_ptr<Gaussian2FactorDynamics> dynamics,
								  Real integralSignificance,
								  boost::shared_ptr<Integrator> integrator )
		: TwoFactorModel( 5 )
		, AffineModel()
		, TermStructureConsistentModel( dynamics->termStructure() )
		, a_( arguments_[0] ), sigma_( arguments_[1] )
		, b_( arguments_[2] ), eta_( arguments_[3] )
		, rho_( arguments_[4] )
		, integralSignificance_( integralSignificance )
		, integrator_( integrator )
		, dynamics_( dynamics )
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

		Vpratio += Bx * By * integrator_->operator()( integrand00, 0, maturity );

		auto integrand11 = [maturity, dynamics]( Time u )
		{
			Real sigma = dynamics->sigma( 1 )(u);
			Real E = dynamics->E( 1, u, maturity );

			return sigma * sigma / E / E;
		};

		Vpratio += Bx * By * integrator_->operator()( integrand11, 0, maturity );

		auto integrand01 = [maturity, dynamics]( Time u )
		{
			return dynamics->sigma( 0 )(u) * dynamics->sigma( 1 )(u) 
				/ dynamics->E( 0, u, maturity ) / dynamics->E( 1, u, maturity );
		};

		Vpratio += 2 * Bx * Bx * rho * integrator_->operator()( integrand01, 0, maturity );

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

	// Brigo Ch. 4.2
	Real GeneralizedG2::swaption( const Swaption::arguments& arg, Real strike ) const
	{
		Date settlement = termStructure()->referenceDate();
		DayCounter dayCounter = termStructure()->dayCounter();
		Time T = dayCounter.yearFraction( settlement,
											  arg.floatingResetDates[0] );
		Real w = (arg.type == VanillaSwap::Payer ? 1 : -1);

		std::vector<Time> t;
		for ( auto fixedPayDate : arg.fixedPayDates )
		{
			t.push_back( dayCounter.yearFraction( settlement,
													fixedPayDate ) );
		}
		Size N_timestep = t.size();
		
		Array cA( N_timestep );
		Array Bx( N_timestep );
		Array By( N_timestep );

		for ( Size i = 0; i < N_timestep; i++ )
		{
			Time tau_i = i == 0 ? t[i] - T : t[i] - t[i - 1];
			Real c = i == N_timestep - 1 ? 1 + strike * tau_i : strike * tau_i;
			cA[i] = c * dynamics_->A( T, t[i] );
			Bx[i] = dynamics_->B( 0, T, t[i] );
			By[i] = dynamics_->B( 1, T, t[i]);
		}

		Real mu_x = dynamics_->meanTforward( 0, T, 0, T );
		Real mu_y = dynamics_->meanTforward( 1, T, 0, T );
		Real sigma_x = sqrt(dynamics_->variance( 0, 0, 0, T ));
		Real sigma_y = sqrt(dynamics_->variance( 1, 1, 0, T ));
		Real rho = dynamics_->rho( 0, 1 )(0.0);
		Real var = dynamics_->variance( 0, 1, 0, T );
		Real rho_xy = rho * var / sigma_x / sigma_y;
		Real rhosqrt = sqrt( 1 - rho_xy * rho_xy );

		auto integrand = [&, N_timestep, w, mu_x, mu_y, sigma_x, sigma_y, rho_xy, rhosqrt]( Real x )
		{
			Real dev = (x - mu_x) / sigma_x;

			Array lambda( N_timestep );
			Array kappa( N_timestep );
			for ( Size i = 0; i < N_timestep; i++ )
			{
				lambda[i] = cA[i] * exp( -Bx[i] * x );
				kappa[i] = -By[i] * (mu_y - 0.5 * rhosqrt * rhosqrt * sigma_y * sigma_y * By[i]
											+ rho_xy * sigma_y * (x - mu_x) / sigma_x);
			}

			auto hyperplane = [N_timestep, &lambda, &By]( Real y )
			{
				Real value = 1.;
				for ( Size i = 0; i < N_timestep; i++ )
				{
					Real val = lambda[i] * exp( -By[i] * y );
					value -= val;
				}

				return value;
			};

			Brent solver;
			solver.setMaxEvaluations( 1000 );
			Real ybar = solver.solve( hyperplane, 1e-6, 0.00, -100.0, 100.0 );

			Real h1 = (ybar - mu_y) / (sigma_y * rhosqrt)
				- rho_xy * (x - mu_x) / (sigma_x * rhosqrt);

			CumulativeNormalDistribution Phi;
			Real val = Phi( -w * h1 );
			for ( Size i = 0; i < N_timestep; i++ )
			{
				Real h2 = h1 + By[i] * sigma_y * rhosqrt;

				val -= lambda[i] * exp( kappa[i] ) * Phi( -w * h2 );
			}

			return exp( -0.5*dev*dev ) * val;
		};

		Real N = arg.nominal;
		Real P0T = termStructure()->discount( T );
		Real upper = mu_x + integralSignificance_ *sigma_x;
		Real lower = mu_x - integralSignificance_ *sigma_x;

		Real val = N * w * P0T * integrator_->operator()( integrand, lower, upper ) / sqrt( 2. * M_PI ) / sigma_x;
		return val;
	}
}