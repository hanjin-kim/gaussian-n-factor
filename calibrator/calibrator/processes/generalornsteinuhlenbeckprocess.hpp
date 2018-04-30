#ifndef CALIBRATOR_PROCESSES_GENERALORNSTEINUHLENBECKPROCESS_HPP
#define CALIBRATOR_PROCESSES_GENERALORNSTEINUHLENBECKPROCESS_HPP

#include <boost/function.hpp>

#include <ql/stochasticprocess.hpp>
#include <ql/models/parameter.hpp>
#include <ql/math/integrals/kronrodintegral.hpp>

#include <calibrator/global.hpp>

namespace HJCALIBRATOR
{
	//! Time-dependelt Ornstein-Uhlenbeck process class
	/*! This class describes the (constrained) time-dependent Ornstein-Uhlenbeck process
	described by
	\f[
	dx(t) = -a(t)x(t)dt + \sigma(t)dW_t
	\f]

	\note The level(\f$ \mu \f$) is constrained to be zero,
	and the speed term \f$ \alpha \f$ must be given as a IntegrableParameter object.
	
	\ingroup processes
	*/
	class GeneralizedOrnsteinUhlenbeckProcess : public StochasticProcess1D
	{
	public :
		GeneralizedOrnsteinUhlenbeckProcess( const Parameter& a, // must be an IntegrableParameter
											 const Parameter& sigma,
											 const Real x0 = 0 )
			: a_(a), sigma_(sigma), x0_(x0)
			, integrator_(GaussKronrodAdaptive( GaussKronrodAdaptive( 1.e-8, 10000 ) ))
		{}

		virtual ~GeneralizedOrnsteinUhlenbeckProcess() {}

		//! \name StochasticProcess1D interface
		//@{
		Real x0() const override { return x0_; }
		Real drift( Time t, Real x ) const override;
		Real diffusion( Time t, Real x ) const override;
		Real variance( Time t0, Real x0, Time dt ) const override;
		Real expectation( Time t0, Real x0, Time dt ) const override;
		Real stdDeviation( Time t0, Real x0, Time dt ) const override;
		//@}

		Parameter a() { return a_; }
		Parameter sigma() { return sigma_; }

		Real a( Time t ) { return a_( t ); }
		Real sigma( Time t ) { return sigma_( t ); }

	private :
		Real E( Time t0, Time t1 ) const;
		Real VrIntegrand( Time t );

		Real x0_;
		Parameter a_;
		Parameter sigma_;

		GaussKronrodAdaptive integrator_;
	};

	// inline definitions

	inline Real GeneralizedOrnsteinUhlenbeckProcess::drift( Time t, Real x ) const
	{
		return  - a_( t ) * x;
	}

	inline Real GeneralizedOrnsteinUhlenbeckProcess::diffusion( Time t, Real x ) const
	{
		return sigma_( t );
	}

	inline Real GeneralizedOrnsteinUhlenbeckProcess::variance( Time t0, Real x0, Time dt ) const
	{
		Time t = t0 + dt;

		auto integrand = [&, t]( Time u )
		{
			Real sigma = sigma_( u );
			Real Etu = E( u, t );

			return sigma * sigma / Etu / Etu;
		};
		
		return integrator_( integrand, t0, t );
	}

	inline Real GeneralizedOrnsteinUhlenbeckProcess::stdDeviation( Time t0, Real x0, Time dt ) const
	{
		return sqrt( variance( t0, x0, dt ) );
	}

	inline Real GeneralizedOrnsteinUhlenbeckProcess::expectation( Time t0, Real x0, Time dt ) const
	{
		/* A part of eq.35 */
		Time s = t0;
		Time t = t0 + dt;

		Real RE = 1 / E( s, t );

		return RE * x0;
	}

	inline Real GeneralizedOrnsteinUhlenbeckProcess::E( Time t0, Time t1 ) const
	{
		const Parameter& a = a_;

		auto integrand = [&a]( Time u )
		{
			return a( u );
		};

		return integrator_( integrand, t0, t1 );
	}

	inline Real GeneralizedOrnsteinUhlenbeckProcess::VrIntegrand( Time t )
	{
		/* Integrand of eq. 37
		\f[
		I(t) =  E^2(t)\sigam^2(t)
		\f]
		*/
		Real Et = E( 0, t );
		Real sigmat = sigma_( t );

		return Et * Et * sigmat * sigmat;
	}
}
#endif // !CALIBRATOR_PROCESSES_GENERALORNSTEINUHLENBECKPROCESS_HPP
