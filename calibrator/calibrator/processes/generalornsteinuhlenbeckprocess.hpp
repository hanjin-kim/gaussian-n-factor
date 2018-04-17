#ifndef CALIBRATOR_PROCESSES_GENERALORNSTEINUHLENBECKPROCESS_HPP
#define CALIBRATOR_PROCESSES_GENERALORNSTEINUHLENBECKPROCESS_HPP

#include <boost/function.hpp>

#include <gsl/gsl_integration.h>

#include <ql/stochasticprocess.hpp>
#include <ql/math/integrals/simpsonintegral.hpp>
#include <ql/math/integrals/kronrodintegral.hpp>

#include <calibrator/global.hpp>
#include <calibrator/models/integrableparameter.hpp>
#include <calibrator/models/shortrate/dynamics/gaussianfactordynamics.hpp>

namespace HJCALIBRATOR
{
	// Use in combination with boost::bind.
	template<class F>
	static double gslFunctionAdapter( double x, void* p )
	{
		// Here I do recover the "right" pointer, safer to use static_cast
		// than reinterpret_cast.
		F* function = static_cast<F*>(p);
		return (*function)(x);
	}

	template<class F>
	gsl_function convertToGslFunction( const F& f )
	{
		gsl_function gslFunction;

		const void* p = &f;
		assert( p != 0 );

		gslFunction.function = &gslFunctionAdapter<F>;
		// Just to eliminate the const.
		gslFunction.params = const_cast<void*>(p);

		return gslFunction;
	}

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
		GeneralizedOrnsteinUhlenbeckProcess( const IntegrableParameter& a, // must be an IntegrableParameter
											 const Parameter& sigma,
											 const Real x0 = 0 );

		GeneralizedOrnsteinUhlenbeckProcess( const Parameter& a, // must be an IntegrableParameter
											 const Parameter& sigma,
											 const Real x0 = 0 )
		{}

		virtual ~GeneralizedOrnsteinUhlenbeckProcess();

		//! \name StochasticProcess1D interface
		//@{
		Real x0() const override { return x0_; }
		Real drift( Time t, Real x ) const override;
		Real diffusion( Time t, Real x ) const override;
		Real variance( Time t0, Real x0, Time dt ) const override;
		Real expectation( Time t0, Real x0, Time dt ) const override;
		Real stdDeviation( Time t0, Real x0, Time dt ) const override;
		//@}

		IntegrableParameter a() { return a_; }
		Parameter sigma() { return sigma_; }

		Real a( Time t ) { return a_( t ); }
		Real sigma( Time t ) { return sigma_( t ); }

	private :
		Real E( Time t0, Time t1 ) const;
		Real VrIntegrand( Time t );

		Real x0_;
		IntegrableParameter a_;
		Parameter sigma_;

		gsl_integration_workspace* int_wrkspcs_;
		gsl_function integrand_Vr_;

		//boost::function<Real( Real )> Vrintegrand_; // variance integrand
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
		Real Et = E( 0, t0 + dt );
		
		Real result, error;
		gsl_integration_qags( &integrand_Vr_, t0, t0 + dt, 0, 1e-7, 1000, int_wrkspcs_, &result, &error );
		
		return result / (Et * Et);
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
		return exp( a_.integral( t0, t1 ) );
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
