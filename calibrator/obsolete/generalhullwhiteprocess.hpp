#ifndef CALIBRATOR_PROCESSES_GENERAL_HULLWHITEPROCESS_HPP
#define CALIBRATOR_PROCESSES_GENERAL_HULLWHITEPROCESS_HPP

#include <boost/bind.hpp>

#include <ql/stochasticprocess.hpp>
#include <ql/processes/forwardmeasureprocess.hpp>
#include <ql/math/integrals/simpsonintegral.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>

#include <calibrator/global.hpp>
#include <calibrator/models/integrableparameter.hpp>

namespace HJCALIBRATOR
{
	//! Time-dependent Hull-White process class
	/*! This class describes the time-dependent Hull-White process governed by
	\f[
	dr(t) = (\theta(t) - a(t)r(t)) dt + \sigma(t) dW_t.
	\f]

	Derivation Reference, Eq. numbering quoted from SSRN-id1514192

	\ingroup processes
	*/
	class GeneralizedHullWhiteProcess : public StochasticProcess1D
	{
	public :
		GeneralizedHullWhiteProcess( const Handle<YieldTermStructure>& h,
									 const IntegrableParameter& a, // must be an IntegrableParameter
									 const Parameter& sigma );
		
		//! \name StochasticProcess1D interface
		//@{
		Real x0() const override { return r0_; }
		Real drift( Time t, Rate r ) const override;
		Real diffusion( Time t, Rate r ) const override;
		Real variance( Time t0, Rate r0, Time dt ) const override;
		Real expectation( Time t0, Rate r0, Time dt ) const override;
		Real stdDeviation( Time t0, Rate r0, Time dt ) const override;
		//@}

		IntegrableParameter a() { return a_; }
		Parameter sigma() { return sigma_; }

		Real a( Time t ) { return a_( t ); }
		Real sigma( Time t ) { return sigma_( t ); }

		Real bondPrice( Time t, Time T, Rate r )const;

	private :
		Real A( Time t, Time T ) const;
		Real B( Time t, Time T ) const;

		Real theta( Time t ) const;
		Real E( Time t, Real multiplier = 1. ) const;
		Real alpha( Time t ) const;

		Real DVIntegrand( Time u, Time t ) const;
		Real VrIntegrand( Time t ) const;
		Real alphaIntegrand( Time u, Time t ) const;

		Handle<YieldTermStructure> termStructure_;
		IntegrableParameter a_;
		Parameter sigma_;
		Real r0_;

		SimpsonIntegral integrator_;
		boost::function<Real( Real )> Vrintegrand_; // variance integrand
		boost::function<Real( Real )> OneOverEintegrand_; // Eq. 31 integrand

	};

	//! Time-dependent Hull-White Forward-Rate process class
	/*
	Derivation Reference, Eq. numbering quoted from SSRN-id1514192

	\ingroup processes
	*/
	/*
	class GeneralizedHullwhiteForwardProcess : public ForwardMeasureProcess1D
	{
		GeneralizedHullwhiteForwardProcess( const Handle<YieldTermStructure>& h,
											const IntegrableParameter& a, // must be an IntegrableParameter
											const Parameter& sigma );

		//! \name StochasticProcess1D interface
		//@{
		Real x0() const override;
		Real drift( Time t, Rate r ) const override;
		Real diffusion( Time t, Rate r ) const override;
		Real variance( Time t0, Rate r0, Time dt ) const override;
		Real expectation( Time t0, Rate r0, Time dt ) const override;
		Real stdDeviation( Time t0, Rate r0, Time dt ) const override;
		//@}

	private :
		boost::shared_ptr<GeneralizedHullWhiteProcess> hwprocess_;
	};
	*/

	// inline definitions

	inline Real GeneralizedHullWhiteProcess::drift( Time t, Rate r ) const
	{
		return theta( t ) - a_( t ) * r;
	}

	inline Real GeneralizedHullWhiteProcess::diffusion( Time t, Rate r ) const
	{
		return sigma_( t );
	}

	inline Real GeneralizedHullWhiteProcess::variance( Time t0, Rate r0, Time dt ) const
	{
		Real Et = E( t0 + dt );

		return integrator_( Vrintegrand_, t0, t0 + dt ) / (Et * Et);
	}

	inline Real GeneralizedHullWhiteProcess::stdDeviation( Time t0, Rate r0, Time dt ) const
	{
		return sqrt( variance( t0, r0, dt ) );
	}

	inline Real GeneralizedHullWhiteProcess::expectation( Time t0, Rate r0, Time dt ) const
	{
		/* A part of eq.35 */
		Time s = t0;
		Time t = t0 + dt;

		Real RE = E( s ) / E( t );

		return RE * r0 + alpha( t ) - RE * alpha( s );
	}

	Real GeneralizedHullWhiteProcess::bondPrice( Time t, Time T, Rate r) const
	{
		return A( t, T ) * exp( -(B( t, T ) * r) );
	}

	Real GeneralizedHullWhiteProcess::A( Time t, Time T ) const
	{
		/* Modefication of the eq. 43 so that we have the general affine form P = Aexp(-Br)
		\f[
		A(t,T) = \frac{P(0,T)}{P(0,t)} + exp\left{B(t,T)f(0,t) - \frac{1}{2}B^2(t,T)V_r(0,t)\right}
		\f]
		*/
		Real discount_t = termStructure_->discount( t );
		Real discount_T = termStructure_->discount( T );
		Real forward = termStructure_->forwardRate( t, t, Continuous, NoFrequency );

		Real BtT = B( t, T );
		Real Vrt = variance( 0, 0, t );

		return (discount_T / discount_t) * exp( BtT*forward - 0.5 * BtT * BtT * Vrt );
	}

	Real GeneralizedHullWhiteProcess::B( Time t, Time T ) const
	{
		/* eq. 31
		\f[
		B(t,T) = E(t) \int_t^T \frac{du}{E(u)}
		\f]
		*/
		return E( t ) * integrator_( OneOverEintegrand_, t, T );
	}

	inline Real GeneralizedHullWhiteProcess::theta( Time t ) const
	{
		const Real dt = 0.000001; // Should it be variable?
		Real f = termStructure_->forwardRate( t, t, Continuous, NoFrequency );
		Real fup = termStructure_->forwardRate( t + dt, t + dt, Continuous, NoFrequency );
		Real f_prime = (fup - f) / dt;

		boost::function<Real( Real )> integrand;
		integrand = boost::bind( &GeneralizedHullWhiteProcess::DVIntegrand, this, _1, t );
		Real IntI = integrator_( integrand, 0, t );

		Real Et = E( t );
		Real at = a_( t );

		Real DV = (2. / Et) * IntI;
		Real D2V = 2.*variance( 0, 0, t ) - (2*at / Et) * IntI;

		/* eq. 39 */
		return f_prime + a_( t ) * f + 0.5 * (D2V + at * DV);
	}

	inline Real GeneralizedHullWhiteProcess::E( Time t, Real multiplier ) const
	{
		/* eq. 30
		\f[
		E(t) = exp(\int_0^t{a(u)du})
		\f]

		
		The multiplier is given to deal with the integrand of the eq. 31 (\f$ B(t,T) \f$)
		*/
		return exp( multiplier * a_.integral( 0, t ) );
	}

	inline Real GeneralizedHullWhiteProcess::DVIntegrand( Time u, Time t ) const
	{
		/* Integrand of \f$ \partial V(0,t) /\partial t \f$ in its derived form
		\f[
		\partial V(0,t) /\partial t = \frac{2}{E(t)\int_0^t \sigma(u,t)\sigma(u)E(u)du
		\f]

		Here, we takes \f$ I(u,t) = \sigma(u,t)\sigma(u)E(u) \f$
		so we have 
		\f[
		\partial V(0,t) /\partial t = \frac{2}{E(t)\int_0^t I(u,t)du
		\partial^2 V(0,t) /\partial t^2 = -\frac{a(t)}{E(t)}\int_0^t I(u,t)du + 2*V_r(0,t)
		\f]

		where \f$ V_r(0,t \f$ is the variance of the short rate.
		\f]
		*/

		Real sigma_u = sigma_( u );

		return (sigma_u * B( u, t )) * sigma_u * E( u );
	}

	Real GeneralizedHullWhiteProcess::VrIntegrand( Time t ) const
	{
		/* Integrand of eq. 37
		\f[
		I(t) =  E^2(t)\sigam^2(t)
		\f]
		*/
		Real Et = E( t );
		Real sigmat = sigma_( t );

		return Et * Et * sigmat * sigmat;
	}

	Real GeneralizedHullWhiteProcess::alpha( Time t ) const
	{
		/* eq. 36 
		\f[
		\alpha = f(0,t) + \frac{1}{E(t)}\int_0^t E(u)\sigma(u)B(u,t)du
		*/
		boost::function<Real( Real )> integrand;
		integrand = boost::bind( &GeneralizedHullWhiteProcess::alphaIntegrand, this, _1, t );
		Real IntI = integrator_( integrand, 0, t );
		Real forward = termStructure_->forwardRate( t, t, Continuous, NoFrequency );

		return forward + IntI / E( t );
	}

	Real GeneralizedHullWhiteProcess::alphaIntegrand( Time u, Time t ) const
	{
		Real sigma_t = sigma_( t );
		
		return E( t ) * sigma_t * sigma_t * B( u, t );
	}
}

#endif // !CALIBRATOR_PROCESSES_GENERAL_HULLWHITEPROCESS_HPP

