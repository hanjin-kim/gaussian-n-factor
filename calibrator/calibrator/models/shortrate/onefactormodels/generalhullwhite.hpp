#ifndef CALIBRATOR_PROCESSES_GENERAL_VASICEK_HPP
#define CALIBRATOR_PROCESSES_GENERAL_VASICEK_HPP

#include <boost/bind.hpp>

#include <ql/models/shortrate/onefactormodel.hpp>

#include <calibrator/global.hpp>
#include <calibrator/processes/generalornsteinuhlenbeckprocess.hpp>

namespace HJCALIBRATOR
{
	//! Time-dependent %Hull-White model class
	/*! This class implements the time-dependent Hull-White model defined by
	\f[
	dr(t) = (\theta(t) - a(t)r(t))dt + \sigma(t) dW_t ,
	\f]
	a risk premium \f$ \lambda \f$ can also be specified.

	\ingroup shortrate
	*/
	class GeneralizedHullWhite : public OneFactorAffineModel, public TermStructureConsistentModel
	{
	public:
		GeneralizedHullWhite( const Handle<YieldTermStructure>& termStructure,
							  const IntegrableParameter& a,
							  const Parameter& sigma );

		// OnFactorModel virtual override
		boost::shared_ptr<Lattice> tree( const TimeGrid& grid ) const override 
		{
			return boost::shared_ptr<Lattice>();
		};
		virtual boost::shared_ptr<ShortRateDynamics> dynamics() const override;

		
		virtual Real discountBondOption( Option::Type type,
										 Real strike,
										 Time maturity,
										 Time bondMaturity ) const;

		Real discountBondOption( Option::Type type, Real strike,
							Time maturity, Time bondStart,
							Time bondMaturity ) const;

		Parameter a() const { return a_; }
		Parameter sigma() const { return sigma_; }

		Real a( Time t ) const { return a_( t ); }
		Real sigma( Time t ) const { return sigma_( t ); }

		/*! Futures convexity bias (i.e., the difference between
		futures implied rate and forward rate) calculated as in
		G. Kirikos, D. Novak, "Convexity Conundrums", Risk
		Magazine, March 1997.
		http://www.powerfinance.com/convexity/

		\note t and T should be expressed in yearfraction using
		deposit day counter, F_quoted is futures' market price.
		*/
		// to be implemented with variable a & \sigma
		/*
		static Rate convexityBias( Real futurePrice,
								   Time t,
								   Time T,
								   Real sigma,
								   Real a );
		*/

		static std::vector<bool> FixedReversion() {
			std::vector<bool> c( 2 );
			c[0] = true; c[1] = false;
			return c;
		}

	protected:
		void generateArguments();

		virtual Real A( Time t, Time T ) const;
		virtual Real B( Time t, Time T ) const;

		Parameter& a_;
		Parameter& sigma_;

	private:
		class Dynamics;
		class FittingParameter;

		boost::shared_ptr<FittingParameter> alpha_;
	};

	class GeneralizedHullWhite::FittingParameter : public TermStructureFittingParameter 
	{
	private:
		class Impl : public Parameter::Impl 
		{
		public:
			Impl( const Handle<YieldTermStructure>& termStructure,
				  const IntegrableParameter& a,
				  const Parameter& sigma );

			Real value( const Array&, Time t ) const override;

			Real B( Time t, Time T ) const;
			Real variance( Time s, Time t ) const;

		private:
			Real E( Time t0, Time t1, Real multiplier = 1. ) const;
			Real integrand( Time u, Time t ) const;
			Real integrandVr( Time t ) const;

			GaussKronrodAdaptive integrator_;
			boost::function<Real( Real )> OneOverEintegrand_; // Eq. 31 integrand
			boost::function<Real( Real )> Vrintegrand_; // Eq. 37 integrand

			Handle<YieldTermStructure> termStructure_;
			IntegrableParameter a_;
			Parameter sigma_;
		};
	public:
		FittingParameter( const Handle<YieldTermStructure>& termStructure,
						  const IntegrableParameter& a,
						  const Parameter& sigma )
			: TermStructureFittingParameter( boost::shared_ptr<Parameter::Impl>(
				new FittingParameter::Impl( termStructure, a, sigma ) ) ) {}

		Real B( Time t, Time T ) const
		{
			return boost::static_pointer_cast<Impl>(implementation())->B( t, T );
		}


		Real variance( Time t, Time T ) const
		{
			return boost::static_pointer_cast<Impl>(implementation())->variance( t, T );
		}
	};

	//! Short-rate dynamics in the time-dependent Hull-White model
	/*! The short-rate follows an time-dependent Hull-White process */
	class GeneralizedHullWhite::Dynamics : public OneFactorModel::ShortRateDynamics {
	public:
		Dynamics( const FittingParameter& fitting,
				  const IntegrableParameter& a,
				  const Parameter& sigma )
			: ShortRateDynamics( boost::shared_ptr<StochasticProcess1D>(
				new GeneralizedOrnsteinUhlenbeckProcess( a, sigma ) ) )
			, fitting_( fitting )
		{}

		virtual Real variable( Time t, Rate r ) const {
			return r - fitting_( t );
		}
		virtual Real shortRate( Time t, Real x ) const {
			return x + fitting_( t );
		}

		FittingParameter fitting_;
	};

	// inline definitions
	inline boost::shared_ptr<OneFactorModel::ShortRateDynamics>
		GeneralizedHullWhite::dynamics() const 
	{
		return boost::shared_ptr<ShortRateDynamics>(
			new Dynamics( *alpha_, static_cast<IntegrableParameter>(a()), sigma() ) );
	}
	inline Real GeneralizedHullWhite::FittingParameter::Impl::value( const Array&, Time t ) const
	{
		Rate forwardRate = termStructure_->forwardRate( t, t, Continuous, NoFrequency );
		Real Et = E( 0, t );

		boost::function<Real( Real )> I_t;
		I_t = boost::bind( &GeneralizedHullWhite::FittingParameter::Impl::integrand, this, _1, t );

		return forwardRate + integrator_( I_t, 0, t ) / Et;
	}

	inline Real GeneralizedHullWhite::FittingParameter::Impl::integrand( Time u, Time t ) const
	{
		/* eq.36 */
		Real sigma_u = sigma_( u );

		return E( 0, u ) * sigma_u * sigma_u * B( u, t );
	}

	inline Real GeneralizedHullWhite::FittingParameter::Impl::integrandVr( Time t ) const
	{
		Real sigma_t = sigma_( t );
		Real E_t = E( 0, t );

		return E_t * E_t * sigma_t * sigma_t;
	}


	inline Real GeneralizedHullWhite::FittingParameter::Impl::E( Time t0, Time t1, Real multiplier ) const
	{
		/* eq. 30 */
		return exp( multiplier * a_.integral( t0, t1 ) );
	}

	inline Real GeneralizedHullWhite::FittingParameter::Impl::B( Time t, Time T ) const
	{
		/* eq. 31 */
		return E( 0, t ) * integrator_( OneOverEintegrand_, t, T );
	}

	inline Real GeneralizedHullWhite::FittingParameter::Impl::variance( Time s, Time t ) const
	{
		/* eq. 37 */
		Real E_t = E( 0, t );

		return integrator_( Vrintegrand_, s, t ) / E_t / E_t;
	}
}

#endif // !CALIBRATOR_PROCESSES_GENERAL_VASICEK_HPP
