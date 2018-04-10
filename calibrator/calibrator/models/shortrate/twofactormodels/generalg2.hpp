#ifndef CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP
#define CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP

#include <boost/bind.hpp>

#include <ql/models/shortrate/twofactormodel.hpp>
#include <ql/instruments/swaption.hpp>

#include <calibrator/global.hpp>
#include <calibrator/processes/generalornsteinuhlenbeckprocess.hpp>
#include <calibrator/models/integrableparameter.hpp>

namespace HJCALIBRATOR
{
	//! Two-additive-factor gaussian model class
	/*! This class implements a two-additive-factor model defined by
	\f[
	dr_t = \varphi(t) + x_t + y_t
	\f]
	where \f$ x_t \f$ and \f$ y_t \f$ are defined by
	\f[
	dx_t = -a x_t dt + \sigma dW^1_t, x_0 = 0
	\f]
	\f[
	dy_t = -b y_t dt + \sigma dW^2_t, y_0 = 0
	\f]
	and \f$ dW^1_t dW^2_t = \rho dt \f$.

	\bug This class was not tested enough to guarantee
	its functionality.

	\ingroup shortrate
	*/
	class GeneralizedG2 : public TwoFactorModel, public AffineModel, public TermStructureConsistentModel
	{
	public:
		GeneralizedG2( const Handle<YieldTermStructure>& termStructure,
					   const IntegrableParameter& a,
					   const Parameter& sigma,
					   const IntegrableParameter& b,
					   const Parameter& eta,
					   const Real rho );

		boost::shared_ptr<ShortRateDynamics> dynamics() const;

		virtual Real discountBond( Time now,
								   Time maturity,
								   Array factors ) const;

		virtual Real discountBondOption( Option::Type type,
										 Real strike,
										 Time maturity,
										 Time bondMaturity ) const override;

		virtual Real discountBondOption( Option::Type type, Real strike,
										 Time maturity, Time bondStart,
										 Time bondMaturity ) const override;

		Parameter a() const { return a_; }
		Parameter sigma() const { return sigma_; }
		Parameter b() const { return b_; }
		Parameter eta() const { return eta_; }
		Real rho() const { return rho_( 0.0 ); }

		Real a( Time t ) const { return a_( t ); }
		Real sigma( Time t ) const { return sigma_( t ); }
		Real b( Time t ) const { return b_( t ); }
		Real eta( Time t ) const { return eta_( t ); }

	protected:
		void generateArguments();

		Real A( Time t, Time T ) const;

		Parameter& a_;
		Parameter& sigma_;
		Parameter& b_;
		Parameter& eta_;
		Parameter& rho_;

	private:
		class Dynamics;
		class FittingParameter;

		boost::shared_ptr<FittingParameter> phi_;
	};


	class GeneralizedG2::FittingParameter : public TermStructureFittingParameter
	{
	private:
		class Impl : public Parameter::Impl
		{
		public:
			Impl( const Handle<YieldTermStructure>& termStructure,
				  const IntegrableParameter& a,
				  const Parameter& sigma,
				  const IntegrableParameter& b,
				  const Parameter& eta,
				  Real rho );

			Real value( const Array&, Time t ) const override;

			Real B( Size i, Time t, Time T ) const;

			Real variance( Time s, Time t ) const;
			Real integrandV_I( Time T, Time u ) const;

		private:
			Real E( Size i, Time t0, Time t1, Real multiplier = 1. ) const;

			// value integrands
			Real integrand_VI1F( Size i, Time u, Time t ) const;
			Real integrand_VI2F( Time u, Time t ) const;

			// variance integrands
			Real integrandVr00( Time t ) const;
			Real integrandVr01( Time t ) const;
			Real integrandVr11( Time t ) const;

			GaussKronrodAdaptive integrator_;
			boost::function<Real( Real )> fnc_E_;
			boost::function<Real( Real )> fnc_OneOverE_;

			boost::function<Real( Real )> OneOverE0integrand_; // Eq. 31 integrand
			boost::function<Real( Real )> OneOverE1integrand_; // Eq. 31 integrand

			boost::function<Real( Real )> Vr00integrand_;
			boost::function<Real( Real )> Vr01integrand_;
			boost::function<Real( Real )> Vr11integrand_;

			Handle<YieldTermStructure> termStructure_;
			IntegrableParameter a_[2];
			Parameter sigma_[2];
			Real rho_;
		};
	public:
		FittingParameter( const Handle<YieldTermStructure>& termStructure,
						  const IntegrableParameter& a,
						  const Parameter& sigma,
						  const IntegrableParameter& b,
						  const Parameter& eta,
						  Real rho )
			: TermStructureFittingParameter( boost::shared_ptr<Parameter::Impl>(
				new FittingParameter::Impl( termStructure, a, sigma, b, eta, rho ) ) ) {}

		Real B0( Time t, Time T ) const
		{
			return boost::static_pointer_cast<Impl>(implementation())->B0( t, T );
		}

		Real B1( Time t, Time T ) const
		{
			return boost::static_pointer_cast<Impl>(implementation())->B1( t, T );
		}

		Real variance( Time t, Time T ) const
		{
			return boost::static_pointer_cast<Impl>(implementation())->variance( t, T );
		}

		Real V_I( Time t, Time T ) const
		{
			return boost::static_pointer_cast<Impl>(implementation())->V_I( t, T );
		}
	};

	//! Short-rate dynamics in the time-dependent Hull-White model
	/*! The short-rate follows an time-dependent Hull-White process */
	class GeneralizedG2::Dynamics : public TwoFactorModel::ShortRateDynamics {
	public:
		Dynamics( const FittingParameter& fitting,
				  const IntegrableParameter& a,
				  const Parameter& sigma,
				  const IntegrableParameter& b,
				  const Parameter& eta,
				  Real rho )
			: ShortRateDynamics( boost::shared_ptr<StochasticProcess1D>( new GeneralizedOrnsteinUhlenbeckProcess( a, sigma ) ),
								 boost::shared_ptr<StochasticProcess1D>( new GeneralizedOrnsteinUhlenbeckProcess( b, eta ) ),
								 rho )
			, fitting_( fitting )
		{}

		virtual Real shortRate( Time t, Real x, Real y ) const {
			return x + y + fitting_( t );
		}

		FittingParameter fitting_;
	};

	Real GeneralizedG2::FittingParameter::Impl::E( Size i, Time t0, Time t1, Real multiplier ) const
	{
		return exp( multiplier * a_[i].integral( t0, t1 ) );
	}

	// inline definitions
	inline boost::shared_ptr<TwoFactorModel::ShortRateDynamics>
		GeneralizedG2::dynamics() const
	{
		return boost::shared_ptr<ShortRateDynamics>(
			new Dynamics( *phi_, static_cast<IntegrableParameter>(a()), sigma(),
						  static_cast<IntegrableParameter>(b()), eta(), rho() ) );
	}

	Real GeneralizedG2::FittingParameter::Impl::value( const Array&, Time t ) const
	{
		Rate forwardRate = termStructure_->forwardRate( t, t, Continuous, NoFrequency );

		Real a_t = a_( t );
		Real b_t = b_( t );
		Real sigma_t = sigma_( t );
		Real eta_t = eta_( t );

		Real E0_t = E( a_, 0, t );
		Real E1_t = E( b_, 0, t );

		boost::function<Real( Real )> I00_t, I01_t, I11_t;
		I00_t = boost::bind( &GeneralizedG2::FittingParameter::Impl::integrand00, this, _1, t );
		I01_t = boost::bind( &GeneralizedG2::FittingParameter::Impl::integrand01, this, _1, t );
		I11_t = boost::bind( &GeneralizedG2::FittingParameter::Impl::integrand11, this, _1, t );

		Real I00 = integrator_( I00_t, 0, t );
		Real I01 = integrator_( I01_t, 0, t );
		Real I11 = integrator_( I11_t, 0, t );

		return forwardRate
			+ 0.5 * (sigma_t * sigma_t - 2 * a_t * I00)
			+ rho_ * (sigma_t * eta_t - (a_t + b_t) * I01)
			+ 0.5 * (eta_t * eta_t - 2 * b_t * I00);
	}

	inline Real GeneralizedG2::FittingParameter::Impl::integrand_VI1F( Size i, Time u, Time t ) const
	{
		Real sigma_u = sigma_[i]( u );
		Real E_u = E( 0, u, t );
		Real B_ut

		return sigma_u * sigma_u / E0_u / E0_u;
	}
	
	inline Real GeneralizedG2::FittingParameter::Impl::integrand_VI2F( Time u, Time t ) const
	{
		return sigma_( u ) * eta_( u ) / E( a_, u, t ) / E( b_, u, t );
	}
	
	Real GeneralizedG2::FittingParameter::Impl::variance( Time s, Time t ) const
	{
		Real E0_t = E( a_, 0, t );
		Real E1_t = E( b_, 0, t );

		Real I00 = integrator_( Vr00integrand_, s, t );
		Real I01 = integrator_( Vr01integrand_, s, t );
		Real I11 = integrator_( Vr11integrand_, s, t );

		return I00 / E0_t / E0_t
			+ 2. * rho_ * I01 / E0_t / E1_t
			+ I11 / E1_t / E1_t;
	}

	inline Real GeneralizedG2::FittingParameter::Impl::integrandV_I( Time T, Time u ) const
	{
		Real sigma = sigma_( u );
		Real eta = eta_( u );
		Real B0 = B( 0, u, T );
		Real B1 = B( 1, u, T );

		Real i_00 = sigma * sigma * B0 * B0;
		Real i_01 = sigma * eta * B0 * B1;
		Real i_11 = eta * eta * B1 * B1;

		return i_00 + 2 * i_01 + i_11;
	}

	inline Real  GeneralizedG2::FittingParameter::Impl::B( Size i, Time t, Time T ) const
	{
		return E( a_, 0, t ) * integrator_( OneOverE0integrand_, t, T );
	}

	Real  GeneralizedG2::FittingParameter::Impl::E( const IntegrableParameter& p, Time t0, Time t1, Real multiplier = 1. ) const
	{
		return exp( multiplier * p.integral( t0, t1 ) );
	}

	Real GeneralizedG2::FittingParameter::Impl::integrandVr00( Time t ) const
	{
		Real E_t = E( a_, 0, t );
		Real sigma_t = sigma_( t );
		
		return sigma_t * sigma_t * E_t * E_t;
	}

	Real GeneralizedG2::FittingParameter::Impl::integrandVr01( Time t ) const
	{
		return sigma_( t ) * sigma_( t ) * E( a_, 0, t ) * E( b_, 0, t );
	}

	Real GeneralizedG2::FittingParameter::Impl::integrandVr11( Time t ) const
	{
		Real E_t = E( b_, 0, t );
		Real eta_t = sigma_( t );

		return eta_t * eta_t * E_t * E_t;
	}
}

#endif // !CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP
