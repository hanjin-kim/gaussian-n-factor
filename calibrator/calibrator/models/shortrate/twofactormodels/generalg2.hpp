#ifndef CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP
#define CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP

#include <ql/models/shortrate/twofactormodel.hpp>

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
	class GeneralizedG2 : public TwoFactorModel, 
		public AffineModel, public TermStructureConsistentModel
	{
	public :
		GeneralizedG2( const Handle<YieldTermStructure>& termStructure,
			const IntegrableParameter& a,
			const Parameter& sigma,
			const IntegrableParameter& b,
			const Parameter& eta,
			const Real rho );

		boost::shared_ptr<ShortRateDynamics> dynamics() const;

		Real discountBond( Time, Time, Rate, Rate ) const;

		Real discountBondOption( Option::Type type,
								 Real strike,
								 Time maturity,
								 Time bondMaturity ) const;

		/*
		Real swaption( const Swaption::arguments& arguments,
					   Rate fixedRate,
					   Real range,
					   Size intervals ) const;
		*/

		Parameter a() const { return a_); }
		Parameter sigma() const { return sigma_; }
		Parameter b() const { return b_; }
		Parameter eta() const { return eta_; }
		Real rho() const { return rho_(0.0); }

		Real a( Time t ) const { return a_( t ); }
		Real sigma( Time t ) const { return sigma_( t ); }
		Real b( Time t ) const { return b_( t ); }
		Real eta( Time t ) const { return eta_( t ); }

	protected :
		void generateArguments();

		virtual Real A( Time t, Time T ) const;
		virtual Real B( Time t, Time T ) const;

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
			IntegrableParameter b_;
			Parameter eta_;
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
}

#endif // !CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP
