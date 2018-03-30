#ifndef CALIBRATOR_MODEL_SHORTRATE_GENERAL_HULLWHITE1F_HPP
#define CALIBRATOR_MODEL_SHORTRATE_GENERAL_HULLWHITE1F_HPP

#include <ql/models/shortrate/onefactormodels/vasicek.hpp>

#include <calibrator/global.hpp>
#include <calibrator/models/shortrate/onefactormodels/general_vasicek.hpp>

namespace HJCALIBRATOR
{
	//! Single-factor, time-dependent Hull-White (extended %Vasicek) model class.
	/*! This class implements the time-dependent single-factor Hull-White model
	defined by
	\f[
	dr_t = (\theta(t) - \alpha r_t)dt + \sigma_t dW_t
	\f]
	where \f$ \alpha \f$ and \f$ \sigma \f$ are constants.

	\test calibration results are tested against cached values

	\bug When the term structure is relinked, the r0 parameter of
	the underlying Vasicek model is not updated.

	\ingroup shortrate
	*/
	class GeneralizedHullWhite : public GeneralizedHullWhite, public TermStructureConsistentModel
	{
	public:
		GeneralizedHullWhite( const Handle<YieldTermStructure>& termStructure,
								const Parameter& alpha,
								const Parameter& sigma );

		boost::shared_ptr<Lattice> tree( const TimeGrid& grid ) const;

		boost::shared_ptr<ShortRateDynamics> dynamics() const;

		// Need to review implementations of the bond option price 
		Real discountBondOption( Option::Type type,
								 Real strike,
								 Time maturity,
								 Time bondMaturity ) const;

		Real discountBondOption( Option::Type type,
								 Real strike,
								 Time maturity,
								 Time bondStart,
								 Time bondMaturity ) const;

		/*! Futures convexity bias (i.e., the difference between
		futures implied rate and forward rate) calculated as in
		G. Kirikos, D. Novak, "Convexity Conundrums", Risk
		Magazine, March 1997.

		\note t and T should be expressed in yearfraction using
		deposit day counter, F_quoted is futures' market price.
		*/
		static Rate convexityBias( Real futurePrice,
								   Time t,
								   Time T,
								   Real sigma,
								   Real a );

		static std::vector<bool> FixedReversion() {
			std::vector<bool> c( 2 );
			c[0] = true; c[1] = false;
			return c;
		}

	protected:
		void generateArguments();

		Real A( Time t, Time T ) const;

	private:
		class Dynamics;
		class FittingParameter;

		Parameter phi_;
	};

	//! Short-rate dynamics in the time-dependent Hull-White model
	/*! The short-rate is here
	\f[
	r_t = \varphi(t) + x_t
	\f]
	where \f$ \varphi(t) \f$ is the deterministic time-dependent
	parameter used for term-structure fitting and \f$ x_t \f$ is the
	state variable following an Ornstein-Uhlenbeck process.
	*/
	class GeneralizedHullWhite::Dynamics : public OneFactorModel::ShortRateDynamics {
	public:
		Dynamics( const Parameter& fitting,
				  const Parameter& a,
				  const Parameter& sigma )
			: ShortRateDynamics( boost::shared_ptr<StochasticProcess1D>(
				new GeneralizedOrnsteinUhlenbeckProcess( a, sigma ) ) ),
			fitting_( fitting ) {}

		Real variable( Time t, Rate r ) const {
			return r - fitting_( t );
		}
		Real shortRate( Time t, Real x ) const {
			return x + fitting_( t );
		}
	private:
		Parameter fitting_;
	};

	//! Analytical term-structure fitting parameter \f$ \varphi(t) \f$.
	/*! \f$ \varphi(t) \f$ is analytically defined by
	\f[
	\varphi(t) = f(t) + \frac{1}{2}[\frac{\sigma(1-e^{-at})}{a}]^2,
	\f]
	where \f$ f(t) \f$ is the instantaneous forward rate at \f$ t \f$.
	*/
	class GeneralizedHullWhite::FittingParameter
		: public TermStructureFittingParameter {
	private:
		class Impl : public Parameter::Impl {
		public:
			Impl( const Handle<YieldTermStructure>& termStructure,
				  const Parameter& a, const Parameter& sigma )
				: termStructure_( termStructure ), a_( a ), sigma_( sigma ) {}

			Real value( const Array&, Time t ) const;
		private:
			Handle<YieldTermStructure> termStructure_;
			Parameter a_;
			Parameter sigma_;
		};
	public:
		FittingParameter( const Handle<YieldTermStructure>& termStructure,
						  const Parameter& a, const Parameter& sigma )
			: TermStructureFittingParameter( boost::shared_ptr<Parameter::Impl>(
				new FittingParameter::Impl( termStructure, a, sigma ) ) ) {}
	};

	// inline definitions

	inline boost::shared_ptr<OneFactorModel::ShortRateDynamics>
		GeneralizedHullWhite::dynamics() const {
		return boost::shared_ptr<ShortRateDynamics>(
			new Dynamics( phi_, a(), sigma() ) );
	}
}

#endif // !CALIBRATOR_MODEL_SHORTRATE_GENERALIZED_HULLWHITE1F_HPP
