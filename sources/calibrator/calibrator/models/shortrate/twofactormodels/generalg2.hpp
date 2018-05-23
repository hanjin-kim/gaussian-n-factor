#ifndef CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP
#define CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP

#include <boost/make_shared.hpp>

#include <ql/models/shortrate/twofactormodel.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/math/solver1d.hpp>
#include <ql/math/integrals/kronrodintegral.hpp>

#include <calibrator/global.hpp>
#include <calibrator/processes/generalornsteinuhlenbeckprocess.hpp>
#include <calibrator/models/shortrate/dynamics/gaussianfactordynamics.hpp>

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

	\todo Tree implementation

	\ingroup shortrate
	*/
	class GeneralizedG2 : public TwoFactorModel, public AffineModel, public TermStructureConsistentModel
	{
		class Dynamics;

		shared_ptr<Gaussian2FactorDynamics> dynamics_;

	public:
		GeneralizedG2( shared_ptr<Gaussian2FactorDynamics> dynamics,
					   Real integralSignificance = 10, 
					   shared_ptr<Integrator> integrator = boost::make_shared<GaussKronrodAdaptive>( GaussKronrodAdaptive( 1.e-8, 10000 ) ) );
		virtual ~GeneralizedG2() {}

		// TwoFactorModel virtual override
		shared_ptr<Lattice> tree( const TimeGrid& grid ) const override
		{
			// todo
			return shared_ptr<Lattice>();
		};

		shared_ptr<ShortRateDynamics> dynamics() const;

		virtual DiscountFactor discount( Time t ) const override
		{
			return termStructure()->discount( t );
		}

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

		virtual Real swaption( const Swaption::arguments& arg, Real strike ) const;

		Parameter a() const { return a_; }
		Parameter b() const { return b_; }
		Parameter sigma() const { return sigma_; }
		Parameter eta() const { return eta_; }
		Parameter rho() const { return rho_; }

	protected:
		// CalibratedModel virtual override
		virtual void generateArguments() override;

		Real A( Time t, Time T ) const;

		Parameter& a_;
		Parameter& sigma_;
		Parameter& b_;
		Parameter& eta_;
		Parameter& rho_;

		Real integralSignificance_;

		shared_ptr<Integrator> integrator_;
	};

	//! Short-rate dynamics in the time-dependent Hull-White model
	/*! The short-rate follows an time-dependent Hull-White process */
	class GeneralizedG2::Dynamics : public TwoFactorModel::ShortRateDynamics {
	public:
		Dynamics( const shared_ptr<Gaussian2FactorDynamics> dynamics )
			: ShortRateDynamics( shared_ptr<StochasticProcess1D>( new GeneralizedOrnsteinUhlenbeckProcess( dynamics->a( 0 ), dynamics->sigma( 0 ) ) ),
								 shared_ptr<StochasticProcess1D>( new GeneralizedOrnsteinUhlenbeckProcess( dynamics->a( 1 ), dynamics->sigma( 1 ) ) ),
								 dynamics->rho( 0, 1 )(0.0) )
			, dynamics_( dynamics )
		{}

		virtual Real shortRate( Time t, Real x, Real y ) const {
			return x + y + dynamics_->phi( t );
		}

		shared_ptr<Gaussian2FactorDynamics> dynamics_;
	};

	// inline definitions
	inline shared_ptr<TwoFactorModel::ShortRateDynamics>	GeneralizedG2::dynamics() const
	{
		return shared_ptr<ShortRateDynamics>( new Dynamics( dynamics_ ) );
	}
}
#endif // !CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP