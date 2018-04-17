#ifndef CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP
#define CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP

#include <ql/models/shortrate/twofactormodel.hpp>

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

	\ingroup shortrate
	*/
	class GeneralizedG2 : public TwoFactorModel, public AffineModel, public TermStructureConsistentModel
	{
		class Dynamics;

		boost::shared_ptr<Gaussian2FactorDynamics> dynamics_;

	public:
		GeneralizedG2( std::shared_ptr<Gaussian2FactorDynamics> dynamics );
		virtual ~GeneralizedG2() {}

		// TwoFactorModel virtual override
		boost::shared_ptr<Lattice> tree( const TimeGrid& grid ) const override
		{
			// todo
			return boost::shared_ptr<Lattice>();
		};

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

		// pure virtual for now (need to obtain a general form)
		virtual Real swaption( Option::Type type, Real strike,
							   Time maturity, Time swapTenor ) const = 0; 

	protected:
		// CalibratedModel virtual override
		virtual void generateArguments() override;

		Real A( Time t, Time T ) const;

		Parameter& a_;
		Parameter& sigma_;
		Parameter& b_;
		Parameter& eta_;
		Parameter& rho_;

		GaussKronrodAdaptive integrator_;
	};

	//! Short-rate dynamics in the time-dependent Hull-White model
	/*! The short-rate follows an time-dependent Hull-White process */
	class GeneralizedG2::Dynamics : public TwoFactorModel::ShortRateDynamics {
	public:
		Dynamics( const boost::shared_ptr<Gaussian2FactorDynamics> dynamics )
			: ShortRateDynamics( boost::shared_ptr<StochasticProcess1D>( new GeneralizedOrnsteinUhlenbeckProcess( dynamics->a( 0 ), dynamics->sigma( 0 ) ) ),
								 boost::shared_ptr<StochasticProcess1D>( new GeneralizedOrnsteinUhlenbeckProcess( dynamics->a( 1 ), dynamics->sigma( 1 ) ) ),
								 dynamics->rho( 0, 1 )(0.0) )
			, dynamics_( dynamics )
		{}

		virtual Real shortRate( Time t, Real x, Real y ) const {
			return x + y + dynamics_->phi( t );
		}

		boost::shared_ptr<Gaussian2FactorDynamics> dynamics_;
	};

	// inline definitions
	inline boost::shared_ptr<TwoFactorModel::ShortRateDynamics>	GeneralizedG2::dynamics() const
	{
		return boost::shared_ptr<ShortRateDynamics>( new Dynamics( dynamics_ ) );
	}
}
#endif // !CALIBRATOR_MODELS_SHORTRATE_TWOFACTORMODELS_GENERALG2_HPP