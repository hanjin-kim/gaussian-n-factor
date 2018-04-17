#ifndef CALIBRATOR_MODELS_SHORTRATE_GNPP_HPP
#define CALIBRATOR_MODELS_SHORTRATE_GNPP_HPP

#include <ql/models/shortrate/onefactormodel.hpp>

#include <calibrator/global.hpp>
#include <calibrator/models/shortrate/dynamics/gaussianfactordynamics.hpp>
#include <calibrator/processes/generalornsteinuhlenbeckprocess.hpp>

namespace HJCALIBRATOR
{
	class GeneralizedG1 : public OneFactorAffineModel, public TermStructureConsistentModel
	{
		class Dynamics;

		std::shared_ptr<Gaussian1FactorDynamics> dynamics_;
		
	public :
		GeneralizedG1( std::shared_ptr<Gaussian1FactorDynamics> dynamics );
		virtual ~GeneralizedG1() {}

		// OneFactorModel virtual override
		boost::shared_ptr<Lattice> tree( const TimeGrid& grid ) const override
		{
			// todo
			return boost::shared_ptr<Lattice>();
		};

		virtual boost::shared_ptr<ShortRateDynamics> dynamics() const override;


		virtual Real discountBondOption( Option::Type type,
										 Real strike,
										 Time maturity,
										 Time bondMaturity ) const override;

		virtual Real discountBondOption( Option::Type type, Real strike,
										 Time maturity, Time bondStart,
										 Time bondMaturity ) const override;

	private :
		// CalibratedModel virtual override
		virtual void generateArguments() override;

		// OneFactorAffineModel virtual override
		virtual Real A( Time t, Time T ) const override;
		virtual Real B( Time t, Time T ) const override;

		Parameter& a_;
		Parameter& sigma_;
	};

	//! Short-rate dynamics in the time-dependent Hull-White model
	/*! The short-rate follows an time-dependent Hull-White process */
	class GeneralizedG1::Dynamics : public OneFactorModel::ShortRateDynamics {
	public:
		Dynamics( std::shared_ptr<Gaussian1FactorDynamics> dynamics )
			: ShortRateDynamics( boost::shared_ptr<StochasticProcess1D>(
				new GeneralizedOrnsteinUhlenbeckProcess( dynamics->a( 0 ), dynamics->sigma( 0 ) ) ) )
			, dynamics_( dynamics )
		{}

		virtual Real variable( Time t, Rate r ) const {
			return r - dynamics_->phi( t );
		}
		virtual Real shortRate( Time t, Real x ) const {
			return x + dynamics_->phi( t );
		}

		std::shared_ptr<Gaussian1FactorDynamics> dynamics_;
	};

	// inline definitions
	inline boost::shared_ptr<OneFactorModel::ShortRateDynamics>	GeneralizedG1::dynamics() const
	{
		return boost::shared_ptr<ShortRateDynamics>( new Dynamics( dynamics_ ) );
	}
}

#endif // !CALIBRATOR_MODELS_SHORTRATE_GNPP_HPP
