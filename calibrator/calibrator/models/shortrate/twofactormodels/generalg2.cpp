#include <calibrator/models/shortrate/twofactormodels/generalg2.hpp>

namespace HJCALIBRATOR
{
	GeneralizedG2::GeneralizedG2( const Handle<YieldTermStructure>& termStructure,
								  const IntegrableParameter& a,
								  const Parameter& sigma,
								  const IntegrableParameter& b,
								  const Parameter& eta,
								  const Real rho )
		: TwoFactorModel( 5 )
		, AffineModel()
		, TermStructureConsistentModel( termStructure )
		, a_( arguments_[0] ), sigma_( arguments_[1] )
		, b_( arguments_[2] ), eta_( arguments_[3] )
		, rho_( arguments_[4] )
	{
		a_ = a;
		b_ = b;
		sigma_ = sigma;
		eta_ = eta;
		rho_ = ConstantParameter( rho, BoundaryConstraint( -1.0, 1.0 ) );

		generateArguments();

		registerWith( termStructure );
	}

	void GeneralizedG2::generateArguments() {

		phi_ = boost::shared_ptr<FittingParameter>( new FittingParameter( termStructure(),
																		  a(), sigma(), b(), eta(), rho() ) );
	}

	Real GeneralizedG2::A( Time t, Time T ) const
	{
		auto termStructure_ = termStructure();
		Real discount_t = termStructure_->discount( t );
		Real discount_T = termStructure_->discount( T );

		Real VtT = phi_->V_I( t, T );
		Real V0T = phi_->V_I( 0, T );
		Real V0t = phi_->V_I( 0, t );

		return (discount_T / discount_t) * exp( -0.5 * (V0T - VtT - V0t) );
	}

	Real GeneralizedG2::V_I( Time t, Time T ) const
	{
		
	}


	Real GeneralizedG2::B( Time t, Time T ) const
	{
		return phi_->B0( t, T );
	}

	Real GeneralizedG2::B1( Time t, Time T ) const
	{
		return phi_->B1( t, T );
	}

	Real GeneralizedG2::discountBond( Time now,
									  Time maturity,
									  Array factors ) const
	{
		Real x = factors[0];
		Real y = factors[1];

		return A( now, maturity ) * exp( -B0( now, maturity ) * x - B1( now, maturity ) * y );
	}
}