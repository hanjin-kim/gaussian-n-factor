#include <ql/pricingengines/blackformula.hpp>

#include <calibrator/models/shortrate/onefactormodels/generalg1.hpp>

namespace HJCALIBRATOR
{
	GeneralizedG1::GeneralizedG1( std::shared_ptr<Gaussian1FactorDynamics> dynamics )
		: OneFactorAffineModel( 2 )
		, TermStructureConsistentModel( dynamics->termStructure() )
		, a_( arguments_[0] ), sigma_( arguments_[1] )
		, dynamics_( dynamics )
	{
		QL_ENSURE( dynamics->dimension() == 1,
				   "The dimension of the dynamics exceeds one" );

		a_ = dynamics->a( 0 );
		sigma_ = dynamics->sigma( 0 );

		registerWith( dynamics->termStructure() );
	}

	void GeneralizedG1::generateArguments()
	{
		dynamics_->a( a_, 0 );
		dynamics_->sigma( sigma_, 0 );
	}


	Real GeneralizedG1::A( Time t, Time T ) const
	{
		return dynamics_->A( t, T );
	}

	Real GeneralizedG1::B( Time t, Time T ) const
	{
		return dynamics_->B( 0, t, T );
	}

	Real GeneralizedG1::discountBondOption( Option::Type type,
												   Real strike,
												   Time maturity,
												   Time bondMaturity ) const
	{
		/*
		From eq. 10 - 11
		Zero-bond put option price corresponds to the black formula :
		\f[
		ZBP(T_F,T_P,X) = XP(0, T_F)\mathcal{N}(d_+) - P(0, T_P)\mathcal{N}(d_-)
		\f]

		where \f$ T_F, T_P, X \f$ correspond to the fixing time (bondStart),
		the paying time (bondMaturity), and the strike, respectively.

		The variance is given as \f$ V_P(0, T_F, T_P) \f$ where :
		\f[
		V_P(t, T_F, T_P) = V_r(t, T_F)B^2(T_F,T_P)
		\f]
		which is the eq. 8.

		\note Unlike the original implemantation, as noted in SSRN-2246054, maturity is not necessary.
		*/

		/* eq. 8 */
		Real Vr_t_TF = dynamics_->variance( 0, maturity );
		Real B_TF_TP = dynamics_->B( 0, maturity, bondMaturity );
		Real Vp_0_TF_TP = Vr_t_TF * B_TF_TP * B_TF_TP;

		/* eq. 11 */
		Real stdDev = sqrt( std::max( Vp_0_TF_TP, 0.0 ) );

		Real f = termStructure()->discount( bondMaturity );
		Real k = termStructure()->discount( maturity )*strike;

		return blackFormula( type, k, f, stdDev );
	}

	Real GeneralizedG1::discountBondOption( Option::Type type,
												   Real strike,
												   Time maturity,
												   Time bondStart,
												   Time bondMaturity ) const
	{
		return discountBondOption( type, strike, bondStart, bondMaturity );
	}

}