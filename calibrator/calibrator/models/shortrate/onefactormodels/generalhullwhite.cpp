#include <ql/pricingengines/blackformula.hpp>

#include <calibrator/models/shortrate/onefactormodels/generalhullwhite.hpp>

namespace HJCALIBRATOR
{
	GeneralizedHullWhite::GeneralizedHullWhite( const Handle<YieldTermStructure>& termStructure,
												const IntegrableParameter& a,
												const Parameter& sigma )
		: OneFactorAffineModel( 2 )
		, TermStructureConsistentModel( termStructure )
		, a_( arguments_[0] ), sigma_( arguments_[1] )
	{
		a_ = a;
		sigma_ = sigma;

		generateArguments();

		registerWith( termStructure );
	}

	void GeneralizedHullWhite::generateArguments() 
	{
		alpha_ = boost::shared_ptr<FittingParameter>( new FittingParameter( termStructure(), a(), sigma() ) );
	}

	Real GeneralizedHullWhite::A( Time t, Time T ) const 
	{
		auto termStructure_ = termStructure();
		Real discount_t = termStructure_->discount( t );
		Real discount_T = termStructure_->discount( T );
		Real forward = termStructure_->forwardRate( t, t, Continuous, NoFrequency );

		Real BtT = alpha_->B( t, T );
		Real Vrt = alpha_->variance( 0, t );

		return (discount_T / discount_t) * exp( BtT*forward - 0.5 * BtT * BtT * Vrt );

	}

	Real GeneralizedHullWhite::B( Time t, Time T ) const 
	{
		return alpha_->B( t, T );
	}

	GeneralizedHullWhite::FittingParameter::Impl::Impl( const Handle<YieldTermStructure>& termStructure,
														const IntegrableParameter& a,
														const Parameter& sigma )
		: termStructure_( termStructure )
		, a_( a )
		, sigma_( sigma )
		, integrator_( GaussKronrodAdaptive( 0.00001, 1000 ) )
		, OneOverEintegrand_( boost::bind( &GeneralizedHullWhite::FittingParameter::Impl::E, this, 0, _1, -1. ) )
		, Vrintegrand_( boost::bind( &GeneralizedHullWhite::FittingParameter::Impl::integrandVr, this, _1 ) )
	{}

	Real GeneralizedHullWhite::discountBondOption( Option::Type type,
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
		Real Vr_t_TF = alpha_->variance( 0, maturity );
		Real B_TF_TP = alpha_->B( maturity, bondMaturity );
		Real Vp_0_TF_TP = Vr_t_TF * B_TF_TP * B_TF_TP;

		/* eq. 11 */
		Real stdDev = sqrt( std::max( Vp_0_TF_TP, 0.0 ) );

		Real f = termStructure()->discount( bondMaturity );
		Real k = termStructure()->discount( maturity )*strike;

		return blackFormula( type, k, f, stdDev );
	}

	Real GeneralizedHullWhite::discountBondOption( Option::Type type,
												   Real strike,
												   Time maturity, // corresponds to now?
												   Time bondStart,
												   Time bondMaturity ) const
	{
		return discountBondOption( type, strike, bondStart, bondMaturity );
	}

	/*
	   Rate GeneralizedHullWhite::convexityBias(Real futuresPrice,
                                  Time t,
                                  Time T,
                                  Real sigma,
                                  Real a) {
        QL_REQUIRE(futuresPrice>=0.0,
            "negative futures price (" << futuresPrice << ") not allowed");
        QL_REQUIRE(t>=0.0,
            "negative t (" << t << ") not allowed");
        QL_REQUIRE(T>=t,
            "T (" << T << ") must not be less than t (" << t << ")");
        QL_REQUIRE(sigma>=0.0,
            "negative sigma (" << sigma << ") not allowed");
        QL_REQUIRE(a>=0.0,
            "negative a (" << a << ") not allowed");

        Time deltaT = (T-t);
        Real tempDeltaT = (1.-std::exp(-a*deltaT)) / a;
        Real halfSigmaSquare = sigma*sigma/2.0;

        // lambda adjusts for the fact that the underlying is an interest rate
        Real lambda = halfSigmaSquare * (1.-std::exp(-2.0*a*t)) / a *
            tempDeltaT * tempDeltaT;

        Real tempT = (1.0 - std::exp(-a*t)) / a;

        // phi is the MtM adjustment
        Real phi = halfSigmaSquare * tempDeltaT * tempT * tempT;

        // the adjustment
        Real z = lambda + phi;

        Rate futureRate = (100.0-futuresPrice)/100.0;
        return (1.0-std::exp(-z)) * (futureRate + 1.0/(T-t));
    }
	*/
}