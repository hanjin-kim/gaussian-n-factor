#ifndef CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCONSTANTDYNAMICS_HPP
#define CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCONSTANTDYNAMICS_HPP

#include <calibrator/models/shortrate/dynamics/gaussianfactor/g++constantmeanreversion.hpp>

namespace HJCALIBRATOR
{
	class GPPConstantDynamics : public GPPConstantMeanReversion
	{
	public :
		GPPConstantDynamics( const Handle<YieldTermStructure>& termStructure,
							 const RealVector& a,
							 const RealVector& sigma,
							 const Matrix& rho )
			: GPPConstantMeanReversion( termStructure, 
										a, RealVectorToParamVector( sigma, PositiveConstraint() ),
										rho )
		{}

		virtual ~GPPConstantDynamics() {}

		virtual Real meanTforward( Size i, Size j, Time T, Time s, Time t ) const override;
		virtual Real integralVariance( Size i, Size j, Time s, Time t ) const override;
		virtual Real variance( Size i, Size j, Time s, Time t ) const override;
		//virtual Real phi( Size i, Size j, Time t ) const override;

	protected:
		GPPConstantDynamics() {}
	};

	class G1ConstantDynamics : public Gaussian1FactorDynamics, public GPPConstantDynamics
	{
	public:
		G1ConstantDynamics( const Handle<YieldTermStructure>& termStructure,
							Real a,
							Real sigma )
			: GaussianFactorDynamics( termStructure,
									  { ConstantParameter( a, NoConstraint() ) },
									  { ConstantParameter( sigma , PositiveConstraint() ) },
									  Matrix( 1, 1, 1 ) )
		{}
		virtual ~G1ConstantDynamics() {}
	};

	class G2ConstantDynamics : public Gaussian2FactorDynamics, public GPPConstantDynamics
	{
	public:
		G2ConstantDynamics( const Handle<YieldTermStructure>& termStructure,
							Real a, Real sigma,
							Real b, Real eta,
							Real rho )
			: GaussianFactorDynamics( termStructure,
									  { ConstantParameter( a, PositiveConstraint() ),  ConstantParameter( b, PositiveConstraint() ) },
									  { ConstantParameter( sigma, PositiveConstraint() ), ConstantParameter( eta, PositiveConstraint() ) },
									  getCorrelationMatrix( rho ) )
		{}
		virtual ~G2ConstantDynamics() {}
	};
}

#endif // !CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCONSTANTDYNAMICS_HPP
