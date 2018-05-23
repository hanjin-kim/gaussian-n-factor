#ifndef CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCONSTANTMEANREVERSION_HPP
#define CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCONSTANTMEANREVERSION_HPP

#include <calibrator/models/shortrate/dynamics/gaussianfactordynamics.hpp>

namespace HJCALIBRATOR
{
	ParamVector&& RealVectorToParamVector( const RealVector& a, Constraint constraint = NoConstraint() );

	class GPPConstantMeanReversion : public virtual GaussianFactorDynamics
	{
	public:
		GPPConstantMeanReversion( const Handle<YieldTermStructure>& termStructure,
								  const RealVector& a,
								  const ParamVector& sigma,
								  const Matrix& rho )
			: GaussianFactorDynamics( termStructure, RealVectorToParamVector( a ), sigma, rho )
		{}

		virtual ~GPPConstantMeanReversion() {}


		virtual Real E( Size i, Time s, Time t ) const;
		virtual Real B( Size i, Time s, Time t ) const;

		virtual Real integralVariance( Size i, Size j, Time s, Time t ) const override;
		virtual Real variance( Size i, Size j, Time s, Time t ) const override;

	protected:
		GPPConstantMeanReversion() {}

		virtual Real phi( Size i, Size j, Time t ) const;
	};

	class G1ConstantMeanReversionDynamics : public Gaussian1FactorDynamics, public GPPConstantMeanReversion
	{
	public:
		G1ConstantMeanReversionDynamics( const Handle<YieldTermStructure>& termStructure,
										 Real a,
										 const Parameter& sigma )
			: GaussianFactorDynamics( termStructure,
									  { ConstantParameter( a, NoConstraint() ) },
									  { sigma },
									  Matrix( 1, 1, 1 ) )
		{}

		virtual ~G1ConstantMeanReversionDynamics() {}
	};
}

#endif // !CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCONSTANTMEANREVERSION_HPP
