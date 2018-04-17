#ifndef CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCONSTANTMEANREVERSION_HPP
#define CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCONSTANTMEANREVERSION_HPP

#include <gsl/gsl_integration.h>

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

	protected:
		GPPConstantMeanReversion() : GaussianFactorDynamics() {}

	private:
		virtual Real E( const Parameter& a, Time s, Time t ) const override;
		virtual Real B( const Parameter& a, Time s, Time t ) const override;
		virtual Real integralVariance( Size i, Size j, Time s, Time t ) const override;
		virtual Real variance( Size i, Size j, Time s, Time t ) const override;
	};

	class G1ConstantMeanReversionDynamics : public GPPConstantMeanReversion, public Gaussian1FactorDynamics
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
