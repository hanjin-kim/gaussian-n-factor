#ifndef CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCMR_PCV_HPP
#define CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCMR_PCV_HPP

#include <calibrator/models/shortrate/dynamics/gaussianfactor/g++constantmeanreversion.hpp>

namespace HJCALIBRATOR
{
	ParamVector convertParamVector( const std::vector<RealVector>& sigma_node,
									const std::vector<RealVector>& initial_sigma );

	class GPPPCMRPCV : public GPPConstantMeanReversion
	{
		std::vector<std::vector<RealVector>> combined_nodes_;

	public:
		GPPPCMRPCV( const Handle<YieldTermStructure>& termStructure,
					const RealVector& a,
					const std::vector<RealVector>& sigma_nodes,
					const std::vector<RealVector>& initial_sigma,
					const Matrix& rho )
			: GPPConstantMeanReversion( termStructure, 
										a, 
										convertParamVector( sigma_nodes, initial_sigma ), rho )
		{
			combineNodes( sigma_nodes );
		}

		virtual ~GPPPCMRPCV() {}

	protected:
		GPPPCMRPCV( const std::vector<RealVector>& sigma_nodes )
		{
			combineNodes( sigma_nodes );
		}

	private:
		void combineNodes( const std::vector<RealVector>& sigma_nodes );


		virtual Real integralVariance( Size i, Size j, Time s, Time t ) const override;
		virtual Real variance( Size i, Size j, Time s, Time t ) const override;
		//virtual Real phi( Size i, Size j, Time t ) const override;
	};

	class G1PPPCMRPCV : public GPPPCMRPCV, public Gaussian1FactorDynamics
	{
	public:
		G1PPPCMRPCV( const Handle<YieldTermStructure>& termStructure,
					 Real a,
					 const RealVector& sigma_node,
					 const RealVector& initial_sigma )
			: GaussianFactorDynamics( termStructure,
									  { ConstantParameter( a, NoConstraint() ) },
									  convertParamVector( { sigma_node }, { initial_sigma } ),
									  Matrix( 1, 1, 1 ) )
			, GPPPCMRPCV( { sigma_node } )
		{}
	};
}

#endif // !CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCMR_PCV_HPP