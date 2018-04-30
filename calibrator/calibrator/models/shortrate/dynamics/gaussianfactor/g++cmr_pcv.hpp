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

		virtual Real meanTforward( Size i, Size j, Time T, Time s, Time t ) const override;
		virtual Real integralVariance( Size i, Size j, Time s, Time t ) const override;
		virtual Real variance( Size i, Size j, Time s, Time t ) const override;

	protected:
		GPPPCMRPCV( const std::vector<RealVector>& sigma_nodes )
		{
			combineNodes( sigma_nodes );
		}

		virtual Real phi( Size i, Size j, Time t ) const override;

	private:
		void combineNodes( const std::vector<RealVector>& sigma_nodes );
	};

	class G1PPPCMRPCV : public Gaussian1FactorDynamics, public GPPPCMRPCV
	{
	public:
		G1PPPCMRPCV( const Handle<YieldTermStructure>& termStructure,
					 Real a,
					 const RealVector& sigma_node,
					 const RealVector& initial_sigma )
			: GaussianFactorDynamics( termStructure,
									  { ConstantParameter( a, PositiveConstraint() ) },
									  convertParamVector( { sigma_node }, { initial_sigma } ),
									  Matrix( 1, 1, 1))
			, GPPPCMRPCV( { sigma_node } )
		{}

		virtual ~G1PPPCMRPCV() {}
	};

	class G2PPPCMRPCV : public Gaussian2FactorDynamics, public GPPPCMRPCV
	{
	public:
		G2PPPCMRPCV( const Handle<YieldTermStructure>& termStructure,
					 Real a,
					 const RealVector& sigma_node,
					 const RealVector& initial_sigma,
					 Real b,
					 const RealVector& eta_node,
					 const RealVector& initial_eta,
					 Real rho )
			: GaussianFactorDynamics( termStructure,
									   { ConstantParameter( a, PositiveConstraint() ), ConstantParameter( b, PositiveConstraint() ) },
									   convertParamVector( { sigma_node, eta_node }, { initial_sigma, initial_eta } ),
									   getCorrelationMatrix( rho ) )
			, GPPPCMRPCV( { sigma_node, eta_node } )
		{}

		virtual ~G2PPPCMRPCV() {}
	};
}

#endif // !CALIBRATOR_MODELS_SHORTRATE_DYNAMICS_GAUSSIANFACTOR_GPPCMR_PCV_HPP