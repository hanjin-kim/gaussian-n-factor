#ifndef CALIBRATOR_MODELS_PARAMETERS_SHORTRATE_GAUSSIANFACTORFITTINGPARAMETER_HPP
#define CALIBRATOR_MODELS_PARAMETERS_SHORTRATE_GAUSSIANFACTORFITTINGPARAMETER_HPP

#include <ql/math/matrix.hpp>
#include <ql/math/integrals/kronrodintegral.hpp>
#include <ql/math/integrals/simpsonintegral.hpp>
#include <ql/math/randomnumbers/sobolrsg.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/models/parameter.hpp>

#include <calibrator/global.hpp>

namespace HJCALIBRATOR
{
	using RealVector = std::vector<Real>;
	using TimeVector = std::vector<Time>;

	using ParamVector = std::vector<Parameter>;
	using ParamMatrix = std::vector<std::vector<Parameter>>;

	using RefParam = std::reference_wrapper<Parameter>;
	using RefParamVector = std::vector<RefParam>;
	using RefParamMatrix = std::vector<std::vector<RefParam>>;

	class GaussianFactorDynamics
	{
	public:
		GaussianFactorDynamics( const Handle<YieldTermStructure>& termStructure,
								const ParamVector& a,
								const ParamVector& sigma,
								const Matrix& rho = Matrix() );

		virtual ~GaussianFactorDynamics() {}

	protected:
		GaussianFactorDynamics() // for virtual inheritance
			: integrator_( GaussKronrodAdaptive( 0.01, 10000 ) )
		{}

	public:
		Size dimension() const { return a_.size(); }

		Parameter a( Size i ) const { return a_[i]; }
		Parameter sigma( Size i ) const { return sigma_[i]; }
		Parameter rho( Size i, Size j ) const;

		Real a( Size i, Time t ) const { return (a_[i])( t ); }
		Real sigma( Size i, Time t ) const { return (sigma_[i])( t ); }

		void a( const Parameter& a, Size i );
		void sigma( const Parameter& sigma, Size i );
		void rho( const Parameter& rho, Size i, Size j );

		Handle<YieldTermStructure> termStructure() const { return termStructure_; }

		Real E( Size i, Time s, Time t ) const { return E( a_[i], s, t ); }
		Real B( Size i, Time s, Time t ) const { return B( a_[i], s, t ); }
		Real meanTforward( Size i, Time T, Time s, Time t ) const;

		virtual Real meanTforward( Size i, Size j, Time T, Time s, Time t ) const;
		virtual Real integralVariance( Size i, Size j, Time s, Time t ) const;
		virtual Real variance( Size i, Size j, Time s, Time t ) const;

		virtual Real A( Time t, Time T ) const;
		
		virtual Real phi( Time t ) const;
		virtual Real variance( Time s, Time t ) const;
		virtual Real integralVariance( Time s, Time t ) const;

	protected:
		virtual Real E( const Parameter& a, Time s, Time t ) const;
		virtual Real E( const Parameter& ai, const Parameter& aj, Time s, Time t ) const;
		virtual Real B( const Parameter& a, Time s, Time t ) const;
		virtual Real B( const Parameter& ai, const Parameter& aj, Time s, Time t ) const;

		virtual Real phi( Size i, Size j, Time t ) const;

		GaussKronrodAdaptive integrator_;
		//SimpsonIntegral integrator_;
	private:
		void setupCorrelMatrix( const Matrix& rho );

		Handle<YieldTermStructure> termStructure_;
		ParamVector a_;
		ParamVector sigma_;
		std::map<std::pair<Size, Size>, Parameter> rho_;
	};

	
	class Gaussian1FactorDynamics : public virtual GaussianFactorDynamics
	{
	public:
		Gaussian1FactorDynamics( const Handle<YieldTermStructure>& termStructure,
								 const Parameter& a,
								 const Parameter& sigma )
			: GaussianFactorDynamics( termStructure, { a }, { sigma }, Matrix( 1, 1, 1 ) )
		{}

		virtual ~Gaussian1FactorDynamics() {}

	protected:
		Gaussian1FactorDynamics() {} // for virtual inheritance
	};

	class Gaussian2FactorDynamics : public virtual GaussianFactorDynamics
	{
	public:
		Gaussian2FactorDynamics( const Handle<YieldTermStructure>& termStructure,
								 const Parameter& a,
								 const Parameter& sigma,
								 const Parameter& b,
								 const Parameter& eta,
								 Real rho )
			: GaussianFactorDynamics( termStructure, { a, b }, { sigma, eta }, getCorrelationMatrix(rho) )
		{}

		virtual ~Gaussian2FactorDynamics() {}

	protected:
		Gaussian2FactorDynamics() {} // for virtual inheritance
		
		Matrix getCorrelationMatrix( Real rho )
		{
			Matrix ret( 2, 2, 1 );
			ret[0][1] = rho;
			ret[1][0] = rho;

			return ret;
		}
	};
	

	inline Real GaussianFactorDynamics::E( const Parameter& ai, const Parameter& aj, Time s, Time t ) const
	{
		return E( ai, s, t ) * E( aj, s, t );
	}

	inline Real GaussianFactorDynamics::B( const Parameter& ai, const Parameter& aj, Time s, Time t ) const
	{
		return B( ai, s, t ) * B( aj, s, t );
	}
}

#endif // !CALIBRATOR_MODELS_PARAMETERS_SHORTRATE_GAUSSIANFACTORFITTINGPARAMETER_HPP
