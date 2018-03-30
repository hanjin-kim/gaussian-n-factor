#ifndef CALIBRATOR_MODELS_PARAMETERS_INTEGRABLECONSTANT_HPP
#define CALIBRATOR_MODELS_PARAMETERS_INTEGRABLECONSTANT_HPP

#include <calibrator/models/integrableparameter.hpp>

namespace HJCALIBRATOR
{
	class IntegrableConstantParameter : public IntegrableParameter
	{
	private:
		class Impl : public IntegrableParameter::Impl
		{
		public:
			Real value( const Array& params, Time ) const {
				return params[0];
			}
			Real integral( const Array& params, Time t0, Time t1 ) const {
				return params[0] * (t1 - t0);
			}
		};
	public:
		IntegrableConstantParameter( const Constraint& constraint )
			: IntegrableParameter( 1,
								   boost::shared_ptr<Impl>( new Impl ),
								   constraint )
		{}

		IntegrableConstantParameter( Real value,
									 const Constraint& constraint )
			: IntegrableParameter( 1,
								   boost::shared_ptr<Impl>( new Impl ),
								   constraint ) 
		{
			params_[0] = value;
			QL_REQUIRE( testParams( params_ ),
						value << ": invalid value" );
		}
	};
}



#endif // !CALIBRATOR_MODELS_PARAMETERS_INTEGRABLECONSTANT_HPP
