#ifndef CALIBRATOR_PARAMETERS_LOGISTIC_FUNCTION_HPP
#define CALIBRATOR_PARAMETERS_LOGISTIC_FUNCTION_HPP

#include <ql/models/parameter.hpp>

#include <calibrator/global.hpp>

namespace HJCALIBRATOR
{
	// Logistic function expression as a Quantlib::Parameter
	// It takes three parameters,
	// p0 : the curve's maximum value
	// p1 : the steepness of the curve
	// p2 : the sigmoid's midpoint
	// 
	// f(x) = p0 / ( 1 + exp( -p1 * ( x - p2 ) ) )
	//	
	class LogisticFunctionParameter : Parameter
	{
		class Impl : public Parameter::Impl
		{
		public :
			Real value( const Array& params, Time t ) const override
			{
				return params[0] / (1 + exp( -params[1] * (t - params[2]) ));
			}
		};
		LogisticFunctionParameter(const Array& params, const Constraint& constraint = NoConstraint())
			: Parameter( params.size(), 
						 boost::shared_ptr<Parameter::Impl>( new LogisticFunctionParameter::Impl() ), constraint )
		{}
	};
}

#endif // !CALIBRATOR_PARAMETERS_LOGISTIC_FUNCTION_HPP

