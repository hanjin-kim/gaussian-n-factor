#ifndef CALIBRATOR_PARAMETERS_ANALYTICPARAMETER_HPP
#define CALIBRATOR_PARAMETERS_ANALYTICPARAMETER_HPP

#include <ql/models/parameter.hpp>

#include <calibrator/global.hpp>

namespace HJCALIBRATOR
{
	//! Base class for integrable model parameters
	class IntegrableParameter : public Parameter
	{
	protected :
		//! Base class for integrable model parameter implementation
		class Impl : public Parameter::Impl
		{
		public :
			virtual ~Impl() {}
			virtual Real value( const Array& params, Time t ) const = 0;
			virtual Real integral( const Array& params, Time t0, Time t1 ) const = 0;
		};
	public :
		IntegrableParameter()
			: Parameter()
		{}

		IntegrableParameter( const Parameter& p )
			: Parameter( p )
		{}

		virtual Real integral( Time t0, Time t1 ) const
		{
			boost::shared_ptr<Impl> impl = boost::static_pointer_cast<Impl>(impl_);
			return impl->integral( params_, t0, t1 );
		}

	protected :
		IntegrableParameter( Size size,
							 const boost::shared_ptr<Impl>& impl,
							 const Constraint& constraint )
			: Parameter( size, impl, constraint )
		{}
	};
}

#endif // !CALIBRATOR_PARAMETERS_ANALYTICPARAMETER_HPP
