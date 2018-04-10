#ifndef CALIBRATOR_PARAMETERS_INTEGRABLEPIECEWISECONSTANT_HPP
#define CALIBRATOR_PARAMETERS_INTEGRABLEPIECEWISECONSTANT_HPP

#include <calibrator/models/integrableparameter.hpp>

namespace HJCALIBRATOR
{
	class IntegrablePiecewiseConstant : IntegrableParameter
	{
	private:
		class Impl : public IntegrableParameter::Impl
		{
		public:
			Impl( const std::vector<Time> nodes = std::vector<Time>() )
				: nodes_( nodes )
			{}

			Real value( const Array& params, Time t ) const 
			{
				Size idx = 0;
				for( Time node : nodes_ )
				{ 
					if ( t > node ) idx++;
					else break;
				}
				
				return params[idx];
			}
			Real integral( const Array& params, Time t0, Time t1 ) const 
			{
				return integral( params, t1 ) - integral( params, t0 );
			}

		private:
			Real integral( const Array& params, const Time t ) const
			{
				Size idx = 0;
				Real sum = 0;
				for ( Size i = 0; i < nodes_.size(); i++ )
				{
					if ( t >= nodes_[i] )
					{
						sum += params[i] * (nodes_[i] - (i > 0 ? nodes_[i - 1] : 0));
					}
					else
					{
						sum += params[i] * (t - (i > 0 ? nodes_[i - 1] : 0));
						break;
					}
				}

				return sum;
			}

			std::vector<Time> nodes_;
			std::vector<Time> diff_t_;
		};
	public:
		IntegrablePiecewiseConstant( const Constraint& constraint )
			: IntegrableParameter( 1,
								   boost::shared_ptr<Impl>( new Impl ),
								   constraint )
		{}

		IntegrablePiecewiseConstant( const std::vector<Real>& values,
									 const std::vector<Time>& nodes,
									 const Constraint& constraint )
			: IntegrableParameter( values.size(),
								   boost::shared_ptr<Impl>( new Impl( nodes ) ),
								   constraint )
		{
			for ( Size i = 0; i < values.size(); i++ )
			{
				params_[i] = values[i];
			}
			QL_REQUIRE( testParams( params_ ), "invalid value" );
		}
	};
}

#endif // !CALIBRATOR_PARAMETERS_INTEGRABLEPIECEWISECONSTANT_HPP

