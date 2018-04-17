#include <calibrator/models/shortrate/dynamics/gaussianfactor/g++constant.hpp>

namespace HJCALIBRATOR
{
	Real GPPConstantDynamics::integralVariance( Size i, Size j, Time s, Time t ) const
	{
		Time dt = t - s;
		Real ai = a( i, 0.0 );
		Real aj = a( j, 0.0 );
		Real sigmai = sigma( i, 0.0 );
		Real sigmaj = sigma( j, 0.0 );
		Real c = sigmai * sigmaj / ai / aj;

		return c * (dt
					 + (1 - exp( -(ai + aj) * dt )) / (ai + aj)
					 - (1 - exp( -ai * dt )) / ai
					 - (1 - exp( -aj * dt )) / aj);
	}

	Real GPPConstantDynamics::variance( Size i, Size j, Time s, Time t ) const
	{
		Real asum = a( i, 0.0 ) + a( j, 0.0 );

		return sigma( i, 0.0 )*sigma( j, 0.0 )*(1 - exp( - asum*(t - s) )) / asum;
	}
}