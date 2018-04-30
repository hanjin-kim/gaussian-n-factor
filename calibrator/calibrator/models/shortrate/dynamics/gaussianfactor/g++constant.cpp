#include <calibrator/models/shortrate/dynamics/gaussianfactor/g++constant.hpp>

namespace HJCALIBRATOR
{
	Real GPPConstantDynamics::meanTforward( Size i, Size j, Time T, Time s, Time t ) const
	{
		Real ai = a( i, 0.0 );
		Real aj = a( j, 0.0 );
		Real sigmai = sigma( i, 0.0 );
		Real sigmaj = sigma( j, 0.0 );

		Real integral = (1 - exp( -ai * (t - s) )) / ai / aj
			- (exp( -aj * (T - t) ) - exp( -aj * T - ai * t + (ai + aj)*s ) ) / (ai + aj) / aj;

		Real val = sigmai * sigmaj * integral;
		return val;
	}

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

		Real val = sigma( i, 0.0 )*sigma( j, 0.0 )*(1 - exp( - asum*(t - s) )) / asum;
		return val;
	}
	Real GPPConstantDynamics::phi( Size i, Size j, Time t ) const
	{
		Real a_i = a( i )(0.0);
		Real a_j = a( j )(0.0);

		Real sigma_i = sigma( i )(0.0);
		Real sigma_j = sigma( j )(0.0);

		Real integral = (1 - exp( -a_i * t )) / a_i / a_i
			- (1 - exp( -(a_i + a_j) * t )) / a_i / a_j
			+ (1 - exp( -a_j * t )) / a_j / a_j;

		return sigma_i * sigma_j * integral;
	}
}