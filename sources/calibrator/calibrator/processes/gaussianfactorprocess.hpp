#ifndef HJCALIBRATOR_PROCESSES_GAUSSIANFACTORPROCESS_HPP
#define HJCALIBRATOR_PROCESSES_GAUSSIANFACTORPROCESS_HPP

#include <ql/stochasticprocess.hpp>

#include <calibrator/models/shortrate/dynamics/gaussianfactordynamics.hpp>

namespace HJCALIBRATOR
{
	class GaussianFactorProcess : public StochasticProcess
	{
		shared_ptr<GaussianFactorDynamics> dynamics_;

	public:		
		GaussianFactorProcess( shared_ptr<GaussianFactorDynamics> dynamics, Array x0 );

		//! \name StochasticProcess interface
		//@{
		Size size() const;
		Disposable<Array> initialValues() const;
		Disposable<Array> drift( Time t, const Array& x ) const;
		Disposable<Matrix> diffusion( Time t, const Array& x ) const;
		Disposable<Array> expectation( Time t0, const Array& x0, Time dt ) const;
		Disposable<Matrix> stdDeviation( Time t0, const Array& x0,
										 Time dt ) const;
		Disposable<Matrix> covariance( Time t0, const Array& x0, Time dt ) const;

		Array x0_;
	};
}

#endif // !HJCALIBRATOR_PROCESSES_GAUSSIANFACTORPROCESS_HPP
