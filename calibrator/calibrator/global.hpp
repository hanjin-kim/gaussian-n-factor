#ifndef CALIBRATOR_GLOBAL_HPP
#define CALIBRATOR_GLOBAL_HPP

#include <boost/shared_ptr.hpp>

namespace HJCALIBRATOR
{
	using namespace QuantLib;

	template <typename T>
	using shared_ptr = boost::shared_ptr<T>;
}

#endif // !CALIBRATOR_GLOBAL_HPP
