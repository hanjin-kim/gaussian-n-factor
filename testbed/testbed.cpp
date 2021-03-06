// testbed.cpp: 콘솔 응용 프로그램의 진입점을 정의합니다.
//

#include "stdafx.h"

#include <iostream>
#include <ctime>

#include <boost/timer.hpp>

#include <ql/instruments/swaption.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/models/calibrationhelper.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/pricingengines/swaption/jamshidianswaptionengine.hpp>
#include <ql/pricingengines/swaption/treeswaptionengine.hpp>
#include <ql/models/shortrate/twofactormodels/g2.hpp>
#include <ql/pricingengines/swaption/g2swaptionengine.hpp>

#include <ql/models/shortrate/calibrationhelpers/swaptionhelper.hpp>
#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/auto_link.hpp>

#include <calibrator/models/shortrate/onefactormodels/generalg1.hpp>
#include <calibrator/models/shortrate/twofactormodels/generalg2.hpp>


#include <calibrator/models/shortrate/dynamics/gaussianfactor/g++constantmeanreversion.hpp>
#include <calibrator/models/shortrate/dynamics/gaussianfactor/g++constant.hpp>
#include <calibrator/models/shortrate/dynamics/gaussianfactor/g++cmr_pcv.hpp>
#include <calibrator/pricingengines/swaption/generalg2swaptionengine.hpp>


using namespace std;
using namespace HJCALIBRATOR;

/*

Size numRows = 15;
Size numCols = 15;

Integer swapLenghts[] = {
	1,     2,     3,     4,     5,	6,	7,	8,	9,	10, 12, 15, 20, 25, 30 }; //swap ������ ����
Volatility swaptionVols[] = { //swaption vol ���� �ɼǸ���, ���� �ɼǸ��� 
	0.20565,0.2153,0.21325,0.20675,0.20585,0.20085,0.19595,0.19215,0.18845,0.18465,0.18465,0.18465,0.18465,0.18465,0.18465,
	0.20175,0.21225,0.20555,0.20035,0.1994,0.1941,0.1888,0.1838,0.1788,0.1738,0.1738,0.1738,0.1738,0.1738,0.1738,
	0.2002,0.20725,0.1974,0.1937,0.1903,0.1859,0.1815,0.1776,0.1737,0.1698,0.1698,0.1698,0.1698,0.1698,0.1698,
	0.2008,0.1964,0.1926,0.1878,0.1842,0.1783,0.1724,0.1706,0.1687,0.1669,0.1669,0.1669,0.1669,0.1669,0.1669,
	0.1987,0.1912,0.1873,0.1839,0.1818,0.1734,0.1651,0.1655,0.1658,0.1662,0.1662,0.1662,0.1662,0.1662,0.1662,
	0.1989,0.1913,0.18625,0.1829,0.1797,0.1728,0.166,0.1662,0.1664,0.1666,0.1666,0.1666,0.1666,0.1666,0.1666,
	0.19905,0.1913,0.1852,0.1819,0.1776,0.1722,0.1669,0.1669,0.1669,0.1669,0.1669,0.1669,0.1669,0.1669,0.1669,
	0.20055,0.19265,0.1857,0.1831,0.1784,0.1742,0.17,0.1696,0.1691,0.1687,0.1687,0.1687,0.1687,0.1687,0.1687,
	0.20205,0.19405,0.1862,0.1844,0.1792,0.1762,0.1732,0.1723,0.1714,0.1705,0.1705,0.1705,0.1705,0.1705,0.1705,
	0.20355,0.1954,0.1867,0.1856,0.18,0.1781,0.1763,0.175,0.1736,0.1723,0.1723,0.1723,0.1723,0.1723,0.1723,
	0.20485,0.199,0.1903,0.1893,0.1838,0.1811,0.1785,0.1788,0.179,0.1793,0.1793,0.1793,0.1793,0.1793,0.1793,
	0.20615,0.2004,0.192,0.191,0.1855,0.1861,0.1866,0.1874,0.1881,0.1889,0.1889,0.1889,0.1889,0.1889,0.1889,
	0.20995,0.19955,0.1943,0.1896,0.1876,0.1872,0.1869,0.1861,0.1853,0.185,0.1845,0.1845,0.1845,0.1845,0.1845,
	0.2133,0.2056,0.19465,0.19,0.1882,0.1878,0.187,0.1869,0.1863,0.186,0.1858,0.1858,0.1858,0.1858,0.1858,
	0.21665,0.21165,0.19495,0.191,0.1887,0.1883,0.188,0.1877,0.1874,0.187,0.1872,0.1872,0.1872,0.1872,0.1872 };
	*/

/*
Size numRows = 4; //swap
Size numCols = 3; // option

Integer swapLenghts[] = { 1,2,5,10 };
Volatility swaptionVols[] = { // Swap X option
	0.095, 0.089, 0.073,
	0.089, 0.082, 0.067,
	0.072, 0.067, 0.053,
	0.056, 0.051, 0.040 };
*/

Size numRows = 7; 
Size numCols = 10; 

Integer swapLenghts[] = { 1,2,3,4,5,6,7,8,9,10};

Volatility swaptionVols[] = {
	0.164,	0.155,	0.143,	0.131,	0.124,	0.119,	0.116,	0.112,	0.11,	0.107,
	0.16,	0.15,	0.139,	0.129,	0.122,	0.119,	0.116,	0.113,	0.11,	0.108,
	0.157,	0.145,	0.134,	0.124,	0.119,	0.115,	0.113,	0.11,	0.108,	0.106,
	0.148,	0.136,	0.126,	0.119,	0.114,	0.112,	0.109,	0.107,	0.105,	0.103,
	0.14,	0.128,	0.121,	0.114,	0.11,	0.107,	0.105,	0.103,	0.102,	0.1,
	0.13,	0.119,	0.113,	0.105,	0.101,	0.099,	0.097,	0.096,	0.095,	0.093,
	0.116,	0.107,	0.1,	0.093,	0.09,	0.089,	0.087,	0.086,	0.085,	0.084
};


int main()
{
	Date todaysDate( 25, April, 2017 ); //���� ��¥
	Calendar calendar = TARGET(); // ��� ���� ó��?
	Date settlementDate( 28, April, 2017 ); //settlementDate
	Settings::instance().evaluationDate() = todaysDate; //evaluationDate
														// flat yield term structure impling 1x10 swap at 2.224%
	boost::shared_ptr<Quote> flatRate( new SimpleQuote( 0.0222400 ) );
	Handle<YieldTermStructure> rhTermStructure(
		boost::shared_ptr<FlatForward>(
			new FlatForward( settlementDate, Handle<Quote>( flatRate ),
							 Actual365Fixed() ) ) ); // YieldTermStructure This abstract class defines the interface of concrete interest rate structures which will be derived from this one

	Frequency fixedLegFrequency = Annual; // fixedLegFrequency
	BusinessDayConvention fixedLegConvention = Unadjusted; //fixedLegConvention
	BusinessDayConvention floatingLegConvention = ModifiedFollowing; //floatingLegConvention
	DayCounter fixedLegDayCounter = Thirty360( Thirty360::European ); //fixedLegDayCounter
	Frequency floatingLegFrequency = Semiannual; //floatingLegFrequency
	VanillaSwap::Type type = VanillaSwap::Payer; //type = VanillaSwap::Payer
	Rate dummyFixedRate = 0.0222400; //FixedRate
	boost::shared_ptr<IborIndex> indexSixMonths( new
												 Euribor6M( rhTermStructure ) );

	Date startDate = calendar.advance( settlementDate, 1, Years,
									   floatingLegConvention );// swap startDate
	Date maturity = calendar.advance( startDate, 10, Years,
									  floatingLegConvention ); // swap maturity
	Schedule fixedSchedule( startDate, maturity, Period( fixedLegFrequency ),
							calendar, fixedLegConvention, fixedLegConvention,
							DateGeneration::Forward, false ); // fixedSchedule
	Schedule floatSchedule( startDate, maturity, Period( floatingLegFrequency ),
							calendar, floatingLegConvention, floatingLegConvention,
							DateGeneration::Forward, false ); // floatSchedule

	boost::shared_ptr<VanillaSwap> swap( new VanillaSwap( //VanillaSwap swap
														  type, 1000.0,
														  fixedSchedule, dummyFixedRate, fixedLegDayCounter,
														  floatSchedule, indexSixMonths, 0.0,
														  indexSixMonths->dayCounter() ) );
	swap->setPricingEngine( boost::shared_ptr<PricingEngine>(
		new DiscountingSwapEngine( rhTermStructure ) ) );
	Rate fixedATMRate = swap->fairRate();
	Rate fixedOTMRate = fixedATMRate * 1.2;
	Rate fixedITMRate = fixedATMRate * 0.8;

	boost::shared_ptr<VanillaSwap> atmSwap( new VanillaSwap(
		type, 1000.0,
		fixedSchedule, fixedATMRate, fixedLegDayCounter,
		floatSchedule, indexSixMonths, 0.0,
		indexSixMonths->dayCounter() ) );
	boost::shared_ptr<VanillaSwap> otmSwap( new VanillaSwap(
		type, 1000.0,
		fixedSchedule, fixedOTMRate, fixedLegDayCounter,
		floatSchedule, indexSixMonths, 0.0,
		indexSixMonths->dayCounter() ) );
	boost::shared_ptr<VanillaSwap> itmSwap( new VanillaSwap(
		type, 1000.0,
		fixedSchedule, fixedITMRate, fixedLegDayCounter,
		floatSchedule, indexSixMonths, 0.0,
		indexSixMonths->dayCounter() ) );

	// defining the swaptions to be used in model calibration
	std::vector<Period> swaptionMaturities; // Period swaptionMaturities
	swaptionMaturities.push_back( Period( 1, Years ) );
	swaptionMaturities.push_back( Period( 2, Years ) );
	swaptionMaturities.push_back( Period( 3, Years ) );
	swaptionMaturities.push_back( Period( 4, Years ) );
	swaptionMaturities.push_back( Period( 5, Years ) );
	//swaptionMaturities.push_back( Period( 6, Years ) );
	swaptionMaturities.push_back( Period( 7, Years ) );
	//swaptionMaturities.push_back( Period( 8, Years ) );
	//swaptionMaturities.push_back( Period( 9, Years ) );
	swaptionMaturities.push_back( Period( 10, Years ) );
	//swaptionMaturities.push_back( Period( 12, Years ) );
	//swaptionMaturities.push_back( Period( 15, Years ) );
	//swaptionMaturities.push_back( Period( 20, Years ) );
	//swaptionMaturities.push_back( Period( 25, Years ) );
	//swaptionMaturities.push_back( Period( 30, Years ) );

	std::vector<boost::shared_ptr<CalibrationHelper> > swaptions; // CalibrationHelper swaptions

																  // List of times that have to be included in the timegrid
	std::list<Time> times;

	Size i;
	for ( i = 0; i<numRows; i++ ) {
		for ( Size j = 0; j < numCols; j++ ) {
			//Size j = numCols - i - 1; // 1x10, 2x8, 3x7, 4x6, 5x5, ...
			Size k = i * numCols + j;
			boost::shared_ptr<Quote> vol( new SimpleQuote( swaptionVols[k] ) );
			swaptions.push_back( boost::shared_ptr<CalibrationHelper>( new
																	   SwaptionHelper( swaptionMaturities[i],
																					   Period( swapLenghts[j], Years ),
																					   Handle<Quote>( vol ),
																					   indexSixMonths,
																					   indexSixMonths->tenor(),
																					   indexSixMonths->dayCounter(),
																					   indexSixMonths->dayCounter(),
																					   rhTermStructure ) ) );
			swaptions.back()->addTimesTo( times ); //Returns a reference to the last element in the vector
		}
	}

	LevenbergMarquardt om;
	//Simplex om( 0.01 );

	boost::shared_ptr<Gaussian2FactorDynamics> dynamics_const( new G2ConstantDynamics( rhTermStructure,
																				 0.2, 0.014,
																				 0.3, 0.014,
																				 -0.75 ) );

	boost::shared_ptr<Gaussian2FactorDynamics> dynamics_CMRPCV( new G2PPPCMRPCV( rhTermStructure,
																		  0.1, {2,5,9.5}, { 0.014,0.014,0.014,0.05 },
																		  0.5, {}, { 0.01 },
																		  -0.75 ) );

	boost::shared_ptr<GeneralizedG2> g2ppmodel_const( new GeneralizedG2( dynamics_const ) );
	boost::shared_ptr<GeneralizedG2> g2ppmodel_CMRPCV( new GeneralizedG2( dynamics_CMRPCV ) );

	//for ( i = 0; i < swaptions.size(); i++ )
	//	swaptions[i]->setPricingEngine( boost::shared_ptr<PricingEngine>(
	//		new GeneralizedG2SwaptionEngine( g2ppmodel_const ) ) );

	//try
	//{
	//	clock_t begin = clock();
	//	g2ppmodel_const->calibrate( swaptions, om,
	//						EndCriteria( 10000, 100, 1.0e-8, 1.0e-8, 1.0e-8 ) );

	//	clock_t end = clock();

	//	cout << "Elapsed time : " << double( end - begin ) / CLOCKS_PER_SEC << endl;
	//	std::cout << "calibrated to:\n"
	//		<< "a     = " << g2ppmodel_const->a()(0) << ", "
	//		<< "sigma = " << g2ppmodel_const->sigma()(0) << "\n"
	//		<< "b     = " << g2ppmodel_const->b()(0) << ", "
	//		<< "eta   = " << g2ppmodel_const->eta()(0) << "\n"
	//		<< "rho   = " << g2ppmodel_const->rho()(0) << ", "
	//		<< std::endl << std::endl;

	//	for ( Size i = 0; i < numRows; i++ ) {
	//		for ( Size j = 0; j < numCols; j++ ) {
	//			/*if ( j != numCols - i - 1 && j != numCols - i - 2 )
	//			{
	//				cout << std::setw( 10 ) << "-";
	//				continue;
	//			}*/

	//			Size k = i * numCols + j;
	//			Real npv = swaptions[k]->modelValue();
	//			Volatility implied = swaptions[k]->impliedVolatility( npv, 1e-4,
	//																  1000, 0.05, 0.50 );
	//			Volatility diff = implied - swaptionVols[k];
	//			cout << std::setw( 10 ) << diff / (swaptionVols[k]) * 100 << "%";
	//		}

	//		cout << endl;
	//	}


	//}
	//catch ( exception& e )
	//{
	//	cout << e.what() << endl;
	//}

	for ( i = 0; i < swaptions.size(); i++ )
		swaptions[i]->setPricingEngine( boost::shared_ptr<PricingEngine>(
			new GeneralizedG2SwaptionEngine( g2ppmodel_CMRPCV ) ) );

	try
	{
		clock_t begin = clock();
		g2ppmodel_CMRPCV->calibrate( swaptions, om,
									EndCriteria( 10000, 100, 1.0e-8, 1.0e-8, 1.0e-8 ) );
		g2ppmodel_CMRPCV->calibrate( swaptions, om,
									 EndCriteria( 10000, 100, 1.0e-8, 1.0e-8, 1.0e-8 ) );

		clock_t end = clock();

		cout << "Elapsed time : " << double( end - begin ) / CLOCKS_PER_SEC << endl;
		std::cout << "calibrated to:\n"
			<< "a     = " << g2ppmodel_CMRPCV->a()(0) << ", "
			<< "sigma = " << g2ppmodel_CMRPCV->sigma()(0) << "\n"
			<< "b     = " << g2ppmodel_CMRPCV->b()(0) << ", "
			<< "eta   = " << g2ppmodel_CMRPCV->eta()(0) << "\n"
			<< "rho   = " << g2ppmodel_CMRPCV->rho()(0) << ", "
			<< std::endl << std::endl;

		for ( Size i = 0; i < numRows; i++ ) {
			for ( Size j = 0; j < numCols; j++ ) {
			/*	if ( j != numCols - i - 1 && j != numCols - i - 2 )
				{
					cout << std::setw( 10 ) << "-";
					continue;
				}*/

				Size k = i * numCols + j;
				Real npv = swaptions[k]->modelValue();
				Volatility implied = swaptions[k]->impliedVolatility( npv, 1e-4,
																	  1000, 0.05, 0.50 );
				Volatility diff = implied - swaptionVols[k];
				cout << std::setw( 10 ) << diff / (swaptionVols[k]) * 100 << "%";
			}

			cout << endl;
		}


	}
	catch ( exception& e )
	{
		cout << e.what() << endl;
	}
 	return 0;
}

