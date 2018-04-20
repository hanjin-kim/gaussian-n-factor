#ifndef HJCALIBRATOR_PRICINGENGINES_SWAPTION_GENERALG2SWAPTIONENGEIN_HPP
#define HJCALIBRATOR_PRICINGENGINES_SWAPTION_GENERALG2SWAPTIONENGEIN_HPP

#include <ql/pricingengines/genericmodelengine.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>

#include <calibrator/models/shortrate/twofactormodels/generalg2.hpp>

namespace HJCALIBRATOR
{
	template <class Model>
	using GenericSwaptionEngine = GenericModelEngine<Model, Swaption::arguments, Swaption::results>;

	class GeneralizedG2SwaptionEngine : public GenericSwaptionEngine<GeneralizedG2>
	{
	public:
		// range is the number of standard deviations to use in the
		// exponential term of the integral for the european swaption.
		// intervals is the number of intervals to use in the integration.
		GeneralizedG2SwaptionEngine( const shared_ptr<GeneralizedG2>& model )
			: GenericSwaptionEngine<GeneralizedG2>( model )
		{}

		void calculate() const {

			QL_REQUIRE( arguments_.settlementType == Settlement::Physical,
						"cash-settled swaptions not priced with G2 engine" );

			// adjust the fixed rate of the swap for the spread on the
			// floating leg (which is not taken into account by the
			// model)
			VanillaSwap swap = *arguments_.swap;
			swap.setPricingEngine( shared_ptr<PricingEngine>(
				new DiscountingSwapEngine( model_->termStructure(), false ) ) );
			Spread correction = swap.spread() *
				std::fabs( swap.floatingLegBPS() / swap.fixedLegBPS() );
			Rate fixedRate = swap.fixedRate() - correction;

			results_.value = model_->swaption( arguments_, fixedRate );
		}
	};
}

#endif // !HJCALIBRATOR_PRICINGENGINES_SWAPTION_GENERALG2SWAPTIONENGEIN_HPP
