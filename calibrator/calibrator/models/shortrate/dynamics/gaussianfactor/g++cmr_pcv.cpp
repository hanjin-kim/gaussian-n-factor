#include <calibrator/models/shortrate/dynamics/gaussianfactor/g++cmr_pcv.hpp>

namespace HJCALIBRATOR
{
	ParamVector convertParamVector( const std::vector<RealVector>& sigma_node,
									const std::vector<RealVector>& initial_sigma )
	{
		QL_REQUIRE( sigma_node.size() == initial_sigma.size(),
					"Sigma parameter dimesion mismatched." );

		ParamVector rv;

		for ( Size i = 0; i < sigma_node.size(); i++ )
		{
			const RealVector& nodes = sigma_node[i];
			const RealVector& initval = initial_sigma[i];

			QL_REQUIRE( nodes.size() + 1 == initval.size(),
						"Requirement not met for " << i << "-th sigma parameter "
						<< ": node size + 1 == init value size" );

			PiecewiseConstantParameter tmpparam( nodes, PositiveConstraint() );

			Real a = tmpparam( 1.2 );

			for ( Size j = 0; j < initval.size(); j++ )
			{
				tmpparam.setParam( j, initval[j] );
			}

			Real b = tmpparam( 1.2 );

			rv.push_back( tmpparam );
		}

		return rv;
	}

	Real GPPPCMRPCV::meanTforward( Size i, Size j, Time T, Time s, Time t ) const
	{
		Real ai = a( i, 0.0 );
		Real aj = a( j, 0.0 );

		const Parameter& sigma_i = sigma( i );
		const Parameter& sigma_j = sigma( j );

		RealVector::const_iterator it_nodes = combined_nodes_[i][j].begin();
		while ( it_nodes != combined_nodes_[i][j].end() && *it_nodes < s )
		{
			it_nodes++;
		}

		Real intsum = 0;

		Real begin = s;

		while ( it_nodes != combined_nodes_[i][j].end() && *it_nodes < t )
		{
			Real end = *it_nodes;
			Real mid = (end + begin) / 2.;

			Real integral = (1 - exp( -ai * (end - begin) )) / ai / aj
				- (exp( -aj * (T - end) ) - exp( -aj * T - ai * end + (ai + aj)*begin )) / (ai + aj) / aj;

			intsum += sigma_i( mid )*sigma_j( mid )*integral;
			begin = end;
			it_nodes++;
		}

		return intsum;
	}

	Real GPPPCMRPCV::integralVariance( Size i, Size j, Time s, Time t ) const
	{
		Real ai = a( i, 0.0 );
		Real aj = a( j, 0.0 );

		const Parameter& sigma_i = sigma( i );
		const Parameter& sigma_j = sigma( j );

		RealVector::const_iterator it_nodes = combined_nodes_[i][j].begin();
		while ( it_nodes != combined_nodes_[i][j].end() && *it_nodes < s )
		{
			it_nodes++;
		}

		Real intsum = 0;

		Real begin = s;
		
		while ( it_nodes != combined_nodes_[i][j].end() && *it_nodes < t )
		{
			Real end = *it_nodes;
			Real mid = (end + begin) / 2.;
			Real dt = end - begin;

			Real sigmai = sigma_i( mid );
			Real sigmaj = sigma_j( mid );
			Real c = sigmai * sigmaj / ai / aj;


			intsum += c * (dt
							+ (1 - exp( -(ai + aj) * dt )) / (ai + aj)
							- (1 - exp( -ai * dt )) / ai
							- (1 - exp( -aj * dt )) / aj);
			begin = end;
			it_nodes++;
		}

		Real end = t;
		Real mid = (t + begin) / 2.;
		Real dt = end - begin;

		Real sigmai = sigma_i( mid );
		Real sigmaj = sigma_j( mid );
		Real c = sigmai * sigmaj / ai / aj;


		intsum += c * (dt
						+ (1 - exp( -(ai + aj) * dt )) / (ai + aj)
						- (1 - exp( -ai * dt )) / ai
						- (1 - exp( -aj * dt )) / aj);

		return intsum;
	}

	Real GPPPCMRPCV::variance( Size i, Size j, Time s, Time t ) const
	{
		Real asum = a( i, 0.0 ) + a( j, 0.0 );

		const Parameter& sigma_i = sigma( i );
		const Parameter& sigma_j = sigma( j );

		RealVector::const_iterator it_nodes = combined_nodes_[i][j].begin();
		while ( it_nodes != combined_nodes_[i][j].end() && *it_nodes < s )
		{
			it_nodes++;
		}

		Real intsum = 0;

		Real begin = s;
		while ( it_nodes != combined_nodes_[i][j].end() && *it_nodes < t )
		{
			Real end = *it_nodes;
			Real mid = (end + begin) / 2.;
			intsum += sigma_i( mid )*sigma_j( mid )*(1 - exp( -asum * (end - begin) )) / asum;
			begin = end;
			it_nodes++;
		}
		
		Real end = t;
		Real mid = (t + begin) / 2.;

		intsum += sigma_i( mid )*sigma_j( mid )*(1 - exp( -asum * (end - begin) )) / asum;

		return intsum;
	}

	void GPPPCMRPCV::combineNodes( const std::vector<RealVector>& sigma_nodes )
	{
		Size dim = dimension();

		combined_nodes_.resize( dim, std::vector<RealVector>( dim ) );

		for ( Size i = 0; i < dim; i++ )
		{
			const RealVector nodes_i = sigma_nodes[i];

			for (Size j = 0; j < dim; j++ )
			{
				const RealVector nodes_j = sigma_nodes[j];

				RealVector& nodeij = combined_nodes_[i][j];
				nodeij.reserve( nodes_i.size() + nodes_j.size() );
				nodeij.insert( nodeij.end(), nodes_i.begin(), nodes_i.end() );
				nodeij.insert( nodeij.end(), nodes_i.begin(), nodes_i.end() );

				std::sort( nodeij.begin(), nodeij.end() );
				nodeij.erase( unique( nodeij.begin(), nodeij.end() ), nodeij.end() );
			}
		}
	}
}