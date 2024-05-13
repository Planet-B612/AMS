#pragma once
#include "../../dSFMT/dSFMT.h"
#include "CommonStruc.h"
#include <memory>
#include <queue>
#include "CommonFunc.h"
#include <cmath>
using namespace std;
class mRRsets
{
private:
	/// __numV: number of nodes in the graph.
	uint32_t __numV;
	/// __numE: number of edges in the graph.
	size_t __numE = 0;
	/// __numRRsets: number of RR sets.
	size_t __numRRsets = 0;
	vector<bool> __vecVisitBool;
	Nodelist __vecVisitNode;
	Nodelist __roots;
	vector<vector<uint32_t>> PosReal;  // Can be replaced by directly revising the input Reverse Graph.
	vector<bool> __Activated;
	uint32_t __numV_left;
	uint32_t __eta_left;

	/// Initialization
	void init_RRcollection()
	{
		__numV = (uint32_t)_graph.size();  
		// for (auto& nbrs : _graph) __numE += nbrs.size();
		_FRsets = FRsets(__numV);
		__vecVisitBool = std::vector<bool>(__numV);
		__vecVisitNode = Nodelist(__numV);
		PosReal=vector<vector<uint32_t>>(__numV);
		__Activated=vector<bool>(__numV, false);
	}

public:
	/// _graph: reverse graph
	Graph & _graph;  // the 'const' modifier is deleted
	/// _FRsets: forward cover sets, _FRsets[i] is the node sets that node i can reach (in reverse graph?)
	/// i.e., the node in RR-sets of _FRsets can influence i in the original graph, i can reach these nodes/RR-sets in the reverse graph.)
	FRsets _FRsets;
	/// _RRsets: reverse cover sets, _RRsets[i] is the node set that can reach node i
	// RRsets _RRsets;
	/// _cascadeModel: the cascade model, default is IC
	//FRsets _FRsets_real; // Forward sets under the same realization

	string _cascadeModel;

	explicit mRRsets(Graph & graph) : _graph(graph) // the 'const' modifier is deleted too.
	{
		init_RRcollection();
	}

	/// Set cascade model
	void set_cascade_model(const string model)
	{
		_cascadeModel = model;
	}

	/// Returns the number of nodes in the graph.
	uint32_t get_nodes() const
	{
		return __numV;
	}

	/// Returns the number of edges in the graph.
	// size_t get_edges() const
	// {
	// 	return __numE;
	// }

	/// Returns the number of RR sets in the graph.
	size_t get_RRsets_num() const
	{
		return __numRRsets;
	}

	/// Get out degree in the original graph
	std::vector<size_t> get_out_degree() const
	{
		std::vector<size_t> outDeg(__numV);
		for (auto& nbrs : _graph)
		{
			for (auto& nbr : nbrs)
			{
				outDeg[get<0>(nbr)]++;
			}
		}
		return outDeg;
	}

	/// Generate a set of n mRR sets
	void build_n_mRRsets(const size_t numSamples)
	{
		if (numSamples > SIZE_MAX)
		{
			std::cout << "Error: the number of RRsets you need is too large" << std::endl;
			exit(1);
		}
		const auto prevSize = __numRRsets;
		__numRRsets = __numRRsets > numSamples ? __numRRsets : numSamples;
		//cout<<"__numRRsets is "<<__numRRsets<<endl;
		double decimal=1.0*__numV_left/__eta_left;
		uint32_t root_num=floor(decimal);
		double prob=decimal-root_num;
		for (auto i = prevSize; i < numSamples; i++)
		{
			// cout<<"----------Generating the "<<i<<"-th RRsets."<<endl;
			build_one_mRRset(i, root_num, prob);
		}
		//cout<<"---------Finished RRsets generation----------"<<endl;
	}

	uint32_t build_one_mRRset(uint32_t RRid, uint32_t root_num, double prob)
	{
		size_t numVisitNode = 0, currNode = 0; 
		root_num += (dsfmt_gv_genrand_open_close() <= prob);		
		while (root_num > 0)	
		{
			uint32_t node = dsfmt_gv_genrand_uint32_range(__numV);
			while (__Activated[node] || __vecVisitNode[node]) node = dsfmt_gv_genrand_uint32_range(__numV);
			__roots.push_back(node);
			--root_num;
			__vecVisitBool[node]=true;
			__vecVisitNode[numVisitNode++]=node;
			_FRsets[node].push_back(RRid);
		}
		while (currNode < numVisitNode)
		{
			const auto expand = __vecVisitNode[currNode++];
			if(_cascadeModel=="IC"||"ic")
			{
				for (auto& nbr : _graph[expand])
				{
					const auto nbrId = get<0>(nbr);
					if (__vecVisitBool[nbrId]||__Activated[nbrId])
						continue;
					const auto randDouble = dsfmt_gv_genrand_open_close();
					if (randDouble > get<1>(nbr))
						continue;
					__vecVisitNode[numVisitNode++] = nbrId;
					//cout<<"The visited node is "<<nbrId<<", ";
					__vecVisitBool[nbrId] = true;
					_FRsets[nbrId].push_back(RRid);
				}
			}
			else if(_cascadeModel=="LT"||"lt")
			{
				if(_graph[expand].size()==0)
				continue;
				uint32_t index=dsfmt_gv_genrand_uint32_range(_graph[expand].size());
				const auto nbrId=get<0>(_graph[expand][index]);
				if (__vecVisitBool[nbrId]||__Activated[nbrId])
					continue;
				__vecVisitNode[numVisitNode++] = nbrId;
				//cout<<"The visited node is "<<nbrId<<", ";
				__vecVisitBool[nbrId] = true;
				_FRsets[nbrId].push_back(RRid);
			}
			else
			assert(false);
		}
		for (size_t i = 0; i < numVisitNode; i++) 
		{
			__vecVisitBool[__vecVisitNode[i]] = false;
			//cout<<__vecVisitNode[i]<<", ";
		}
		//cout<<endl;
		//_RRsets.push_back(RRset(__vecVisitNode.begin(), __vecVisitNode.begin() + numVisitNode));
		return 0;
	}
	/// Refresh the RRsets
	void refresh_RRsets()
	{
		// for (auto i = __numRRsets; i--;)
		// {
		// 	RRset().swap(_RRsets[i]);
		// }
		// RRsets().swap(_RRsets);
		for (auto i = __numV; i--;)
		{
			FRset().swap(_FRsets[i]);
		}
		__numRRsets = 0;
	}

	/// Release memory
	void release_memory()
	{
		refresh_RRsets();
		std::vector<bool>().swap(__vecVisitBool);
		Nodelist().swap(__vecVisitNode);
		FRsets().swap(_FRsets);
		Nodelist().swap(__roots);
	}

	/// Generate one node with probabilities according to their weights for the LT cascade model
	static inline size_t gen_random_node_by_weight_LT(const Edgelist& edges)
	{
		const double weight = dsfmt_gv_genrand_open_close();
		size_t minIdx = 0, maxIdx = edges.size() - 1;
		if (weight < get<1>(edges.front())) return 0; // First element
		if (weight > get<1>(edges.back())) return edges.size() + 1; // No element
		while (maxIdx > minIdx)
		{
			const size_t meanIdx = (minIdx + maxIdx) / 2;
			const auto meanWeight = get<1>(edges[meanIdx]);
			if (weight <= meanWeight) maxIdx = meanIdx;
			else minIdx = meanIdx + 1;
		}
		return maxIdx;
	}
	
};

using TmRRsets = mRRsets;
using PmRRsets = std::shared_ptr<TmRRsets>;
