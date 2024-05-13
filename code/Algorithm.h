#pragma once
#include "mRRsets.h"
#include "CommonStruc.h"
#include "CommonFunc.h"
#include <memory>
#include <algorithm>
#include "Memory.h"
using namespace std;

class Algorithm
{
	private:
		uint32_t __numV=0;
		//size_t __numE=0;
		size_t num_RRsets = 0;
		mRRsets RR;
		Nodelist Seeds;
		vector<double> __cost;
		double __eta=1;
		double __delta=0.1;
		double __epsilon_Inf=0.1;
		double __delta_Inf=0.01;
		Graph __F_graph;
		Graph __R_graph;
		string __model;
		string __result_dir;
		vector<vector<bool>> __vecCover;  // record whether an FRset_real is covered by some see

	public:

		Algorithm(Graph& graph) : RR(graph)
		{
			__numV = RR.get_nodes();
			//__numE = RR.get_edges();
			__R_graph=graph;
		}
		~Algorithm()
		{ }

void release_mem()
{
	RR.release_memory();
}

vector<uint32_t> max_ratio_lazy(const double Q)
{
	// The first element in tuple is the node ID, the second is the hyper-degree coverage, the third is the ratio, the fourth is the indicator whether the node has been updated in this round.
	vector<tuple<uint32_t, double, double, uint32_t>> ratio(__numV);
	double total_cost=0.0;
	for (uint32_t i = __numV; i--;)
	{
		const uint32_t deg = RR._FRsets[i].size();  // The number of RR-sets covered by i.
		get<0>(ratio[i]) = i;
		get<1>(ratio[i])=deg;
		get<2>(ratio[i])=1.0*deg/__cost[i];
		get<3>(ratio[i])=0;
	}
	make_max_heap(ratio);

	size_t total_deg = 0;  // Record the number of _RRsets covered by seed set.
	vector<bool> RR_Mark(num_RRsets, false);  // check if an edge is removed
	vector<bool> veriRR_state(num_RRsets, false);
	Seeds.clear();
	double expInf=0.0;
	for(uint32_t i=0;i<__eta;i++)  // The i-th seed we are selecting. After we selecting each seed, we should put it to the i-th place.
	{
		// if(std::fmod(i,15)==0) 
		// cout<<"selecting "<<i<<"-th seed"<<endl;
		uint32_t nodeId;//=get<0>(ratio[0]);
		while(get<3>(ratio[0])!=i)  // since we will delete and insert the node, the j-th node would be different.
		{
			nodeId=get<0>(ratio[0]);
			//uint32_t nodeDeg=get<1>(ratio[0]);
			uint32_t nodeDeg=RR._FRsets[nodeId].size();
			for(uint32_t RRId: RR._FRsets[nodeId])
			{
				if(RR_Mark[RRId]==true) 	--nodeDeg;
			}
			//uint32_t num_covered=count(RR._FRsets[nodeId].begin(),RR._FRsets[nodeId].end(),true);
			//nodeDeg=nodeDeg-num_covered;
			tuple<uint32_t, double, double, uint32_t> updated_node=make_tuple(nodeId, nodeDeg, 1.0*nodeDeg/__cost[nodeId],i);
			max_heap_replace_max_value(ratio, updated_node);
		}
		nodeId=get<0>(ratio[0]);
		Seeds.push_back( nodeId );
		total_cost=total_cost+__cost[nodeId];

		if(__model=="1")	// 1: MC
		{
			total_deg=total_deg+get<1>(ratio[0]);  // No prob estimation, only evaluate based on expected influence.
			expInf=1.0*__numV*total_deg/num_RRsets; 
			if(expInf>=389237)  //__eps_MC==lambda
			//if(est_prob>=__prob)  //__eps_MC==lambda
			{
				//cout<<"The total cost is "<<total_cost<<endl;
				cout<<"The expected inf is: "<<expInf<<endl;
				return Seeds;
			}
		}
		
		for(auto rr:RR._FRsets[nodeId])
		{
			if(RR_Mark[rr]) continue;
			RR_Mark[rr]=true;
		}
		tuple<uint32_t, double, double, uint32_t> disable_node=make_tuple(nodeId, -1.0, -1.0, i);
		max_heap_replace_max_value(ratio, disable_node);
	}
	cout<<"Probably an error exists in seed selection, this line should not appear."<<endl;
	cout<<"The expInf is "<<expInf<<endl;
	return {};
	
}

void set_cascade_model(const string model)
{
	__model=model;
	RR.set_cascade_model(model);
	//FF.set_cascade_model(model);
}

void set_parameters(const string model, vector<double> cost, Graph F_graph, double eta, double eps_Inf, double delta_Inf, string result_dir)
{
	RR.set_cascade_model(model);
	__cost=cost;
	__F_graph=F_graph;
	__epsilon_Inf=eps_Inf;
	__delta_Inf=delta_Inf;
	__result_dir=result_dir;
}

void set_F_graph(const Graph F_graph)
{
	__F_graph=F_graph;
}

Nodelist mine()
{
	return {};
}

pair<double, double> actual_prob_eval(const vector<uint32_t> & vecSeed, uint32_t simulations)
{
	//Timer prob_eval_time("prob_eval_time");
	//fstream result_bk(__result_dir, ios::app);
	uint32_t nodeId;
	uint32_t exceed=0; // Record the num of simulations that exceed Q
	queue<uint32_t> Que;
	vector<uint32_t> vecActivated;
	double spread=0.0;
	//const auto numV = __F_graph.size();
	bool* activated = (bool *)calloc(__numV, sizeof(bool));
	uint32_t* visited = (uint32_t *)calloc(__numV, sizeof(uint32_t));  // visited is used to record in which simulation the node is visited. If not visited in the current simulation, then visit it and let the the node's value in visited be the number of current simulation.
	vector<double> vecThr(__numV);  // Threshold of LT
	vector<double> vecActivateWeight(__numV, 0.0);  // total weight of active neighbors in LT
	for (auto seedId : vecSeed) activated[seedId] = true;
	for (uint32_t i = 0; i < simulations; i++)
	{
		for (auto seed : vecSeed)
		{
			Que.push(seed);
		}

		// BFS traversal
		if (__model == "IC"||"ic")
		{
			while (!Que.empty())
			{
				nodeId = Que.front();
				Que.pop();
				for (auto& nbr : __F_graph[nodeId])
				{
					if (activated[get<0>(nbr)]) continue;
					if (dsfmt_gv_genrand_open_close() <= get<1>(nbr))
					{
						activated[get<0>(nbr)] = true;
						vecActivated.push_back(get<0>(nbr)); //Records which nodes are activated by the node.
						Que.push(get<0>(nbr));
					}
				}
			}
		}
		else if (__model == "LT"||"lt")
		{
			while (!Que.empty())
			{
				nodeId = Que.front();
				Que.pop();
				for (auto& nbr : __F_graph[nodeId])
				{
					if (activated[get<0>(nbr)]) continue;
					if (visited[get<0>(nbr)] < i + 1)
					{
						// First time visit this node
						visited[get<0>(nbr)] = i + 1;
						vecThr[get<0>(nbr)] = dsfmt_gv_genrand_open_close();
						vecActivateWeight[get<0>(nbr)] = 0.0;
					}
					vecActivateWeight[get<0>(nbr)] += get<1>(nbr);
					if (vecActivateWeight[get<0>(nbr)] >= vecThr[get<0>(nbr)])
					{
						// Activation weight is greater than threshold
						activated[get<0>(nbr)] = true;
						vecActivated.push_back(get<0>(nbr));
						Que.push(get<0>(nbr));
					}
				}
			}
		}
		uint32_t active = count(activated, activated+__numV, true);
		if(active >=__eta) { exceed=exceed+1; }
		spread += active;
		// if(i<101)
		// {
		// 	cout<<"the spread is :"<<active<<", ";
		// 	result_bk<<"the spread is :"<<active<<", ";
		// }
		
		// if(count(activated, activated+__numV, true) >=__eta) { exceed=exceed+1; }
		// spread += vecActivated.size();
		for (auto activatedNode : vecActivated) activated[activatedNode] = false;
		vecActivated.clear();
	}
	free(activated);
	free(visited);
	//auto simu_time=prob_eval_time.get_total_time();
	//cout<<"The time for simulation is "<<simu_time<<endl;
	return make_pair(1.0*exceed/simulations, 1.0*spread/simulations);
	//return make_pair(exceed, spread);
}

	

};//cls


using TAlg = Algorithm;
using PAlg = std::shared_ptr<TAlg>;
