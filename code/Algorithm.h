#pragma once
#include "RRsets.h"
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
		size_t __numE=0;
		size_t num_RRsets = 0;
		RRsets RR;
		Nodelist Seeds;
		vector<double> __cost;
		double __eta=1;
		double __delta=0.1;
		double __epsilon_Inf=0.1;
		double __delta_Inf=0.01;
		Graph __F_graph;
		Graph __R_graph;
		CascadeModel __model;
		string __result_dir;
		vector<vector<bool>> __vecCover;  // record whether an FRset_real is covered by some see

	public:

		Algorithm(Graph& graph) : RR(graph)
		{
			__numV = RR.get_nodes();
			__numE = RR.get_edges();
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

		if(__mode==1)	// 1: MC
		{
			total_deg=total_deg+get<1>(ratio[0]);  // No prob estimation, only evaluate based on expected influence.
			expInf=1.0*__numV*total_deg/num_RRsets;
			if(std::fmod(i, __dist)==0) est_prob=prob_est_MC(Seeds);  
			if(est_prob>=__prob+__eps_MC)  //__eps_MC==lambda
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

void set_cascade_model(const CascadeModel model)
{
	__model=model;
	RR.set_cascade_model(model);
	//FF.set_cascade_model(model);
}

void set_parameters(const CascadeModel model, vector<double> cost, uint8_t mode, double prob, Graph F_graph, double eta, double eps_MC, double delta_MC, double tau, double delta_1, double lambda, double delta_2, double eps_Inf, double delta_Inf, double veri_epsilon_Inf, string result_dir, uint E_PCG, double RR_ratio, uint32_t dist, uint16_t THD)
{
	RR.set_model_mode(model, mode);
	__cost=cost;
	__mode=mode;
	__prob=prob;
	__eta=eta;
	__eps_MC=eps_MC;
	delta_MC=delta_MC;
	__tau=tau;
	__delta_1=delta_1;
	__lambda=lambda;
	__delta_2=delta_2;
	__F_graph=F_graph;
	__epsilon_Inf=eps_Inf;
	__delta_Inf=delta_Inf;
	__veri_epsilon_Inf=veri_epsilon_Inf;
	__result_dir=result_dir;
	__E_PCG=E_PCG;
	__RR_ratio=RR_ratio;
	__dist=dist;
	__THD=THD;
}

void set_F_graph(const Graph F_graph)
{
	__F_graph=F_graph;
}

double effic_inf_valid_algo(const double delta = 1e-3, const double eps = 0.01)
{
	return RR.effic_inf_valid_algo(Seeds, delta, eps);
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
		if (__model == IC)
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
		else if (__model == LT)
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

	/// Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
	double effic_inf_valid_algo(const vector<uint32_t>& vecSeed, const double delta = 1e-3, const double eps = 0.01)
	{
		const double c = 2.0 * (exp(1.0) - 2.0);
		const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
		//const double LambdaL = 2100000;
		//cout<<"The num of RR-sets is "<<LambdaL<<endl;
		size_t numHyperEdge = 0;
		size_t numCoverd = 0;
		std::vector<bool> vecBoolSeed(__numV);
		for (auto seed : vecSeed) vecBoolSeed[seed] = true;
		std::vector<bool> __vecVisitBool(__numV, false);
		std::vector<uint32_t> __vecVisitNode(__numV);

		while (numCoverd < LambdaL)
		{
			numHyperEdge++;
			size_t numVisitNode = 0, currIdx = 0;
			const auto uStart = dsfmt_gv_genrand_uint32_range(__numV);
			if (vecBoolSeed[uStart])
			{
				// Stop, this sample is covered
				numCoverd++;
				continue;
			}
			__vecVisitNode[numVisitNode++] = uStart;
			__vecVisitBool[uStart] = true;
			while (currIdx < numVisitNode)
			{
				const auto expand = __vecVisitNode[currIdx++];
				if (__model == IC)
				{
					for (auto& nbr : __R_graph[expand])
					{
						const auto nbrId = get<0>(nbr);
						if (__vecVisitBool[nbrId])
							continue;
						const auto randDouble = dsfmt_gv_genrand_open_close();
						if (randDouble > get<1>(nbr))
							continue;
						if (vecBoolSeed[nbrId])
						{
							// Stop, this sample is covered
							numCoverd++;
							goto postProcess;
						}
						__vecVisitNode[numVisitNode++] = nbrId;
						__vecVisitBool[nbrId] = true;
					}
				}
			}
		postProcess:
			for (auto i = 0; i < numVisitNode; i++)
				__vecVisitBool[__vecVisitNode[i]] = false;
		}
		return 1.0 * numCoverd * __numV / numHyperEdge;
	}

};//cls


using TAlg = Algorithm;
using PAlg = std::shared_ptr<TAlg>;
