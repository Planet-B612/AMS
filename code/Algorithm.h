#pragma once
#include "mRR-sets.h"
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
		TRRcollection RR, veriRR;
		Nodelist Seeds;
		vector<double> __cost;
		double __eta=1;
		double __eps_MC=0.1;
		double __delta_MC=0.1;
		double __tau=0.1;
		double __delta_1=0.1;
		double __lambda=__eps_MC;
		double __delta_2=0.1;
		uint8_t __mode=0; // 0: no prob estimation, 1: MC, 2: RR
		double __prob=0.1;
		double __epsilon_Inf=0.1;
		double __delta_Inf=0.01;
		double __veri_epsilon_Inf=0.05;
		Graph __F_graph;
		Graph __R_graph;
		CascadeModel __model=IC;
		string __result_dir;
		uint __E_PCG;
		double __RR_ratio=1.0;
		uint32_t __dist=1;
		uint16_t __THD=32;

		size_t __kappa=0;
		size_t __theta=0;
		vector<vector<bool>> __vecCover;  // record whether an FRset_real is covered by some seed
		vector<uint32_t> __vecCount;  // record how many FRset_real in a realization has been covered by the seeds
		vector<bool> __vecEta;	// record whether this realization has reached eta
		double __freq=0.0;

	public:

		Algorithm(Graph& graph) : RR(graph), veriRR(graph)//, FF(graph)    // The 'const' modifier if deleted
		{
			__numV = RR.get_nodes();
			__numE = RR.get_edges();
			//__cost.resize(cost.size());
			//__cost=cost;
			__R_graph=graph;
			//num_RRsets = num;
			//RR.build_n_RRsets(num);
		}
		~Algorithm()
		{ }

bool comp(tuple<uint32_t, uint32_t, double, uint32_t> a, tuple<uint32_t, uint32_t, double, uint32_t> b)
{
	return get<2>(a) > get<2>(b);
}; // sort in descending order.

void release_mem()
{
	RR.release_memory();
	veriRR.release_memory();
}

vector<uint32_t> max_ratio_lazy(const double Q)
{
	// The first element in tuple is the node ID, the second is the hyper-degree coverage, the third is the ratio, the fourth is the indicator whether the node has been updated in this round.
	vector<tuple<uint32_t, double, double, uint32_t>> ratio(__numV);
	double total_cost=0.0;
	//double maxDeg = 0.0;
	for (uint32_t i = __numV; i--;)
	{
		const uint32_t deg = RR._FRsets[i].size();  // The number of RR-sets covered by i.
		get<0>(ratio[i]) = i;
		get<1>(ratio[i])=deg;
		get<2>(ratio[i])=1.0*deg/__cost[i];
		get<3>(ratio[i])=0;
		//if (deg > maxDeg) maxDeg = deg;
	}
	make_max_heap(ratio);
	//sort(ratio.begin(),ratio.end(),comp);  // in descending order.
	// RRsets degMap(__numV); // degMap: map degree to the nodes with this degree
	// cout<<"maxDeg is: "<<maxDeg<<endl;
	// for (auto i = __numV; i--;)
	// {
	// 	if (coverage[i] == 0) continue;
	// 	degMap[coverage[i]].push_back(i);  // Record the nodes that have the same coverage. Map coverage to nodes.
	// }

	size_t total_deg = 0;  // Record the number of _RRsets covered by seed set.
	vector<bool> RR_Mark(num_RRsets, false);  // check if an edge is removed
	//FRsets veriFR=RR._FRsets;
	//veriRR.build_n_RRsets(num_RRsets);
	vector<bool> veriRR_state(num_RRsets, false);

	Seeds.clear();
	double expInf=0.0;
	double est_prob=0.0;
	// uint32_t last_seed=0;
	// bool equal=false;
	// uint32_t eq_cnt=0;
	// vector<uint32_t> eq_rec_node;
	// vector<double> eq_rec_cost;
	// vector<uint32_t> eq_rec_cnt;
	//cout<<"The seeds are ";
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
		else if (__mode==2)  // 2: RR
		{
			total_deg=total_deg+get<1>(ratio[0]);  // No prob estimation, only evaluate based on expected influence.
			expInf=1.0*__numV*total_deg/num_RRsets;
			est_prob=prob_est_RR(nodeId); 
			if(est_prob>=__prob+__lambda)  //lambda==__eps_MC
			//if(est_prob>=__prob)
			{
				cout<<"The expected inf is: "<<expInf<<", the estimated prob is: "<<est_prob<<endl;
				return Seeds;
			}
		}
		else
		{
			total_deg=total_deg+get<1>(ratio[0]);  // No prob estimation, only evaluate based on expected influence.
			//if(nodeId!=get<0>(ratio[0])) { cout<<"The node inserted is different from the node added to deg"; exit(1); }
			expInf=1.0*__numV*total_deg/num_RRsets;
			if(expInf>=Q)
			{
				cout<<"No estimation is needed. The expInf is "<<expInf<<endl;
				
				//veriFR.build_n_RRsets(500000);
				//uint32_t veri_num_RRsets=RR.get_RR_sets_size();
				
				//cout<<1.0*__numV*count(veriRR_state.begin(),veriRR_state.end(),true)/num_RRsets<<endl;
				return Seeds;
			}
		}

		//uint32_t cnt=0;
		// for(auto rr: veriRR._FRsets[nodeId])
		// {
		// 	if(veriRR_state[rr]==true) continue;
		// 	//cnt++;
		// 	veriRR_state[rr]=true;
		// }
		// if(cnt!=get<1>(ratio[0])) exit(1);
		// uint32_t veriDeg = count(veriRR_state.begin(),veriRR_state.end(),true);
		// if (((1.0*total_deg/veriDeg>=2)||(1.0*total_deg/veriDeg<=0.5))&&(i>2000))
		// {
		// 	cout<<"i="<<i<<". ";
		// 	cout<<veriDeg<<", "<<total_deg<<endl;
		// 	exit(1);
		// }
		
		for(auto rr:RR._FRsets[nodeId])
		{
			if(RR_Mark[rr]) continue;
			RR_Mark[rr]=true;
		}
		tuple<uint32_t, double, double, uint32_t> disable_node=make_tuple(nodeId, -1.0, -1.0, i);
		max_heap_replace_max_value(ratio, disable_node);
		// last_seed=nodeId;
		//vector<tuple<uint32_t, uint32_t, double, uint32_t>>::iterator it = find(ratio.begin(),ratio.end(),disable_node)-ratio.begin();
		// auto it = find(ratio.begin(),ratio.end(),disable_node);
		// int pos=distance(ratio.begin(), it);
		// cout<<pos<<", ";
	}
	cout<<"Probably an error exists in seed selection, this line should not appear."<<endl;
	cout<<"The expInf is "<<expInf<<endl;
	// for(int i=0;i<eq_rec_node.size();i++)
	// {
	// 	check<<eq_rec_node[i]<<"	"<<eq_rec_cost[i]<<"	"<<eq_rec_cnt[i]<<endl;
	// }
	return {};

}

vector<uint32_t> max_ratio_lazy_ECG_mine(const double Q)
{
	vector<tuple<uint32_t, double, double, uint32_t>> ratio(__numV);
	double total_cost=0.0;
	double RR_Deg_Thresh=1.0*(1+__epsilon_Inf)*Q*num_RRsets/__numV;
	//double veriRR_Deg_Thresh=1.0*(1+__veri_epsilon_Inf)*Q*__veri_num_RRsets/__numV;
	double veriRR_Deg_Thresh=1.0*Q*__veri_num_RRsets/__numV;
	for (uint32_t i = __numV; i--;)
	{
		const uint32_t deg = RR._FRsets[i].size();  // The number of RR-sets covered by i.
		get<0>(ratio[i]) = i;
		get<1>(ratio[i])=deg;
		get<2>(ratio[i])=1.0*deg/__cost[i];
		get<3>(ratio[i])=0;
	}
	make_max_heap(ratio);

	double total_deg = 0.0;  // Record the number of _RRsets covered by seed set.
	double veri_total_deg=0.0;
	double pre_total_deg=0.0;
	double eps_current=0.0;
	double w=log(__numV);
	double c_UB=0.0001;
	double expInf=0.0;
	double c_i=0.0001;
	fstream result_bk(__result_dir, ios::app);
	vector<bool> RR_Mark(num_RRsets, false);  // check if an edge is removed
	vector<bool> veriRR_Mark(__veri_num_RRsets, false);
	Seeds.clear();
	for(uint32_t i=0;i<__eta;i++) 
	{
		uint32_t nodeId;//=get<0>(ratio[0]);
		while(get<3>(ratio[0])!=i)
		{
			nodeId=get<0>(ratio[0]);
			//uint32_t nodeDeg=get<1>(ratio[0]);
			double nodeDeg=RR._FRsets[nodeId].size();
			for(uint32_t RRId: RR._FRsets[nodeId])
			{
				if(RR_Mark[RRId]==true) 	nodeDeg-=1.0;
			}
			nodeDeg=min(1.0*nodeDeg, RR_Deg_Thresh-1.0*total_deg);
			tuple<uint32_t, double, double, uint32_t> updated_node=make_tuple(nodeId, nodeDeg, 1.0*nodeDeg/__cost[nodeId],i);
			max_heap_replace_max_value(ratio, updated_node);
		}
		nodeId=get<0>(ratio[0]);
		Seeds.push_back( nodeId );
		c_i=__cost[nodeId];
		total_cost=total_cost+c_i;
		total_deg=total_deg+get<1>(ratio[0]);
		expInf=1.0*__numV*total_deg/num_RRsets;
		pre_total_deg=veri_total_deg;
		//pre_total_deg=total_deg;
		for(auto rr: veriRR._FRsets[nodeId])
		{
			if(veriRR_Mark[rr]) continue;
			else
			{
				veriRR_Mark[rr]=true;
				veri_total_deg+=1.0;
			}
		}
		veri_total_deg=min(1.0*veri_total_deg, veriRR_Deg_Thresh);
		// if(std::fmod(i+1,10)==0)
		// {
			// eps_current = 1.0*pre_total_deg/(power((sqrt(pre_total_deg+2.0*w/9.0))-sqrt(w/2.0),2) - sqrt(w/18.0))-1;
			// if(eps_current<__veri_epsilon_Inf) eps_current= __veri_epsilon_Inf;
		// }
		// c_UB=max( c_UB, 1.0*c_i*(__eta-(1.0+eps_current)*pre_total_deg)/((1.0+eps_current)*(veri_total_deg-pre_total_deg)+2.0*eps_current*pre_total_deg/c_i) );
		// if(i==0)
		// {
		// 	eps_current = 1.0*total_deg/(power((sqrt(total_deg+2.0*w/9.0))-sqrt(w/2.0),2) - w/18.0)-1.0;
			// if(eps_current<=0) 
			// { 
			// 	cout<<i<<", "<<total_deg<<", "<<eps_current<<endl; 
			// 	cout<<sqrt(total_deg+2.0*w/9.0)<<", "<<w<<", "<<(sqrt(total_deg+2.0*w/9.0))-sqrt(w/2.0)<<", "<<power((sqrt(total_deg+2.0*w/9.0))-sqrt(w/2.0),2)<<", "<<w/18.0<<", "<<(power((sqrt(total_deg+2.0*w/9.0))-sqrt(w/2.0),2) - w/18.0)<<", "<<1.0*total_deg/(power((sqrt(total_deg+2.0*w/9.0))-sqrt(w/2.0),2) - w/18.0)<<endl;
			// 	exit(1); 
			// }
			// ASSERT(eps_current>0);
		// 	if(eps_current<__epsilon_Inf) eps_current=1.0;
		// 	c_UB=max( c_UB, c_i*__eta/((1+eps_current)*total_deg));
		// }
		if((pre_total_deg>2.0*w/3.0)&&(__E_PCG==1)) // To make sure eps_current is positive.
		{
			eps_current = 1.0*pre_total_deg/(power((sqrt(pre_total_deg+2.0*w/9.0))-sqrt(w/2.0),2) - w/18.0)-1.0;
			//if(eps_current<=0) { cout<<i<<", "<<pre_total_deg<<", "<<eps_current<<endl;}
			//if(eps_current<__epsilon_Inf) eps_current= 1.0;
			c_UB=max( c_UB, 1.0*c_i*(__eta*__veri_num_RRsets/__numV-(1.0+eps_current)*pre_total_deg)/((1.0+eps_current)*(total_deg-pre_total_deg)+2.0*eps_current*pre_total_deg/c_i) );
		}

		if(veri_total_deg>=veriRR_Deg_Thresh-0.01)
		{
			cout<<"expInf = "<<expInf<<", total cost = "<<total_cost<<", c_UB = "<<c_UB<<". approx_ratio= "<<1.0*total_cost/c_UB<<endl;
			if(__E_PCG==1) result_bk<<"numV = "<<__numV<<", eta = "<<1.0*__eta/__numV<<", expInf = "<<expInf<<", veri_total_deg = "<<veri_total_deg<<", veriRR_Deg_Thresh = "<<veriRR_Deg_Thresh<<". approx_ratio = "<<1.0*total_cost/c_UB<<endl;
			result_bk.close();
			return Seeds;
		}
		for(auto rr:RR._FRsets[nodeId])
		{
			if(RR_Mark[rr]) continue;
			RR_Mark[rr]=true;
		}
		tuple<uint32_t, double, double, uint32_t> disable_node=make_tuple(nodeId, 0.0, -1.0, i);
		max_heap_replace_max_value(ratio, disable_node);
	}
	cout<<"Probably an error exists in seed selection, this line should not appear."<<endl;
	cout<<"veri_total_deg = "<<veri_total_deg<<", veriRR_Deg_Thresh = "<<veriRR_Deg_Thresh<<endl;
	return {};
}

vector<uint32_t> max_ratio_lazy_TKDE(const double Q)
{
	vector<tuple<uint32_t, double, double, uint32_t>> ratio(__numV);
	double total_cost=0.0;
	double Deg_Thresh=1.0*Q*num_RRsets/__numV;
	for (uint32_t i = __numV; i--;)
	{
		double deg = RR._FRsets[i].size();  // The number of RR-sets covered by i.
		deg=min(1.0*deg, Deg_Thresh);
		get<0>(ratio[i]) = i;
		get<1>(ratio[i])=deg;
		get<2>(ratio[i])=1.0*deg/__cost[i];
		get<3>(ratio[i])=0;
		//if (deg > maxDeg) maxDeg = deg;
	}
	make_max_heap(ratio);

	double total_deg = 0.0;  // Record the number of _RRsets covered by seed set.
	vector<bool> RR_Mark(num_RRsets, false);  // check if an edge is removed

	Seeds.clear();
	double expInf=0.0;
	for(uint32_t i=0;i<__eta;i++)  // The i-th seed we are selecting. After we selecting each seed, we should put it to the i-th place.
	{
		uint32_t nodeId=get<0>(ratio[0]);
		while(get<3>(ratio[0])!=i)  // since we will delete and insert the node, the j-th node would be different.
		{
			nodeId=get<0>(ratio[0]);
			double nodeDeg=1.0*RR._FRsets[nodeId].size();  // The num of RR-sets influenced by nodeId exclusively.
			for(uint32_t RRId: RR._FRsets[nodeId])
			{
				if(RR_Mark[RRId]==true) 	nodeDeg-=1.0;
			}
			nodeDeg=min(1.0*nodeDeg, Deg_Thresh-1.0*total_deg);
			tuple<uint32_t, double, double, uint32_t> updated_node=make_tuple(nodeId, nodeDeg, 1.0*nodeDeg/__cost[nodeId],i);
			max_heap_replace_max_value(ratio, updated_node);
		}
		nodeId=get<0>(ratio[0]);
		Seeds.push_back( nodeId );
		total_cost=total_cost+__cost[nodeId];

		total_deg=total_deg+get<1>(ratio[0]);  // No prob estimation, only evaluate based on expected influence.
		expInf=1.0*__numV*total_deg/num_RRsets;
		//if(i<3)	cout<<"The influence of the first "<<i<<" nodes is "<<expInf<<", The coverage is "<<total_deg<<endl;
		if(total_deg>=Deg_Thresh-0.01) // in case of calculation precision
		{
			cout<<"The expInf is "<<expInf<<". Q is "<<Q<<endl;
			return Seeds;
		}

		for(auto rr:RR._FRsets[nodeId])
		{
			if(RR_Mark[rr]) continue;
			RR_Mark[rr]=true;
		}
		tuple<uint32_t, double, double, uint32_t> disable_node=make_tuple(nodeId, -1.0, -1.0, i);
		max_heap_replace_max_value(ratio, disable_node);
		//last_seed=nodeId;
	}
	cout<<"Probably an error exists in seed selection, this line should not appear."<<endl;
	cout<<"The expInf is "<<expInf<<endl;
	return {};
}

vector<uint32_t> max_ratio_lazy_count_duplicate_seeds(const double Q)
{
	// The first element in tuple is the node ID, the second is the hyper-degree coverage, the third is the ratio, the fourth is the indicator whether the node has been updated in this round.
	vector<tuple<uint32_t, double, double, uint32_t>> ratio(__numV);
	double total_cost=0.0;
	//double maxDeg = 0.0;
	for (uint32_t i = __numV; i--;)
	{
		const uint32_t deg = RR._FRsets[i].size();  // The number of RR-sets covered by i.
		get<0>(ratio[i]) = i;
		get<1>(ratio[i])=deg;
		get<2>(ratio[i])=1.0*deg/__cost[i];
		get<3>(ratio[i])=0;
		//if (deg > maxDeg) maxDeg = deg;
	}
	make_max_heap(ratio);
	//sort(ratio.begin(),ratio.end(),comp);  // in descending order.
	// RRsets degMap(__numV); // degMap: map degree to the nodes with this degree
	// cout<<"maxDeg is: "<<maxDeg<<endl;
	// for (auto i = __numV; i--;)
	// {
	// 	if (coverage[i] == 0) continue;
	// 	degMap[coverage[i]].push_back(i);  // Record the nodes that have the same coverage. Map coverage to nodes.
	// }

	size_t total_deg = 0;  // Record the number of _RRsets covered by seed set.
	vector<bool> RR_Mark(num_RRsets, false);  // check if an edge is removed

	Seeds.clear();
	double expInf=0.0;
	double est_prob=0.0;
	uint32_t last_seed=0;
	ofstream check("~/graphInfo/check_info.txt", ios::app);
	bool equal=false;
	uint32_t eq_cnt=0;
	vector<uint32_t> eq_rec_node;
	vector<double> eq_rec_cost;
	vector<uint32_t> eq_rec_cnt;
	//cout<<"The seeds are ";
	for(uint32_t i=0;i<__numV;i++)  // The i-th seed we are selecting. After we selecting each seed, we should put it to the i-th place.
	{
		// if(std::fmod(i,15)==0) 
		// cout<<"selecting "<<i<<"-th seed"<<endl;
		uint32_t nodeId=get<0>(ratio[0]);
		while(get<3>(ratio[0])!=i)  // since we will delete and insert the node, the j-th node would be different.
		{
			nodeId=get<0>(ratio[0]);
			uint32_t nodeDeg=get<1>(ratio[i]);
			uint32_t num_covered=count(RR._FRsets[nodeId].begin(),RR._FRsets[nodeId].end(),true);
			nodeDeg=nodeDeg-num_covered;
			tuple<uint32_t, double, double, uint32_t> updated_node=make_tuple(nodeId, nodeDeg, 1.0*nodeDeg/__cost[nodeId],i);
			max_heap_replace_max_value(ratio, updated_node);
		}
		Seeds.push_back( nodeId );
		total_cost=total_cost+__cost[nodeId];
		if(last_seed==nodeId)
		{
			equal=true;
			eq_cnt++;
			if(i==__numV-1)
			{
				eq_rec_node.push_back(nodeId);
				eq_rec_cnt.push_back(eq_cnt);
				eq_rec_cost.push_back(__cost[nodeId]);
			}
		}
		else
		{
			if(eq_cnt!=0)
			{
				eq_rec_node.push_back(nodeId);
				eq_rec_cnt.push_back(eq_cnt);
				eq_rec_cost.push_back(__cost[nodeId]);
			}
			equal=false;
			eq_cnt=0;
		}
		last_seed=nodeId;

		if(__mode==1)	// 1: MC
		{
			total_deg=total_deg+get<1>(ratio[i]);  // No prob estimation, only evaluate based on expected influence.
			expInf=1.0*__numV*total_deg/num_RRsets;
			est_prob=prob_est_MC(Seeds);  
			if(est_prob>=__prob+__eps_MC)  //__eps_MC==lambda
			{
//				cout<<"reached prob"<<endl;
				//cout<<"The total cost is "<<total_cost<<endl;
				cout<<"The expected inf is: "<<expInf<<endl;
				return Seeds;
			}
		}
		else if (__mode==2)  // 2: RR
		{
			total_deg=total_deg+get<1>(ratio[i]);  // No prob estimation, only evaluate based on expected influence.
			expInf=1.0*__numV*total_deg/num_RRsets;
			est_prob=prob_est_RR(nodeId);  
			if(est_prob>=__prob+__lambda)  //lambda==__eps_MC
			{
				cout<<"reached prob"<<endl;
				cout<<"The expected inf is: "<<expInf<<endl;
				return Seeds;
			}
		}
		else
		{
			total_deg=total_deg+get<1>(ratio[i]);  // No prob estimation, only evaluate based on expected influence.
			expInf=1.0*__numV*total_deg/num_RRsets;
			if(expInf>=Q)
			{
				cout<<"No estimation is needed. The expInf is "<<expInf<<endl;
				return Seeds;
			}
		}

		for(auto rr:RR._FRsets[nodeId])
		{
			if(RR_Mark[rr]) continue;
			RR_Mark[rr]=true;
		}
		tuple<uint32_t, double, double, uint32_t> disable_node=make_tuple(nodeId, 0, 0.0, i);
		max_heap_replace_max_value(ratio, disable_node);
	}
	// cout<<"Probably an error exists in seed selection, this line should not appear."<<endl;
	// cout<<"The expInf is "<<expInf<<endl;
	for(int i=0;i<eq_rec_node.size();i++)
	{
		check<<eq_rec_node[i]<<"	"<<eq_rec_cost[i]<<"	"<<eq_rec_cnt[i]<<endl;
	}
	return {};

}

double max_cover_lazy(const double Q)
{
	//cout<<"The amplifier is :"<<amp<<endl;
	vector<double> coverage(__numV, 0);
	size_t maxDeg = 0;
	for (uint32_t i = __numV; i--;)
	{
		const double deg = RR._FRsets[i].size();  // The number of RR-sets covered by i.
		coverage[i] = deg;
		if (deg > maxDeg) maxDeg = deg;
	}
	RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
	//cout<<"maxDeg is: "<<maxDeg<<endl;
	for (auto i = __numV; i--;)
	{
		if (coverage[i] == 0) continue;
		degMap[coverage[i]].push_back(i);  // Record the nodes that have the same coverage.
	}

	size_t total_deg = 0;  // Record the number of _RRsets covered by seed set.
	vector<bool> edgeMark(num_RRsets, false);  // check if an edge is removed

	Seeds.clear();
	double expInf=0.0;
	//cout<<"The seeds are ";
	for (auto deg = maxDeg; deg > 0; deg--) // In each round, try every node that may cover current deg, i.e., nodes in degMap[deg].
	{
		auto& vecNode = degMap[deg];  // vecNode is the set of nodes that covers deg _RRsets. With &, the operation on vecNode would affect degMap as well.
		for (auto idx = vecNode.size(); idx--;)
		{
			auto argmaxIdx = vecNode[idx];
			const auto currDeg = coverage[argmaxIdx];  // During seed selection, the coverage of nodes in vecNode would decrease.
			if (deg > currDeg)  // Then, the current node argmaxIdx should not be in degMap[deg], but in degMap[currDeg].
			{
				degMap[currDeg].push_back(argmaxIdx);
				continue;
			}
			Seeds.push_back(argmaxIdx);
			// cout<<argmaxIdx<<"; ";
			total_deg += currDeg;
			//cout<<endl<<__numV<<", "<<num_RRsets<<endl;
			expInf = 1.0 * total_deg * __numV / num_RRsets;
			//expInf = 1.0 * total_deg * Q / num_RRsets;
			// if (expInf >= Q)
			// {
			// 	// Top-k influential nodes constructed
			// 	//expInf = 1.0 * total_deg * __numV / num_RRsets;
			// 	//expInf = 1.0 * total_deg * Q / num_RRsets;
			// 	cout<<"The # of seeds is "<<Seeds.size()<<", the total_deg is "<<total_deg;
			// 	std::cout << "  >>>[greedy-lazy] expInf: " << expInf<<endl;
			// 	return expInf;
			// }
			double prob=prob_est_MC(Seeds);  // MC by default
			if(prob>=__prob)
			{
				return Seeds.size();
			}
			coverage[argmaxIdx] = 0;
			for (auto edgeIdx : RR._FRsets[argmaxIdx])
			{
				if (edgeMark[edgeIdx]) continue;
				edgeMark[edgeIdx] = true;
				for (auto nodeIdx : RR._RRsets[edgeIdx])
				{
					if (coverage[nodeIdx] == 0) continue; // This node is seed, skip
					coverage[nodeIdx]--;  // Refresh the coverage of nodes.
				}
			}
		}
		degMap.pop_back(); // why?
	}
	cout<<"???"<<endl;
	return 1.0 * __numV; // All RR sets are covered.
}

double max_cover_topk(const double Q)
{
	FRset coverage(__numV, 0);
	// Obtain maxDeg
	size_t maxDeg = 0;
	for (auto i = __numV; i--;)
	{
		const auto deg = RR._FRsets[i].size();
		coverage[i] = deg;
		if (deg > maxDeg) maxDeg = deg;
	}
	// Construct degMap
	RRsets degMap(maxDeg + 1); // degMap: map degree to the nodes with this degree
	for (auto i = __numV; i--;)
	{
		//if (coverage[i] == 0) continue;
		degMap[coverage[i]].push_back(i);
	}
	Nodelist sortedNode(__numV); // sortedNode: record the sorted nodes in ascending order of degree
	Nodelist nodePosition(__numV); // nodePosition: record the position of each node in the sortedNode
	Nodelist degreePosition(maxDeg + 2); // degreePosition: the start position of each degree in sortedNode
	uint32_t idxSort = 0;
	size_t idxDegree = 0;
	for (auto& nodes : degMap)
	{	// degreePosition is not initialized??
		degreePosition[idxDegree + 1] = degreePosition[idxDegree] + (uint32_t)nodes.size();
		idxDegree++;
		for (auto& node : nodes) // Initialize nodePosition, sortedNode as there position in DegMap
		{
			nodePosition[node] = idxSort;
			sortedNode[idxSort++] = node;
		}
	}
	// check if an edge is removed
	std::vector<bool> edgeMark(num_RRsets, false);
	// record the total of top-k marginal gains, calculate sumtopk.
	// size_t sumTopk = 0;
	// for (auto deg = maxDeg + 1; deg--;)
	// {
	// 	if (degreePosition[deg] <= __numV - targetSize)
	// 	{
	// 		sumTopk += deg * (degreePosition[deg + 1] - (__numV - targetSize));
	// 		break;
	// 	}
	// 	sumTopk += deg * (degreePosition[deg + 1] - degreePosition[deg]);
	// }
	// __boundMin = 1.0 * sumTopk;
	Seeds.clear();
	size_t total_deg = 0;
	/*
	* sortedNode: position -> node
	* nodePosition: node -> position
	* degreePosition: degree -> position (start position of this degree)
	* coverage: node -> degree
	* e.g., swap the position of a node with the start position of its degree
	* swap(sortedNode[nodePosition[node]], sortedNode[degreePosition[coverage[node]]])
	*/
	double expInf = 0.0;
	while(expInf<Q)
	{
		const auto seed = sortedNode.back(); // The last node in sortedNode
		sortedNode.pop_back();
		//const auto newNumV = sortedNode.size();
		//sumTopk += coverage[sortedNode[newNumV - targetSize]] - coverage[seed];
		total_deg += coverage[seed];
		Seeds.push_back(seed);
		coverage[seed] = 0;
		for (auto edgeIdx : RR._FRsets[seed])
		{
			if (edgeMark[edgeIdx]) continue;
			edgeMark[edgeIdx] = true;
			for (auto nodeIdx : RR._RRsets[edgeIdx])
			{
				if (coverage[nodeIdx] == 0) continue; // This node is seed, skip
				const auto currPos = nodePosition[nodeIdx]; // The current position
				const auto currDeg = coverage[nodeIdx]; // The current degree
				const auto startPos = degreePosition[currDeg]; // The start position of this degree
				const auto startNode = sortedNode[startPos]; // The node with the start position
				// Swap this node to the start position with the same degree, and update their positions in nodePosition
				std::swap(sortedNode[currPos], sortedNode[startPos]);
				nodePosition[nodeIdx] = startPos;
				nodePosition[startNode] = currPos;
				// Increase the start position of this degree by 1, and decrease the degree of this node by 1
				degreePosition[currDeg]++;
				coverage[nodeIdx]--;
				// If the start position of this degree is in top-k, reduce topk by 1
				//if (startPos >= newNumV - targetSize) sumTopk--;
			}
		}
		expInf = 1.0 * total_deg * Q / num_RRsets;
	}
	// std::cout << "  >>>[greedy-topk] influence: " << expInf << '\n';
	return expInf;
}

// double max_cover(const double Q)
// {
// 	//if (Q >= 500000) return max_cover_topk(Q);
// 	return max_cover_lazy(Q);
// }

void set_cascade_model(const CascadeModel model)
{
	__model=model;
	RR.set_cascade_model(model);
	//FF.set_cascade_model(model);
}

// double prob_est(const vector<uint32_t> vecSeed, double Q)
// {
// 	if(__mode==2) // RR
// 	{
// 		return prob_est_RR(vecSeed);
// 	}
// 	if(__mode==1) // MC
// 	{
// 		return prob_est_MC(vecSeed);
// 	}
// 	else
// 	{
// 		cout<<"This mode value should not invoke prob_est()"<<endl;
// 		exit(1);
// 	}
// }

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

bool TEST(Nodelist A, double kappa, double Gamma, double beta, uint32_t L)
{
	uint32_t n=__numV;
	uint32_t l=ceil(2*(1+kappa)*Gamma/((2+kappa)*n) + 8*(3+2*kappa)*(1+kappa)*log(2/beta)/(3*kappa*kappa));
	uint32_t M=floor((2+kappa)*n*l/(2*(1+kappa)*Gamma));
	//RRsets U;
	bool Pass=false;
	uint32_t Z=0;
	if(L<=M)
	{
		//uint32_t a=L+RR.get_RR_sets_size();
		RR.build_n_RRsets(L+RR._RRsets.size());
		//num_RRsets=RR.get_RR_sets_size();
		return Pass;
	}
	for(uint32_t i=0;i<M;i++)
	{
		//RRcollection Ui(R_graph);
		RR.build_n_RRsets(RR._RRsets.size()+1);
		//num_RRsets++;
		//RR._RRsets.push_back(Ui._RRsets[0]);
		int count=0;
		for(auto node: RR._RRsets[RR._RRsets.size()-1])  // Considering the latest RR set.
		{
			for(auto ai:A)
			{
				if(node==ai)
				count++;
			}
		}
		Z+=min(count,1);
		if(Z==l)
		{
			Pass=true;
			//cout<<"Pass becomes true. The num of RR sets is "<<RR.get_RR_sets_size()<<endl;
			return Pass;
		}
	}
	return Pass;
}

vector<uint32_t> tegc(const double delta, const double alpha, const double sigma, const double gamma)
{
	//uint32_t numV=R_graph.size();
	//double D=power(exp(1)*numV/floor((1-alpha)*eta),floor((1-alpha)*eta) );
	//double D=power(numV, )
	//double D=floor((1-alpha)*__eta);
	double D=__numV;
	double W1[3]={(1-alpha)*__eta, gamma/(1-alpha), D};
	double W2[3]={__eta, sigma, delta/6};
	double ut=2*__numV*(3+W1[1])* ( log(6/delta)+8*log(D) ) /(3*W1[1]*W1[1]*W1[0]);
	double lt=2*__numV*log(1/W2[2])/(W2[1]*W2[1]*W2[0]);
	auto T=ceil(max(ut,lt));
	//cout<<"The max-num of RR-sets needed is "<<T<<endl;
	double theta=delta/3;
	//RRsets R;
	Nodelist vecSeed;
	while(RR._RRsets.size()<=T)
	{
		uint32_t R0=min( T, ceil(2*__numV*log(3/theta)/(sigma*sigma*__eta)) );
		if(RR._RRsets.size()<R0)
		{
			RR.build_n_RRsets(R0); // The function will buuild RR-sets to R0.
			//num_RRsets+=R0-RR._RRsets.size();
		}
		num_RRsets=RR.get_RR_sets_size();
		max_ratio_lazy_TKDE((1-alpha+gamma)*__eta);
		vecSeed=Seeds;
		if(RR._RRsets.size()==T)
		{
			cout<<"The num of RR-sets in TEGC is "<<RR._RRsets.size()<<endl;
			return vecSeed;
		}
		bool Pass=TEST(vecSeed, gamma/(2*(1-alpha)),(1-alpha)*__eta, 2*theta/3, T-RR._RRsets.size());
		if(Pass==true)
		{
			cout<<"The num of RR-sets in TEGC is "<<RR._RRsets.size()<<endl;
			return vecSeed;
		}
		//R=R.push_back(R) has been done in TEST().
		theta=theta/2;
	}
	cout<<"The num of RR-sets in TEGC is "<<RR._RRsets.size()<<endl;
	return vecSeed;
}

vector<uint32_t> bcgc(const double delta, const double alpha, const double sigma, const double gamma)
{
	//uint32_t numV=R_graph.size();
	//uint32_t mu=power(numV, 8);
	double Lambda=(1-alpha+gamma)*__eta;
	double W1[3]={(1-alpha)*__eta, gamma/(1-alpha), delta/(2)};
	double W2[3]={__eta, sigma, delta/2};
	double ut=2*__numV*(3+W1[1])*(  log(1/W1[2]) +  8*log(__numV) )/(3*W1[1]*W1[1]*W1[0]);  //floor((1-alpha)*__eta)*log(exp(1)*__numV/floor((1-alpha)*__numV))
	double lt=2*__numV*log(1/W2[2])/(W2[1]*W2[1]*W2[0]);
	auto T=ceil(max(ut,lt));
	//double xi=0.002*__eta;
	//T=(2.0*(1+1.0*__epsilon_Inf/3)*log(1.0/__delta_Inf)*__numV/(__epsilon_Inf*__epsilon_Inf*xi));
	cout<<"The num of RR-sets needed is "<<T<<endl;
	// Timer bcgc_RR_MC("BCGC-RRset");
	RR.build_n_RRsets(T);
	// auto bcgc_RR_time=bcgc_RR_MC.get_total_time();
	// cout<<"The time for BCGC RR-sets is "<<bcgc_RR_time<<endl;
	num_RRsets=RR.get_RR_sets_size();  // num_RRsets must be updated, since it is used in max_cover_lazy() which does not invoke get_RR_sets_size().
	//max_ratio_lazy_Trun(Lambda);
	max_ratio_lazy_TKDE(Lambda);
	//max_ratio_lazy((1+alpha)*__eta);
	//max_cover_lazy(Lambda);
	return Seeds;
}

vector<uint32_t> mine_ECG(const double delta, const double sigma, const double gamma)
{
	double xi=0.00017*__eta;
	//cout<<xi<<endl;
	num_RRsets=(2.0*(1+1.0*__epsilon_Inf/3)*log(1.0/__delta_Inf)*__numV/(__epsilon_Inf*__epsilon_Inf*xi));
	//num_RRsets=65800000;
	if(__E_PCG==1) num_RRsets=ceil(num_RRsets*__RR_ratio);
	cout<<"The num of RR-sets needed is "<<num_RRsets<<endl;
	// Timer mine_RR_MC("mine-RRset");
	RR.build_n_RRsets(num_RRsets);
	__veri_num_RRsets=2.0*(1+1.0*__veri_epsilon_Inf/3)*log(1.0/__delta_Inf)*__numV/(__veri_epsilon_Inf*__veri_epsilon_Inf*__eta);
	if(__E_PCG==1) __veri_num_RRsets=ceil(num_RRsets*__RR_ratio);
	//__veri_num_RRsets=num_RRsets;
	veriRR.build_n_RRsets(__veri_num_RRsets);
	//cout<<"Current memory usage after generating RR-sets is "<<getProcMemory()<<endl;
	// auto mine_RR_time=mine_RR_MC.get_total_time();
	// cout<<"The time for Mine RR-sets is "<<mine_RR_time<<endl;
	//num_RRsets=RR.get_RR_sets_size();  // num_RRsets must be updated, since it is used in max_cover_lazy() which does not invoke get_RR_sets_size().
	//cout<<"Current memory usage after generating MFFcol is "<<getProcMemory()<<endl;
	max_ratio_lazy_ECG_mine((1.0+__veri_epsilon_Inf)*__eta);
	//max_ratio_lazy_TKDE((1-0.04/3)*__eta);
	//max_cover_lazy(Lambda);
	return Seeds;
}

vector<uint32_t> mine_PCG(const double delta, const double sigma, const double gamma)
{
	// double Lambda=(1+gamma)*__eta;
	// double W1[3]={__eta, gamma, delta/2};
	// double W2[3]={__eta, sigma, delta/2};
	//double ut=2*__numV*(3+W1[1])* (log(1/W1[2])+8*log(__numV)) /(3*W1[1]*W1[1]*W1[0]);
	//double ut=2*__numV*(3+W1[1])*(  log(1/W1[2]) + floor((1)*__eta)*log(exp(1)*numV/floor((1-alpha)*numV))  )/(3*W1[1]*W1[1]*W1[0]);
	//double ut=2*__numV*(3+W1[1])*(  log(1/W1[2]) + floor(__eta)*log(exp(1)*__numV/floor(__numV))  )/(3*W1[1]*W1[1]*W1[0]);
	//double lt=2*__numV*log(1/W2[2])/(W2[1]*W2[1]*W2[0]);
	//uint32_t T=2*ceil(max(ut,lt));
	//T=ceil(min(ut,lt));  // for fast test
	//Timer mine_RR_MC("mine-RRset");
	double xi=0.00017*__eta; // Parameter setting particularlly for friendster
	//double xi=0.002*__eta;
	//cout<<xi<<endl;
	num_RRsets=(2.0*(1+1.0*__epsilon_Inf/3)*log(1.0/__delta_Inf)*__numV/(__epsilon_Inf*__epsilon_Inf*xi));
	RR.build_n_RRsets(num_RRsets);
	cout<<"The num of RR-sets needed is "<<num_RRsets<<endl;
	cout<<"Current memory usage after generating RR-sets is "<<getProcMemory()<<endl;
	//auto mine_RR_time=mine_RR_MC.get_total_time();
	//cout<<"The time for Mine RR-sets is "<<mine_RR_time<<endl;
	//num_RRsets=RR.get_RR_sets_size();  // num_RRsets must be updated, since it is used in max_cover_lazy() which does not invoke get_RR_sets_size().
	__kappa=ceil( log(1/__delta_2)/(2*__lambda*__lambda) );  // # realizations
	__theta=ceil( log(1/__delta_1)/(2*__tau*__tau) );  // The num of roots
	__vecCover.resize(__kappa, vector<bool>(__theta, false));
	__vecEta.resize(__kappa, false);
	__vecCount.resize(__kappa, 0.0);
	//cout<<"The num of realization is "<<__kappa<<", The num of roots in each realization is "<<__theta<<endl;
	Timer MINE_MFFcol("MINE_MFFcol");
	RR.generate_MFFcol_20231029(__kappa, __theta);
	auto MINE_MFFcol_time=MINE_MFFcol.get_total_time();
	cout<<"The time for MINE MFFcol generation is "<<MINE_MFFcol_time<<endl;
	//Timer MINE_Greedy("MINE_Greedy");
	cout<<"Current memory usage after generating MFFcol is "<<getProcMemory()<<endl;
	max_ratio_lazy(__eta);
	//auto MINE_Greedy_MC_time=MINE_Greedy.get_total_time();
	//cout<<"The time for MINE greedy is "<<MINE_Greedy_MC_time<<endl;
	//max_cover_lazy(Lambda);
	return Seeds;
}

vector<uint32_t> kdd(const double delta, const double sigma, const double gamma)  // alpha should be 0 in kdd. The setting of T assures that the deviation is at most gamma (sigma), according to Lemma 1 in TKDE, by simple Chernoff bound.
{
	//double Lambda=(1+gamma)*__eta;
	double Lambda=__eta;
	double W1[3]={__eta, gamma, delta/2};
	double W2[3]={__eta, sigma, delta/2};
	//double ut=2*__numV*(3+W1[1])*(  log(1/W1[2]) + floor(__eta)*log(exp(1)*__numV/floor(__numV))  )/(3*W1[1]*W1[1]*W1[0]);
	double ut=2*__numV*(3+W1[1])* (log(1/W1[2])+8*log(__numV)) /(3*W1[1]*W1[1]*W1[0]);
	double lt=2*__numV*log(1/W2[2])/(W2[1]*W2[1]*W2[0]);
	auto T=ceil(max(ut,lt));  		// The determination of RR-set size could also be the same as TEGC.
	cout<<"The num of RR-sets needed is "<<T<<endl;
	//T=min(ut,lt); // for quick test
	Timer kdd_RR_MC("KDD-RRset");
	RR.build_n_RRsets(T);
	auto kdd_RR_time=kdd_RR_MC.get_total_time();
	cout<<"The time for KDD RR-sets is "<<kdd_RR_time<<endl;
	num_RRsets=RR.get_RR_sets_size();
	//cout<<"RR-sets for KDD have been generated."<<endl;
	//cout<<"The num of simulations in each MC is "<<ceil( log(2/__delta_MC)/(2*__eps_MC*__eps_MC) )<<endl;
	Timer kdd_greedy_MC("KDD");
	max_ratio_lazy(Lambda);
	auto kdd_time=kdd_greedy_MC.get_total_time();
	//cout<<"The time for kdd to run greedy+MC is "<<kdd_time<<endl;
	return Seeds;
}

/// @brief Estimate the probability with RR-sets
/// @return the estimated probability
/// Different from prob_est_MC, the delta here is an absolute value instead of implying 1/n^\delta
/// Note that we delete the identifier "const" in front of Graph& graph, since we want to revise the edges. Will the changes be automatically applied to the original vec.rev.graph file?
/// The construction of RR-sets seems to be unnecessary, since we only want to know the spacific number of RR-sets covered by the seeds.
double prob_est_RR(const uint32_t seed)
{
	//Timer RR_time("RR_time");
	// size_t kappa=ceil( 0.3*log(1/__delta_2)/(2*__lambda*__lambda) );  // # realizations
	// size_t theta=2*ceil( log(1/__delta_1)/(2*__tau*__tau) );  // The num of roots
	//RR.generate_MFFcol(kappa, theta);
	// double freq=0.0;
	//for(FRsets FRsets_real: RR._MFFcol)
	for(uint32_t i=0;i<__kappa;i++)
	{
		//vector<bool> vecCover(__theta,false);
		//double freq=0.0;
		//uint32_t cnt=0;  // This is wrong. In each time of estimation, the covered RR-sets of each realization is set to be 0, which however should be accumulative.
		//for(uint32_t FRset_real:FRsets_real[seed])
		if(__vecEta[i])	continue;  // If this realization has already reached eta, there is no need to update it again.
		for(uint32_t FRset_real:RR._MFFcol[i][seed])
		{
			if(__vecCover[i][FRset_real]) { continue; }
			else
			{
				__vecCover[i][FRset_real]=true;
				__vecCount[i]+=1;
			}
		}
		//if(1.0*__vecCount[i]*__numV/__theta - 1.0*__tau*__numV >= __eta)	
		if(1.0*__vecCount[i]*__numV/__theta >= __eta)	
		{
			__freq+=1;
			__vecEta[i]=true;
		}
		// if(1.0*count(vecCover.begin(),vecCover.end(),true)*__numV/__theta - 1.0*__tau*__numV >= __eta)
		// {	//cout<<freq<<" ";	
		// freq+=1.0;  }
	}
	//const auto total_time=RR_time.get_total_time();
	return 1.0*__freq/__kappa;
}

/// Evaluate influence spread using MC simulation.
// This function should be carried out on forward graph!!
double prob_est_MC(const vector<uint32_t>& vecSeed)
{
	ASSERTT(!__F_graph.empty(), "The input forward graph __F_graph for prob_est_MC() is empty.");
	Timer MC_time("MC_time");
	// cout << "  >>>Evaluating probability...\n";
	uint32_t nodeId, currProgress = 0;
	double prob=0.0;
	queue<uint32_t> Que;
	vector<uint32_t> vecActivated;
	//const auto numV = graph.size();
	uint32_t simulations=ceil( log(2/__delta_MC)/(2*__eps_MC*__eps_MC) );
	ASSERTT(simulations>0, "The value of MC simulation times is invalid!");
	//cout<<"We need to do "<<simulations<<" times of MC"<<endl;
	//double spread = (double)vecSeed.size();
	bool* activated = (bool *)calloc(__numV, sizeof(bool));
	uint32_t* visited = (uint32_t *)calloc(__numV, sizeof(uint32_t));  // visited is used to record in which simulation the node is visited. If not visited in the current simulation, then visit it and let the the node's value in visited be the number of current simulation.
	vector<double> vecThr(__numV);  // Threshold of LT
	vector<double> vecActivateWeight(__numV, 0.0);  // total weight of active neighbors in LT
	for (auto seedId : vecSeed) activated[seedId] = true;
	for (uint32_t i = 0; i < simulations; i++)
	{
		// if (i * 100 >= simulations * currProgress)
		// {
		// 	const auto evalTime = MC_time.get_operation_time();
		// 	//if (evalTime > 100)
		// 		//std::cout << "\tMC-Progress at: " << currProgress << "%, " << "time used: " << evalTime << std::endl;
		// 	std::cout << "\tMC-Progress at: " << currProgress << "%, " << model<< " time used: " << evalTime<<std::endl;
		// 	currProgress += 20;
		// }
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
		uint32_t num_activated=count(activated, activated+__numV, true);
		if(count(activated, activated+__numV, true) >=__eta) { prob=prob+1; }
		//cout<<"The num of influenced nodes is "<<num_activated<<endl;
		//spread += (double)vecActivated.size() / simulations;
		for (auto activatedNode : vecActivated) activated[activatedNode] = false;
		vecActivated.clear();
	}
	free(activated);
	free(visited);
	//cout << "  >>>MC-Influence spread: " << spread << ", time used (sec): " << EvalTimer.get_total_time() << endl;
	// cout<<prob<<" times of spread meets the standard"<<endl;
	//cout<<"Current probability is "<<1.0*prob/simulations<<endl;
	return 1.0*prob/simulations;
}

pair<double, double> multi_THD_eval(const vector<uint32_t> & vecSeed, uint32_t simulations)
{
	Timer prob_eval_time("prob_eval_time");
	double prob=0.0;
	double spread=0.0;
	pair<double,double> res;
	vector<future<pair<double, double>>> futures(__THD);
	uint32_t avg_load = simulations/__THD +1;
	vector<uint32_t> actual_load(__THD, 0);
	for (uint16_t id = 0; id < __THD; id++)
	{
		uint32_t start = id*avg_load;
		uint32_t end = start+avg_load;
		//cout<<start<<", "<<end<<endl;
		end=min(end, simulations);
		actual_load[id]=end-start;
		//cout<<actual_load[id]<<endl;
		futures[id]=std::async(std::launch::async, &Algorithm::actual_prob_eval, this,
                                  vecSeed, actual_load[id]);
	}
	std::for_each(futures.begin(), futures.end(), std::mem_fn(&std::future<pair<double, double>>::wait));

	for (uint16_t i = 0; i < __THD; i++)
	{
		res=futures[i].get();
		prob+=1.0*res.first*actual_load[i];
		spread+=1.0*res.second*actual_load[i];
	}
	prob=1.0*prob/simulations;
	spread=1.0*spread/simulations;
	auto simu_time=prob_eval_time.get_total_time();
	cout<<"The time for simulation is "<<simu_time<<endl;
	return make_pair(prob, spread);
}

// pair<double, double> single_THD_eval(uint16_t tid, const vector<uint32_t> & vecSeed, uint32_t simulations)
// {

// }

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
