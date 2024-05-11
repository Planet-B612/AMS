#include "../dSFMT/dSFMT.h"
#include "graph.h"
#include <iostream>
#include <vector>
#include "CommonStruc.h"
#include <cstring>
#include "Timer.h"
#include "Memory.h"
#include "MemoryUsage.h"

using namespace std;

int main(int argn, char **argv)
{
    float eta_0=0.1;
    uint eta_start=0;
    uint eta_end=1;
    string model="IC";
    bool Rnd_cost=false;
    int simRnd=500;
    uint format_graph=0;  // 0: do not format graph, 1: form the forward graph, 2: form the reverse graph.
    vector<string> dataset={"facebook", "livejournal", "pokec", "dblp", "friendster"};
    vector<string> alg_arr={"BCGC", "CLEAR", "SCORE", "ASM", "MINE"};
    vector<int> algs={0};
	vector<int> data={0};
    for (int i = 0; i < argn; i++)
    {
        if (argv[i] == string("-dataset"))
			data = {stoi(argv[i + 1])};
		if (argv[i] == string("-alg"))
			algs = {stoi(argv[i + 1])};
        if (argv[i] == string("-model"))
			model = argv[i + 1]; 
        if (argv[i] == string("-format_graph"))
			format_graph = stoi(argv[i + 1]);
        if (argv[i] == string("-simRnd"))
			simRnd = stoi(argv[i + 1]);
        if (argv[i] == string("-Rnd_cost"))
			Rnd_cost = stoi(argv[i + 1]);
    }
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
    string parameter_str="# The parameters are: Rnd_cost="+to_string(Rnd_cost)+", model="+model;
    cout<<parameter_str<<endl;
    string result_dir="./results/backup.txt";
    std::fstream result_bk(result_dir, ios::app);
	assert(!result_bk.fail());
	result_bk<<parameter_str<<endl;
    for (auto k:data)
    {
        string dataset_dir = "../graphInfo/"+dataset[k];
        if(format_graph!=0)
        {
            GraphBase::format_graph(dataset_dir, format_graph-1);
            return 0;
        }
        Graph R_graph = GraphBase::load_graph(dataset_dir, 1);
	    uint32_t numV=R_graph.size();
        vector<double> cost(numV,1.0);
        double total_cost=0.0;
        vector<double> reg_prob(numV,0.0);
        string cost_file;
		if(Rnd_cost)
		{
			cost_file="../graphInfo/"+dataset[k]+"_cost.txt";
		}
		else
		{
			cost_file="../graphInfo/"+dataset[k]+"_cost_001DEG.txt";
		}
		std::ifstream inFile;
		inFile.open(cost_file);
		if(!inFile)
		{
			cout<<"cannot open the cost file at "<<cost_file<<endl;
			exit(1);
		}
    	inFile.seekg(0, std::ios_base::beg);
		for(size_t i=0;i<numV;i++)
		{
			inFile>>cost[i];
		}
		inFile.close();
        for(usint i=eta_start;i<eta_end;i++)
        {
            double eta=(eta_0+i*0.02)*numV;
            for(auto alg:algs)
            {
                cout<<"Running alg: "<<alg_arr[alg]<<", at eta = "<<eta_0+i*0.02<<", dataset = "<<dataset[k]<<", # node = "<<numV<<", eta = "<<eta<<endl;
                TAlg Alg(R_graph);
                Alg.set_parameters(model);
                vector<uint32_t> seeds;
                Timer Alg_time("Alg_time");
		        if(alg==0)
		        // BCGC in TKDE
		        seeds = Alg.mine();
                auto run_time=Alg_time.get_total_time();
		        auto memory=getProcMemory();
                ofstream out_seeds("./results/seeds"+alg_arr[alg]+dataset[k]+to_string(eta_0+0.02*i)+model+to_string(Rnd_cost)+".txt", ios::out);
                assert((!out_seeds.fail()));
                for(auto node : seeds)
                {  out_seeds<<node<<endl;  }
                out_seeds.close();
                for(auto node : seeds)	
                {
                    total_cost+=cost[node];
                }
                Graph F_graph = GraphBase::load_graph(dataset_dir, 0);
                Alg.set_F_graph(F_graph);
                pair<double,double> simul_Results=Alg.actual_prob_eval(seeds, simRnd);	
                double actual_prob=simul_Results.first;
                double actual_inf=simul_Results.second;
                string results;
                results="("+dataset[k]+", "+to_string(eta_0+0.02*i)+", "+alg_arr[alg]+", "+to_string(total_cost)+", "+to_string(actual_prob)+", "+to_string(run_time)+", "+to_string(actual_inf)+", "+to_string(memory)+")";
                result_bk<<results<<endl;
        	    cout<<results<<endl;
            }
        }

    }
    
	result_bk.close();
	return 0;
}