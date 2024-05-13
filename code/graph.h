#pragma once
#include "FileCtrl.h"
#include "serialize.h"
#include "CommonStruc.h"
#include <fstream>
#include <iostream>
#include <tuple>
using std::tuple;
using namespace std;
class GraphBase
{
public:
	/// Format the input for future computing, which is much faster for loading. Vector serialization is used.
	static void format_graph(const string filename, const bool isReverse)
	{
		size_t numV, numE;
		uint32_t srcId, dstId;
		float weight = 0.0;
		ifstream infile(filename.c_str());
		if (!infile.is_open())
		{
			cout << "The file \"" + filename + "\" can NOT be opened\n";
			return;
		}
		infile >> numV >> numE;
		if(isReverse)
		{
			Graph vecGRev(numV);
			vector<size_t> vecInDeg(numV);
			for (auto i = numE; i--;)
			{
				infile >> srcId >> dstId;
				vecGRev[dstId].push_back(Edge(srcId, weight));
			}
			infile.close();
			for (size_t idx = 0; idx < numV; idx++)
			{
				vecInDeg[idx] = vecGRev[idx].size();
			}
			auto idx = 0;
			for (auto& inNbrs : vecGRev)
			{
			if (inNbrs.empty()) continue; // Skip if there is no in-neighbor.
			weight = (float)1.0 / vecInDeg[idx++];
			for (auto& inNbr : inNbrs)
			{
				get<1>(inNbr) = weight;
			}
			}
			TIO::save_graph_struct(filename, vecGRev, true);
		}//if Reverse
		else  // not reverse
		{
			Graph vecG(numV);
			vector<size_t> vecInDeg(numV);
			for (auto i = numE; i--;)
			{
				infile >> srcId >> dstId;
				pair<uint32_t, float> e=make_pair(dstId, weight);
				vecG[srcId].push_back(e);
				vecInDeg[dstId]+=1;
			}
			infile.close();
			auto idx = 0;
			for (auto& outNbrs : vecG)
			{
			if (outNbrs.empty()) continue; // Skip if there is no in-neighbor.
			weight = (float)1.0 / vecInDeg[idx++];
			for (auto& outNbr : outNbrs)
			{
				get<1>(outNbr) = weight;
			}
			}
			TIO::save_graph_struct(filename, vecG, false);
		}

		cout << "The graph is formatted!" <<  endl;
	}

	/// Load graph via vector deserialization.
	static Graph load_graph(const  string graphName, const bool isReverse)
	{
		Graph graph;
		TIO::load_graph_struct(graphName, graph, isReverse);

    if (isReverse)
    {
      // Reverse graph
      for (auto& nbrs : graph)
      {
        for (auto& nbr : nbrs)
        {
          get<1>(nbr) = float(1.0 / nbrs.size());
        }
      }
    }
    else
    {
      // Forward graph
      vector<uint32_t> vecInDeg(graph.size());
      for (auto& nbrs : graph)
      {
        for (auto& nbr : nbrs)
        {
          vecInDeg[get<0>(nbr)]++;
        }
      }
      for (auto& nbrs : graph)
      {
        for (auto& nbr : nbrs)
        {
          get<1>(nbr) = (float)1.0 / vecInDeg[get<0>(nbr)];
        }
      }
    }
		return graph;
	}

};//cls
