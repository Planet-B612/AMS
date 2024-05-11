#pragma once

#if !defined(DSFMT_MEXP)
#ifdef __GNUC__
#define DSFMT_MEXP 19937
#endif
#endif
/// Node list
typedef std::vector<uint32_t> Nodelist;
/// Edge structure: neighbor id, the edge weight
typedef std::pair<uint32_t, float> Edge;
/// Edgelist structure from one source/target node
typedef std::vector<Edge> Edgelist;
/// Graph structure
typedef std::vector<Edgelist> Graph;
/// One forward reachable set
typedef std::vector<uint32_t> FRset;
/// A set of forward reachable sets
typedef std::vector<FRset> FRsets;
/// One reverse reachable set
typedef std::vector<uint32_t> RRset;
/// A set of reverse reachable sets
typedef std::vector<RRset> RRsets;
typedef std::vector<FRsets> FRcollection;
/// Define the way of storing results
typedef std::tuple<double, double, double, double> Res;
typedef std::vector<Res> vecRes;

using namespace std;
// size_t is "unsigned (long) int", 0 -- 4 294 967 295
// uint32_t is the same as size_t, while size_t is more suitable to sizeof operator
typedef unsigned int uint;
typedef unsigned short int usint;
typedef unsigned char uint8;
typedef char int8;
typedef long long int64;
typedef unsigned long long uint64;
typedef pair<int, int> ipair;
typedef pair<double, double> dpair;

/// Cascade models: IC, LT
enum CascadeModel { IC, LT };
