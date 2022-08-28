/*************************************************************************
    > File Name: algo.h
    > Author: anryyang
    > Mail: anryyang@gmail.com 
    > Created Time: Thu 28 Sep 2017 02:27:17 PM
 ************************************************************************/

#include<vector>
#include <unordered_map>
#include <map>
#include "graph.h"
#include "sampler.h"

double powerIter(uint src, uint tgt, std::vector<double>& svec, std::vector<double>& tvec, double& smax, double& tmax, const Graph& graph, Config& config, uint& ellb);
double tunePowerIter(uint src, uint tgt, std::vector<double>& svec, std::vector<double>& tvec, double& smax, double& tmax, const Graph& graph, Config& config, uint& ellb);
double monteCarlo(uint src, uint tgt, uint len_walk, uint64 n_walk, const Graph& graph);
double smm(uint ell, uint src, uint tgt, const Graph& graph);
double adaptiveMonteCarlo(uint src, uint tgt, std::vector<double>& svec, std::vector<double>& tvec, const Graph& graph, Config& config);
double newAdaptiveMonteCarlo(uint src, uint tgt, std::vector<double>& svec, std::vector<double>& tvec, const Graph& graph, Config& config);
void runTPC(uint src, uint tgt, uint len_walk, double eps, const Graph& graph);
double runMC2(uint src, uint tgt, uint64 n_walk, const Graph& graph);
double runHAY(uint src, uint tgt, uint64 n_walk, const Graph& graph);