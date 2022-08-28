/*************************************************************************
    > File Name: algo.cc
    > Author: anryyang
    > Mail: anryyang@gmail.com 
    > Created Time: Thu 28 Sep 2017 02:26:32 PM
 ************************************************************************/

#include "algo.h"
#include <math.h>
#include <iostream>
#include <queue>
#include <set>
#include <list>
#include <algorithm>
#include "mtwist.h"
#include <random>
// #include "sparsehash/dense_hash_map"

using namespace std;

double monteCarlo(uint src, uint tgt, uint len_walk, uint64 n_walk, const Graph& graph){
	fastSrand();
	double x=0;

    for(uint len=0; len<len_walk; len++){
	    for(uint64 i=0; i<n_walk; i++){
	    	uint cur = src;
	        for(uint j=0; j<len; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint k = fastRand()%deg;
	        	cur = graph.m_edges[cur][k];
	        }
	        if(cur==tgt){
	        	x += 1.0;
	    	}
	    }
	}
	return x*1.0/n_walk/graph.getDeg(tgt);
}

void runTPC(uint src, uint tgt, uint len_walk, double eps, const Graph& graph){
	fastSrand();
    for(uint len=1; len<len_walk; len++){
    	double beta = pow(0.95,len); // this value is unknown
    	uint64 n_walk = 20000*(pow(len_walk,1.5)*sqrt(beta)/eps + pow(len_walk,3)*pow(beta,1.5)/eps/eps);
	    for(uint64 i=0; i<n_walk; i++){
	    	uint cur = src;
	        for(uint j=0; j<len/2; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint k = fastRand()%deg;
	        	cur = graph.m_edges[cur][k];
	        }
	    }
	    for(uint64 i=0; i<n_walk; i++){
	    	uint cur = tgt;
	        for(uint j=0; j<len/2; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint k = fastRand()%deg;
	        	cur = graph.m_edges[cur][k];
	        }
	    }
	    for(uint64 i=0; i<n_walk; i++){
	    	uint cur = src;
	        for(uint j=0; j<len/2; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint k = fastRand()%deg;
	        	cur = graph.m_edges[cur][k];
	        }
	    }
	    for(uint64 i=0; i<n_walk; i++){
	    	uint cur = tgt;
	        for(uint j=0; j<len/2; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint k = fastRand()%deg;
	        	cur = graph.m_edges[cur][k];
	        }
	    }
	}
}

double runMC2(uint src, uint tgt, uint64 n_walk, const Graph& graph){
	fastSrand();
	double x=0;
	uint source=src;
	uint target=tgt;
	if(graph.getDeg(src)>graph.getDeg(tgt)){
		source = tgt;
		target = src;
	}

    for(uint64 i=0; i<n_walk; i++){
    	uint cur = source;
        while(1){
        	uint deg = graph.m_deg[cur];
        	uint k = fastRand()%deg;
        	uint cur2 = graph.m_edges[cur][k];
        	if(cur2==target){
        		if(cur==source){
        			x++;
        		}
        		break;
        	}
        	cur = cur2;
        }
    }
	return x*1.0/n_walk;
}

double runHAY(uint src, uint tgt, uint64 n_walk, const Graph& graph){
	uint root = 0;

    std::vector<uint> nodes;
    for(uint i=0; i<graph.getN(); i++){
		if(i!=root){
        	nodes.push_back(i);
		}
    }

	int x=0;

	for(uint64 i=0;i<n_walk;i++){
		// cout << i << " walk" << endl;
		std::unordered_map<uint, int> Ti;
		Ti[root]=1;

		int flag=0;

		auto rng = std::default_random_engine {};
		std::shuffle(nodes.begin(), nodes.end(), rng);
		for(const auto& v: nodes){
			// cout << v << " node" << endl;
			std::vector<uint> path;
			uint cur = v;
			path.push_back(cur);
			int j=0;
			while(1){
				if(Ti.find(cur) != Ti.end()){
					if(cur!=v){
						// for(const auto& u: path){
						// 	Ti[u]=1;
						// }
						for(int z=0;z<path.size();z++){
							uint u = path[z];
							Ti[u]=1;
							if(z<path.size()-1){
								if(u==src && path[z+1]==tgt){
									flag=1;
								}
								if(u==tgt && path[z+1]==src){
									flag=1;
								}
								
							}
						}
					}
					break;
				}
				j++;
				uint deg = graph.m_deg[cur];
				uint k = fastRand()%deg;
				uint next = graph.m_edges[cur][k];
				path.push_back(next);
				cur=next;
				if(j>=graph.getN()-1){
					break;
				}
			}
		}

		x+=flag;

		// if(Ti[src]==1 && Ti[tgt]==1){
		// 	x++;
		// }
	}

	return x*1.0/n_walk;
}

double adaptiveMonteCarlo(uint src, uint tgt, std::vector<double>& svec, std::vector<double>& tvec, const Graph& graph, Config& config){
	fastSrand();
	uint len_walk = config.lenwalk;
	uint tau = config.tau;

	uint64 eta_star_s = config.eta_star_s;
	uint64 eta_star_t = config.eta_star_t;
	uint64 eta_s = uint64(max(eta_star_s*1.0/pow(2,tau-1), 1.0));
	uint64 eta_t = uint64(max(eta_star_t*1.0/pow(2,tau-1), 1.0));

	double mu_ss = 0;
	double mu_st = 0;
	double mu_tt = 0;

	for(uint i=0; i<tau; i++){

		mu_ss = 0;
		mu_st = 0;
		double sigma_ss = 0;
		double sigma_st = 0;
	    for(uint64 k=0; k<eta_s; k++){
	    	uint cur = src;
	    	double xss = 0;
	    	double xst = 0;
	        for(uint j=0; j<len_walk; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint l = fastRand()%deg;
	        	cur = graph.m_edges[cur][l];
	        	xss += svec[cur];
	        	xst += tvec[cur];
	        }
	        if(k==0){
	        	sigma_ss=0;
	        	sigma_st=0;
	        }
	        else{
	        	sigma_ss = (k-1)*sigma_ss*1.0/k + pow(xss-mu_ss,2)*1.0/(k+1);
	        	sigma_st = (k-1)*sigma_st*1.0/k + pow(xst-mu_st,2)*1.0/(k+1);
	    	}
	        mu_ss = (k*mu_ss+xss)*1.0/(k+1);
	        mu_st = (k*mu_st+xst)*1.0/(k+1);
	    }

	    sigma_ss = (eta_s-1)*sigma_ss/eta_s;
	    sigma_st = (eta_s-1)*sigma_st/eta_s;
	    eta_s = eta_s*2;

		mu_tt = 0;
		double sigma_tt = 0;
	    for(uint64 k=0; k<eta_t; k++){
	    	uint cur = tgt;
	    	double xtt = 0;
	        for(uint j=0; j<len_walk; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint l = fastRand()%deg;
	        	cur = graph.m_edges[cur][l];
	        	xtt += tvec[cur];
	        }
	        if(k==0){
	        	sigma_tt=0;
	        }
	        else{
	        	sigma_tt = (k-1)*sigma_tt*1.0/k + pow(xtt-mu_tt,2)*1.0/(k+1);
	        }
	        mu_tt = (k*mu_tt+xtt)*1.0/(k+1);
	    } 

	    sigma_tt = (eta_t-1)*sigma_tt/eta_t;
	    eta_t = eta_t*2;

	    double eps_f = config.calErr(eta_s, eta_t, sigma_ss, sigma_st, sigma_tt);
	    cout << "error:" << eps_f << endl;
	    if(eps_f<config.epsilon/2.0){
	    	break;
	    }
	}

	double rf = mu_ss/graph.getDeg(src) + mu_tt/graph.getDeg(tgt) -2*mu_st/graph.getDeg(tgt);
	return rf;
}


double newAdaptiveMonteCarlo(uint src, uint tgt, std::vector<double>& svec, std::vector<double>& tvec, const Graph& graph, Config& config){
	fastSrand();
	uint len_walk = config.lenwalk;
	uint tau = config.tau;

	uint64 eta_star_all = config.eta_star_all;
	uint64 eta_all = uint64(max(eta_star_all*1.0/pow(2,tau-1), 1.0));

	double mu_all = 0;

	for(uint i=0; i<tau; i++){
		mu_all = 0;
		double sigma_all = 0;
	    for(uint64 k=0; k<eta_all; k++){
	    	uint cur = src;
	    	double xss = 0;
	    	double xst = 0;
	        for(uint j=0; j<len_walk; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint l = fastRand()%deg;
	        	cur = graph.m_edges[cur][l];
	        	xss += svec[cur]*1.0;
	        	xst += tvec[cur]*1.0;
	        }

	    	cur = tgt;
	    	double xtt = 0;
			double xts = 0;
	        for(uint j=0; j<len_walk; j++){
	        	uint deg = graph.m_deg[cur];
	        	uint l = fastRand()%deg;
	        	cur = graph.m_edges[cur][l];
	        	xtt += tvec[cur]*1.0;
				xts += svec[cur]*1.0;
	        }

			double mu_update = xss/graph.m_deg[src] + xtt/graph.m_deg[tgt] - xst/graph.m_deg[tgt]-xts/graph.m_deg[src];
			mu_all += mu_update;
			sigma_all=sigma_all+pow(mu_update,2);
	    }

		mu_all = mu_all/eta_all;
		sigma_all = sigma_all/eta_all-pow(mu_all,2);
	    eta_all = eta_all*2;

	    double eps_f = config.calNewErr(eta_all, sigma_all);
	    // cout << "error:" << eps_f << endl;
	    if(eps_f<config.epsilon/2.0){
	    	break;
	    }
	}

	double rf = mu_all;
	return rf;
}


double powerIter(uint src, uint tgt, std::vector<double>& svec, std::vector<double>& tvec, double& smax, double& tmax, const Graph& graph, Config& config, uint& ellb){
    svec[src] = 1.0;
    std::set<uint> S;
    S.insert(src);

	tvec[tgt] = 1.0;
    std::set<uint> T;
    T.insert(tgt);

    ellb=0;
    double pi=0;

    uint ell = config.ell;
    double epsilon = config.epsilon;
    uint ds = graph.getDeg(src);
    uint dt = graph.getDeg(tgt);
    double delta = config.delta;
    uint tau = config.tau;


    for(uint i=0;i<ell;i++){

    		smax = 0;
    		tmax = 0;
    		set<uint> tmpS;
    		set<uint> tmpT;

	    	for(std::set<uint>::iterator it=S.begin(); it!=S.end(); ++it){
		        uint v = *it;

		        double residue = svec[v];
		        if(v==src) pi+=residue/(double)graph.m_deg[v];
				svec[v]=0;
				for(const auto& u: graph.m_edges[v]){
					tmpS.insert(u);
					double update = residue/(double)graph.m_deg[u];
					svec[u] += update;
					if(svec[u]>smax){
						smax = svec[u];
					}
				}
	    	}
	    	S = tmpS;

	    	for(std::set<uint>::iterator it=T.begin(); it!=T.end(); ++it){
		        uint v = *it;

		        double residue = tvec[v];
		        if(v==tgt) pi+=residue/(double)graph.m_deg[v];
		        if(v==src) pi-=2*residue/(double)graph.m_deg[v];
				tvec[v]=0;
				for(const auto& u: graph.m_edges[v]){
					tmpT.insert(u);
					double update = residue/(double)graph.m_deg[u];
					tvec[u] += update;
					if(tvec[u]>tmax){
						tmax = tvec[u];
					}
				}
	    	}
	    	T = tmpT;

	    	ellb+=1;

	    	double picost=0;
	    	for(std::set<uint>::iterator it=S.begin(); it!=S.end(); ++it){
	    		uint v = *it;
	    		picost+=(double)graph.m_deg[v];
	    	}

	    	for(std::set<uint>::iterator it=T.begin(); it!=T.end(); ++it){
	    		uint v = *it;
	    		picost+=(double)graph.m_deg[v];
	    	}

	    	uint lenwalk = ell-ellb;
	    	smax = 2*smax;
	    	tmax = 2*tmax;
	    	int64 eta_star_s = int64(8*log(3*(tau+1)/delta)*lenwalk*lenwalk*pow(max(smax/ds,tmax/dt),2)/epsilon/epsilon);
        	int64 eta_star_t = int64(8*log(3*(tau+1)/delta)*lenwalk*lenwalk*tmax*tmax/epsilon/epsilon/dt/dt);
        	double mccost = (eta_star_s + eta_star_t);

	    	if(picost>mccost){
	    		break;
	    	}
    }

    return pi;
}


double tunePowerIter(uint src, uint tgt, std::vector<double>& svec, std::vector<double>& tvec, double& smax, double& tmax, const Graph& graph, Config& config, uint& ellb){
    svec[src] = 1.0;
    std::set<uint> S;
    S.insert(src);

	tvec[tgt] = 1.0;
    std::set<uint> T;
    T.insert(tgt);

    // ellb=0;
    double pi=0;

    // uint ell = config.ell;
    double epsilon = config.epsilon;
    uint ds = graph.getDeg(src);
    uint dt = graph.getDeg(tgt);
    double delta = config.delta;
    uint tau = config.tau;

	// uint max_ellb = 0; 

    for(uint i=0;i<ellb;i++){

			// if(max_ellb!=0 && ellb>=max_ellb){
			// 	break;
			// }

    		smax = 0;
    		tmax = 0;
    		set<uint> tmpS;
    		set<uint> tmpT;

	    	for(std::set<uint>::iterator it=S.begin(); it!=S.end(); ++it){
		        uint v = *it;

		        double residue = svec[v];
		        if(v==src) pi+=residue/(double)graph.m_deg[v];
				svec[v]=0;
				for(const auto& u: graph.m_edges[v]){
					tmpS.insert(u);
					double update = residue/(double)graph.m_deg[u];
					svec[u] += update;
					if(svec[u]>smax){
						smax = svec[u];
					}
				}
	    	}
	    	S = tmpS;

	    	for(std::set<uint>::iterator it=T.begin(); it!=T.end(); ++it){
		        uint v = *it;

		        double residue = tvec[v];
		        if(v==tgt) pi+=residue/(double)graph.m_deg[v];
		        if(v==src) pi-=2*residue/(double)graph.m_deg[v];
				tvec[v]=0;
				for(const auto& u: graph.m_edges[v]){
					tmpT.insert(u);
					double update = residue/(double)graph.m_deg[u];
					tvec[u] += update;
					if(tvec[u]>tmax){
						tmax = tvec[u];
					}
				}
	    	}
	    	T = tmpT;

	    	// ellb+=1;

	    	// double picost=0;
	    	// for(std::set<uint>::iterator it=S.begin(); it!=S.end(); ++it){
	    	// 	uint v = *it;
	    	// 	picost+=(double)graph.m_deg[v];
	    	// }

	    	// for(std::set<uint>::iterator it=T.begin(); it!=T.end(); ++it){
	    	// 	uint v = *it;
	    	// 	picost+=(double)graph.m_deg[v];
	    	// }

	    	// uint lenwalk = ell-ellb;
	    	// smax = 2*smax;
	    	// tmax = 2*tmax;
	    	// int64 eta_star_s = int64(8*log(3*(tau+1)/delta)*lenwalk*lenwalk*pow(max(smax/ds,tmax/dt),2)/epsilon/epsilon);
        	// int64 eta_star_t = int64(8*log(3*(tau+1)/delta)*lenwalk*lenwalk*tmax*tmax/epsilon/epsilon/dt/dt);
        	// double mccost = (eta_star_s + eta_star_t);

	    	// if(picost>mccost && max_ellb==0){
	    	// 	// break;
			// 	max_ellb = ellb+4;
	    	// }
    }

    return pi;
}


double smm(uint ell, uint src, uint tgt, const Graph& graph){
    // std::unordered_map<int,double> sppr;
    // std::unordered_map<int,double> tppr;
    std::vector<double> svec(graph.getN(),0);
    std::vector<double> tvec(graph.getN(),0);

    svec[src] = 1.0;
    std::set<uint> S;
    S.insert(src);

	tvec[tgt] = 1.0;
    std::set<uint> T;
    T.insert(tgt);

    double pi=0;

    for(uint i=0;i<ell;i++){
    		// cout << i << " " << S.size() << " " << T.size() << endl;

    		set<uint> tmpS;
    		set<uint> tmpT;

	    	for(std::set<uint>::iterator it=S.begin(); it!=S.end(); ++it){
		        uint v = *it;

		        double residue = svec[v];
		        if(v==src) pi+=residue/(double)graph.m_deg[v];
		        if(v==tgt) pi-=2*residue/(double)graph.m_deg[v];
				svec[v]=0;
				for(const auto& u: graph.m_edges[v]){
					tmpS.insert(u);
					double update = residue/(double)graph.m_deg[v];
					svec[u] += update;
				}
	    	}
	    	S = tmpS;

	    	for(std::set<uint>::iterator it=T.begin(); it!=T.end(); ++it){
		        uint v = *it;

		        double residue = tvec[v];
		        if(v==tgt) pi+=residue/(double)graph.m_deg[tgt];
				tvec[v]=0;
				for(const auto& u: graph.m_edges[v]){
					tmpT.insert(u);
					double update = residue/(double)graph.m_deg[v];
					tvec[u] += update;
				}
	    	}
	    	T = tmpT;
    }

    return pi;
}