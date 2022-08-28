/*************************************************************************
    > File Name: bere.cc
    > Author: anryyang
    > Mail: anryyang@gmail.com 
    > Created Time: Thu 28 Sep 2021 02:19:39 PM
 ************************************************************************/

#include<iostream>
#include<boost/program_options.hpp>
#include<string>

#include "graph.h"
#include "algo.h"

using namespace std;
namespace po = boost::program_options;

typedef pair<uint, uint> ipair;
#define MP make_pair

namespace { 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
 
} // namespace

Config parseParams(int argc, char** argv){
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("data-folder,f", po::value<string>()->required(), "graph data folder")
        ("graph-name,g", po::value<string>()->required(), "graph file name")
        ("algo,a", po::value<string>()->required(), "algorithm name")
        ("edge-seed,s", po::value<int>()->default_value(0), "use edge seeds")
        ("epsilon,e", po::value<double>()->default_value(0.5), "epsilon")
        ("delta,d", po::value<double>()->default_value(0.01), "failure probability")
        ("tau,t", po::value<uint>()->default_value(1), "tau")
        ("num,n", po::value<uint>()->default_value(5), "number of queries")
    ;

    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc),  vm); // can throw 
    po::notify(vm);

    Config config;

    if (vm.count("help")){
        cout << desc << '\n';
        exit(0);
    }
    if (vm.count("data-folder")){
        config.strFolder = vm["data-folder"].as<string>();
    }
    if (vm.count("graph-name")){
        config.strGraph = vm["graph-name"].as<string>();
    }
    if (vm.count("algo")){
        config.strAlgo = vm["algo"].as<string>();
    }
    if (vm.count("edge-seed")){
        config.edgeSeed = vm["edge-seed"].as<int>();
    }
    if (vm.count("epsilon")){
        config.epsilon = vm["epsilon"].as<double>();
    }
    if (vm.count("delta")){
        config.delta = vm["delta"].as<double>();
    }
    if (vm.count("tau")){
        config.tau = vm["tau"].as<uint>();
    }
    if (vm.count("num")){
        config.numQuery = vm["num"].as<uint>();
    }

    return config;
}

void loadSeed(string folder, string file_name, int stype, int count, const Graph& graph, vector<ipair>& vecSeeds, vector<double>& vecER){
    string strSeed="";
    if(stype>0){
        strSeed = "/edge_seeds.txt";
    }
    else{
        strSeed = "/rand_seeds.txt";
    }

    FILE *fin = fopen((folder + "/" + file_name + strSeed).c_str(), "r");
    uint s, t;
    double st;
    // vector<ipair> seeds;
    int i=0;
    while (fscanf(fin, "%d %d %lf", &s, &t, &st) != EOF) {
        if(graph.getDeg(s)<=1||graph.getDeg(t)<=1){
            continue;
        }
        ipair p = MP(s,t);
        vecSeeds.push_back(p);
        vecER.push_back(st);
        i++;
        if(i>=count)
            break;
    }
    fclose(fin);
    // return seeds;
}

int main(int argc, char **argv){
    Config config;
    try{
        config = parseParams(argc, argv);
        config.check();
    }
    catch (const exception &ex){
        cerr << ex.what() << '\n';
        return ERROR_IN_COMMAND_LINE;
    }

    config.display();

    int query_count = config.numQuery;
    Graph graph(config.strFolder, config.strGraph);
    vector<ipair> seeds;
    vector<double> exact_ers;
    loadSeed(config.strFolder, config.strGraph, config.edgeSeed, query_count, graph, seeds, exact_ers);
    cout << "read " << seeds.size() << " seeds Done!" << endl;

    config.lambda = graph.getLambda();

    vector<double> errors(exact_ers.size(),0);

    if(config.strAlgo==GEER){
        cout << "start bppr with power-method!" << endl;
        Timer tm(1, "BERE");
        std::vector<double> svec(graph.getN(),0);
        std::vector<double> tvec(graph.getN(),0);
        int q=0;
        for(const auto& p: seeds){
            uint s = p.first;
            uint t = p.second;
            cout << "-------------processing " << s << "," << t << endl;
            cout << graph.getDeg(s) << " " << graph.getDeg(t) << endl;

            svec.clear();
            tvec.clear();

            config.setSMM(graph.getDeg(s), graph.getDeg(t));
            uint ellb=0;
            double smax=0;
            double tmax=0;
            double rb = powerIter(s, t, svec, tvec, smax, tmax, graph, config, ellb);


            double rf = 0;
            if(config.ell > ellb){
                uint ellf = config.ell - ellb;
                cout << "ell_f:" << ellf << endl;
                config.setGEER(smax, tmax, ellf, graph.getDeg(s), graph.getDeg(t));
                // rf = adaptiveMonteCarlo(s, t, svec, tvec, graph, config);
                rf = newAdaptiveMonteCarlo(s, t, svec, tvec, graph, config);
            }

            double er = rf+rb;
            errors[q] = abs(exact_ers[q]-er);
            cout << "ER:" << er << " EXACT-ER: " << exact_ers[q] << " Err: " << errors[q] << endl;
            q++;
        }
    }
    else if(config.strAlgo==TP){
        Timer tm(1, "TP");
        config.setTP();
        cout << "lenwalk: " << config.lenwalk << endl;
        cout << "numwalks: " << config.numwalks << endl;
        srand(time(NULL));
        int q=0;
        for(const auto& p: seeds){
            uint s = p.first;
            uint t = p.second;
            cout << "-------------processing " << s << "," << t << endl;
            double ss = monteCarlo(s, s, config.lenwalk, config.numwalks, graph);
            double st = monteCarlo(s, t, config.lenwalk, config.numwalks, graph);
            double tt = monteCarlo(t, t, config.lenwalk, config.numwalks, graph);
            double er = ss+tt-2*st;
            // cout << "ER:" << er << endl;
            errors[q] = abs(exact_ers[q]-er);
            cout << "ER:" << er << " EXACT-ER: " << exact_ers[q] << " Err: " << errors[q] << endl;
            q++;
        }
    }
    else if(config.strAlgo==TPC){
        Timer tm(1, "TPC");
        config.setTPC();
        cout << "lenwalk: " << config.lenwalk << endl;
        srand(time(NULL));
        for(const auto& p: seeds){
            uint s = p.first;
            uint t = p.second;
            cout << "-------------processing " << s << "," << t << endl;
            runTPC(s, t, config.lenwalk, config.epsilon, graph);
        }
    }
    else if(config.strAlgo==MC2){
        Timer tm(1, "MC2");
        double gamma=100;
        for(int i=0; i<exact_ers.size(); i++){
            if(exact_ers[i]<gamma){
                gamma = exact_ers[i]/2;
            }
        }
        // gamma =max(gamma-config.epsilon,1.0/graph.getM());
        // double gamma = 1.0/graph.getM();
        config.setMC2(gamma);
        cout << "numwalks: " << config.numwalks << endl;
        srand(time(NULL));
        int q=0;
        for(const auto& p: seeds){
            uint s = p.first;
            uint t = p.second;
            cout << "-------------processing " << s << "," << t << endl;
            double er = runMC2(s, t, config.numwalks, graph);
            // cout << "ER:" << er << endl;
            errors[q] = abs(exact_ers[q]-er);
            cout << "ER:" << er << " EXACT-ER: " << exact_ers[q] << " Err: " << errors[q] << endl;
            q++;
        }
    }
    else if(config.strAlgo==AMC){
        Timer tm(1, "AMC");
        std::vector<double> svec(graph.getN(),0);
        std::vector<double> tvec(graph.getN(),0);
        srand(time(NULL));
        int q=0;
        for(const auto& p: seeds){
            uint s = p.first;
            uint t = p.second;
            cout << "-------------processing " << s << "," << t << endl;
            cout << graph.getDeg(s) << " " << graph.getDeg(t) << endl;
            config.setAMC(graph.getDeg(s), graph.getDeg(t));
            svec[s]=1;
            tvec[t]=1;
            double er = svec[s]/graph.getDeg(s) + tvec[t]/graph.getDeg(t)-2*svec[t]/graph.getDeg(t);
            // er += adaptiveMonteCarlo(s, t, svec, tvec, graph, config);
            er += newAdaptiveMonteCarlo(s, t, svec, tvec, graph, config);
            svec[s]=0;
            tvec[t]=0;
            // cout << "ER:" << er << endl;
            errors[q] = abs(exact_ers[q]-er);
            cout << "ER:" << er << " EXACT-ER: " << exact_ers[q] << " Err: " << errors[q] << endl;
            q++;
        }
    }
    else if(config.strAlgo==SMM){
        Timer tm(1, "SMM");
        srand(time(NULL));
        int q=0;
        for(const auto& p: seeds){
            uint s = p.first;
            uint t = p.second;
            cout << "-------------processing " << s << "," << t << endl;
            cout << graph.getDeg(s) << " " << graph.getDeg(t) << endl;
            config.setSMM(graph.getDeg(s), graph.getDeg(t));
            cout << "ell: " << config.ell << endl;
            double er = smm(config.ell, s, t, graph);
            // cout << "ER:" << er << endl;
            errors[q] = abs(exact_ers[q]-er);
            cout << "ER:" << er << " EXACT-ER: " << exact_ers[q] << " Err: " << errors[q] << endl;
            q++;
        }
    }
    else if(config.strAlgo==HAY){
        Timer tm(1, "HAY");
        config.setHAY(graph.getM());
        cout << "#walks: " << config.numwalks << endl;
        srand(time(NULL));
        int q=0;
        for(const auto& p: seeds){
            uint s = p.first;
            uint t = p.second;
            cout << "-------------processing " << s << "," << t << endl;
            double er = runHAY(s, t, config.numwalks, graph);
            errors[q] = abs(exact_ers[q]-er);
            cout << "ER:" << er << " EXACT-ER: " << exact_ers[q] << " Err: " << errors[q] << endl;
            q++;
        }
    }

    double sum_err=0;
    double min_err=100;
    double max_err=0;
    for(int q=0; q<errors.size(); q++){
        if(errors[q]>max_err){
            max_err=errors[q];
        }
        if(errors[q]<min_err){
            min_err=errors[q];
        }
        sum_err+=errors[q];
    }

    cout << "min-error: " << min_err << endl;
    cout << "max-error: " << max_err << endl;
    cout << "avg-error: " << sum_err/errors.size() << endl;
    
    cout << Timer::used(1)*1000/query_count << " milli-seconds per query" << endl;
    Timer::show();

    return SUCCESS;
}
