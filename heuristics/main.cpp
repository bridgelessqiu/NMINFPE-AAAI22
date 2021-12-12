#include "other_function.h"
#include "heuristic.h"
#include "syds.h"
#include<algorithm>
#include<iostream>
#include<cstdlib>
#include<unordered_map>
#include<string>
#include<vector>
#include <fstream>

using namespace std;

int main(int argc, char *argv[])
{   
    //--------------------------------
    //    Command line arguments     -
    //--------------------------------
    string network_name = argv[1]; // google_plus, power_2...

    int herusitc_num = atoi(argv[2]); // 1: the main herusitic, 2: the second one, 3: random select

    string exp_type = argv[3]; // uniform, random

    string k;

    if(exp_type == "uniform") k = argv[4]; //k-uniform threshold

    //-----------------------------------------------------------------
    //    Construct the name of the network file and threshold file   -
    //-----------------------------------------------------------------
    string network_file = "../networks/real/" + network_name + "/" + network_name + ".edges";

    string thresh_file; // the name of the threshold file

    if(exp_type == "uniform")
    {
        thresh_file = "../networks/real/" + network_name + "/" + network_name + "_" + k + "_uniform_thresh.txt";
    }
    else if(exp_type == "random")
    {
       thresh_file = "../networks/real/" + network_name + "/" + network_name + "_random_thresh.txt";
    }
    else cout<<"Unknown experiment type"<<endl;

    string herusitc_name;

    // -------------------------- 
    //     Call the herustic    -
    // --------------------------
    Syds S;
    S = construct_network(network_file, thresh_file);

    vector<int> A_star;

    if(herusitc_num == 4)
    {
        A_star = GreedyFull(S.neighbors, S.threshold);    
        herusitc_name = "GreedyFull";
    }
    else if(herusitc_num == 2)
    {
        A_star = GreedyNP(S.neighbors, S.threshold);    
        herusitc_name = "GreedyNP";
    }
    else if(herusitc_num == 3)
    {
        A_star = GreedySub(S.neighbors, S.threshold);
        herusitc_name = "GreedySub";
    }

    else if(herusitc_num == 1)
    {
        A_star = GreedyThresh(S.neighbors, S.threshold);
        herusitc_name= "GreedyThresh";
    }

    else
    {
        cout<<"unknown heruistics"<<endl;
        return 0;
    }

    cout<<"The Hamming weight: "<<A_star.size()<<endl;


    // ---------------------------------------
    //      The name of the result file      -
    // ---------------------------------------
    ofstream result;
    string file_name;

    if(exp_type == "uniform")
    {    
        file_name = "results/" + exp_type + "_thresh/" + herusitc_name + "/"  + herusitc_name + "_" + network_name + "_" +  k  + ".txt";
    }
    else if(exp_type == "random")
    {
        file_name = "results/" + exp_type + "_thresh/" + herusitc_name + "/" + herusitc_name + "_" + network_name + ".txt";
    }

    result.open(file_name);
    result<<A_star.size()<<"\n";
    result.close();

    return 0;
}
