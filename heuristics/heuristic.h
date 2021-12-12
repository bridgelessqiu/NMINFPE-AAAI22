#ifndef HEURISTIC_H_06191147
#define HEURISTIC_H_06191147 

#include<string>
#include<unordered_map>
#include<vector>

// GreedyFull
std::vector<int> GreedyFull_sub(const std::vector<std::vector<int>>&, std::vector<int>, int, int);
std::vector<int> GreedyFull(const std::vector<std::vector<int>>& , const std::vector<int>&);

// GreedyNP
std::vector<int> GreedyNP_sub(const std::vector<std::vector<int>>& , std::vector<int> , int , int);
std::vector<int> GreedyNP(const std::vector<std::vector<int>>& , const std::vector<int>& );

// GreedySub
std::vector<int> GreedySub_sub(const std::vector<std::vector<int>>& neighbors, std::vector<int> residual_threshold, int u, int current_opt);
std::vector<int> GreedySub(const std::vector<std::vector<int>>& neighbors, const std::vector<int>& threshold);

// GreedyThresh
std::vector<int> GreedyThresh_sub(const std::vector<std::vector<int>>& neighbors, std::vector<int> residual_threshold, int u, int current_opt);
std::vector<int> GreedyThresh(const std::vector<std::vector<int>>& neighbors, const std::vector<int>& threshold);

#endif
