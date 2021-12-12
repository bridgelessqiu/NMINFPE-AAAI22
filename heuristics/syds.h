#ifndef SYDS_H_06210354
#define SYDS_H_06210354

#include<vector>
#include<unordered_map>
#include<string>

struct Syds
{
    std::vector<std::vector<int>> neighbors;

    std::vector<int> threshold;
    
    std::unordered_map<std::string, int> name_id_mapping; //[name : id]
};

#endif
