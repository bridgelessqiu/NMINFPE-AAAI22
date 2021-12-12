#include "heuristic.h"
#include<list>
#include<queue>
#include<unordered_map>
#include<cstdlib>
#include<random>
#include<numeric>
#include<algorithm>
#include<vector>
#include<iostream>

using namespace std;

// --------------- MAIN HEURISITC ------------------

/* ------------------------------ *
 * The subrountine for GreedyFull *
 * ------------------------------ */
vector<int> GreedyFull_sub(const vector<vector<int>>& neighbors, vector<int> residual_threshold, int u, int current_opt)
{
    int n = neighbors.size(); // the number of nodes
    vector<bool> state_one(n, 0); // state_one[u] = 1 iff u is in state 1
    vector<int> thresh = residual_threshold; // This records the original thresholds
    vector<int> num_unsati_nei(n, 0); // num_unsati_nei[u] : number of neighbors of u who are selected, but unsatisfied

    // u must be state 1
    state_one[u] = 1;
    int hamming_weight = 0; // will update later

    // Determine the set of nodes that are passively set to 1
    queue<int> que; // contains vertices that are NEWLY set to state 1
    que.push(u);

    while(!que.empty())
    {
        hamming_weight++;
        int v = que.front();
        que.pop();

        for(auto w : neighbors[v]) // closed neighborhood!
        {
            residual_threshold[w]--;
            if(residual_threshold[w] <= 0 && !state_one[w])
            {
                state_one[w] = 1;
                que.push(w);
            }
        }
    }

    int delta = max(0, residual_threshold[u]); // the additional number of vertices to satisfy u
    
    if((delta + hamming_weight) >= current_opt)
    {
        return vector<int> (n + 1, 0);
    }
    
    if(residual_threshold[u] > 0)
    {
        for(auto v : neighbors[u]) // only neigbors of u is considered, since all other vertices in A are satisfied
        {
            num_unsati_nei[v]++; // Kinda sloppy here
        }
    }

    // ------------------------------------------------------------------------
    
    while(delta > 0) // each iteration, we select one vertex to set to state 1
    {
        int v_star = -1; // the next vertex to set to state 1
        int lowest_score = 2 * n;

        for(int v = 0; v != n; ++v)
        {
            if(state_one[v] == 1 || num_unsati_nei[v] <= 0) continue;
            
            if(residual_threshold[v] + hamming_weight >= current_opt) continue;

            if(thresh[v] > neighbors[v].size()) continue; // this should never be called

            int rho_v = 0; // the number of vertices being passively set to 1
            int eps_v = 0; // the decrease of delta

            vector<bool> temp_state_one = state_one;
            vector<int> temp_residual_threshold = residual_threshold;

            que.push(v); //we can reuse que becasue it is always empty after each loop

            temp_state_one[v] = 1;

            while(!que.empty())
            {
                int x = que.front();
                que.pop();

                for(auto y : neighbors[x]) // closed neighborhood!
                {
                    if(y != v && state_one[y] == 1 && temp_residual_threshold[y] > 0) eps_v++;
                    
                    temp_residual_threshold[y]--;

                    if(temp_state_one[y] == 0 && temp_residual_threshold[y] <= 0)
                    {
                        que.push(y);
                        temp_state_one[y] = 1;
                        rho_v++;
                    }
                }
            }

            int sigma_v = max(0, temp_residual_threshold[v]); // the number of additional vertices that need to be selected if we selected v
            if(lowest_score > (rho_v + sigma_v - eps_v))
            {
                v_star = v;
                lowest_score = rho_v + sigma_v - eps_v;
            }
               
        } //end for(auto v : B)

        if(v_star == -1)
        {
            return vector<int> (n + 1, 0); // no candiate vertices are selected
        }
        
        // Set v_star to 1
        vector<bool> temp_state_one = state_one;
        temp_state_one[v_star] = 1;

        que.push(v_star); //we can reuse que becasue it is always empty after each loop

        int eps_star = 0;
        
        while(!que.empty())
        {
            hamming_weight++; // each element in que is a new state-1 vertex
            int x = que.front();
            que.pop();

            for(auto y : neighbors[x]) // closed neighborhood!
            {
                if(y != v_star && state_one[y] == 1 && residual_threshold[y] > 0)
                {
                    eps_star++; 
                }

                residual_threshold[y]--;

                if(y != v_star && state_one[y] == 1 && residual_threshold[y] + 1 > 0 && residual_threshold[y]<=0) // y just got satisfied
                {
                    for(auto z : neighbors[y]) num_unsati_nei[z]--;  
                }

                if(!temp_state_one[y] && residual_threshold[y] <= 0)
                {
                    que.push(y);
                    temp_state_one[y] = 1;
                }
            }
        }

        state_one = temp_state_one;

        int sigma_star = max(0, residual_threshold[v_star]); // the number of additional vertices that need to be selected if we selected v

        delta = delta + sigma_star - eps_star; 

        //if(hamming_weight + delta >= current_opt) // This part might need to be changed
        if(hamming_weight >= current_opt) // This part might need to be changed
        {
            return vector<int> (n + 1, 0);
        }
       
   // -------------------------------------------------------

        if(residual_threshold[v_star] > 0)
        {
            for(auto x : neighbors[v_star])
            {
                num_unsati_nei[x]++;
            }
        }
    } //end while(delta != 0)

    // Construct A
    vector<int>A(hamming_weight);

    int ind = 0;
    for(int w = 0; w < n; ++w)
    {
        if(state_one[w] == 1)
        {
            A[ind] = w;
            ind++;
        }
    }

    return A;
} // end function


/* -------------- *
 *   GreedyFull   *
 * ---------------*/ 
vector<int> GreedyFull(const vector<vector<int>>& neighbors, const vector<int>& threshold)
{
    cout<<"---------------------"<<endl;
    cout<<"-    GreedyFull     -"<<endl;
    cout<<"---------------------"<<endl;

    cout<<"Running ..."<<endl;

    int n = neighbors.size(); // the numbe of nodes
    int obj = n + 1; // the optimal hamming weight
    vector<int> A_star;
    vector<int> C;

    // order vertices by thresholds
    vector<int> sorted_vectices(n);
    iota(sorted_vectices.begin(), sorted_vectices.end(), 0); //[0, 1, ..., n]
    sort(sorted_vectices.begin(), sorted_vectices.end(), [&](int i, int j){return threshold[i] < threshold[j];} );
    
    // Enumerate over the sorted order
    for(auto u : sorted_vectices)
    //for(int u = 0; u < n; ++u)
    {   
        if(threshold[u] > neighbors[u].size()) continue; // u will always be in state 0, this should never happen by default
        if(threshold[u] >= obj) continue;
            
        C = GreedyFull_sub(neighbors, threshold, u, obj); // the optimal fixed point that includes u

        if(obj > C.size())
        {
            obj = C.size();
            A_star = C;
        }
    }
    
    cout<<"Vertices to be set to state one are:"<<endl;
    cout<<"[";
    for(auto v : A_star) cout<<v<<", ";
    cout<<"]"<<endl;
    cout<<"----------\n";

    // ----------------------------------------------
    // -    Check if it is indeed a fixed point     -
    // ----------------------------------------------
    vector<int> state_one(n, 0);
    for(auto u : A_star) state_one[u] = 1;

    for(int u = 0; u < n; ++u)
    {
        int inp = 0;
        for(auto v : neighbors[u])
        {
            inp += state_one[v];
        }
        if(state_one[u] == 1 && inp < threshold[u])
        {
            cout<<"not a fixed point"<<endl;       
        }
        if(state_one[u] == 0 && inp >= threshold[u])
        {
            cout<<"not a fixed point"<<endl;            
        }
    }

    return A_star; // or we can return A_star
}
      
/* -----------------------------*
 * The subrountine for GreedyNP *
 * -----------------------------*/
vector<int> GreedyNP_sub(const vector<vector<int>>& neighbors, vector<int> residual_threshold, int u, int current_opt)
{
    int n = neighbors.size(); // the number of nodes
    vector<bool> state_one(n, 0); // state_one[u] = 1 iff u is in state 1
    vector<int> num_unsati_nei(n, 0); // num_unsati_nei[u] : number of neighbors of u who are selected, but unsatisfied
    vector<int> thresh = residual_threshold;

    // u must be state 1
    state_one[u] = 1;
    int hamming_weight = 0; // will update later

    // Determine the set of nodes that are passively set to 1
    queue<int> que; // contains vertices that are NEWLY set to state 1
    que.push(u);

    while(!que.empty())
    {
        hamming_weight++;
        int v = que.front();
        que.pop();

        for(auto w : neighbors[v]) // closed neighborhood!
        {
            residual_threshold[w]--;
            if(residual_threshold[w] <= 0 && !state_one[w])
            {
                state_one[w] = 1;
                que.push(w);
            }
        }
    }

    int delta = max(0, residual_threshold[u]); // the additional number of vertices to satisfy u
    
    if((delta + hamming_weight) >= current_opt)
    {
        return vector<int> (n + 1, 0);
    }
    
    vector<int> B(n, 0); // candidate vertices. IMPORTNAT: only state-0 vertices with unsatisfied state-1 neighbors can be cadidates.

    if(residual_threshold[u] > 0)
    {
        for(auto v : neighbors[u]) // only neigbors of u is considered, since all other vertices in A are satisfied
        {
            num_unsati_nei[v]++; // Kinda sloppy here
            if(state_one[v] == 0)
            {
                B[v] = 1;
            }
        }
    }

    // ------------------------------------------------------------------------

    while(delta > 0) // each iteration, we select one vertex to set to state 1
    {
        int v_star = -1; // the next vertex to set to state 1
        int lowest_score = 2 * n;

        for(int v = 0; v != B.size(); ++v)
        {
            if(B[v] == 0) continue; // only consider candidate vertices

            if(thresh[v] > neighbors[v].size()) continue; // this should never be called
            
            if(residual_threshold[v] + hamming_weight >= current_opt) continue;

            int sigma_v = residual_threshold[v] - 1; // the number of additional vertices that need to be selected if we selected v 
            int eps_v = num_unsati_nei[v]; // the number of vertices that can satisfy

            if(lowest_score > (sigma_v - eps_v))
            {
                v_star = v;
                lowest_score = sigma_v - eps_v;
            }
               
        } //end for(auto v : B)

        if(v_star == -1)
        {
            return vector<int> (n + 1, 0); // no candiate vertices are selected 
        }
        
        // Set v_star to 1
        vector<bool> temp_state_one = state_one;
        temp_state_one[v_star] = 1;
        B[v_star] = 0; // Important
        que.push(v_star); //we can reuse que becasue it is always empty after each loop

        int eps_star = 0;
        
        while(!que.empty())
        {
            hamming_weight++; // each element in que is a new state-1 vertex
            int x = que.front();
            que.pop();

            for(auto y : neighbors[x]) // closed neighborhood!
            {
                if(y != v_star && state_one[y] == 1 && residual_threshold[y] > 0)
                {
                    eps_star++; 
                }

                residual_threshold[y]--;

                if(y != v_star && state_one[y] == 1 && residual_threshold[y] + 1 > 0 && residual_threshold[y] <= 0) // y just got satisfied
                {
                    for(auto z : neighbors[y]) num_unsati_nei[z]--;  
                }

                if(!temp_state_one[y] && residual_threshold[y] <= 0)
                {
                    que.push(y);
                    temp_state_one[y] = 1;
                    B[y] = 0;
                }
            }
        }

        state_one = temp_state_one;

        int sigma_star = max(0, residual_threshold[v_star]); // the number of additional vertices that need to be selected if we selected v
        
        delta = delta + sigma_star - eps_star; 
        
        if(hamming_weight + delta >= current_opt)
        {
            return vector<int> (n + 1, 0);
        }
       
        // -------------------------------------------------------------------
        // - Remove candiate vertices in B whose neighbors are all satisfied -
        // -------------------------------------------------------------------
        
        for(int w = 0; w < B.size(); ++w)
        {
            if(B[w] == 0) continue;

            if(num_unsati_nei[w] <= 0)
            {
                B[w] = 0;
            }
        }

   // -------------------------------------------------------

        if(residual_threshold[v_star] > 0)
        {
            for(auto x : neighbors[v_star])
            {
                if(!state_one[x]) B[x] = 1;
                num_unsati_nei[x]++;
            }
        }

    } //end while(delta != 0)
    
    // Construct A
    vector<int>A(hamming_weight);

    int ind = 0;
    for(int w = 0; w < n; ++w)
    {
        if(state_one[w] == 1)
        {
            A[ind] = w;
            ind++;
        }
    }
    return A;
} // end function


/* ----------- *
 *  GreedyNP   * 
 * ----------- */
vector<int> GreedyNP(const vector<vector<int>>& neighbors, const vector<int>& threshold)
{

    cout<<"-------------------"<<endl;
    cout<<"-    GreedyNP     -"<<endl;
    cout<<"-------------------"<<endl;

    cout<<"Running ..."<<endl;

    int n = neighbors.size(); // the numbe of nodes
    int obj = n + 1; // the optimal hamming weight
    vector<int> A_star;
    vector<int> C;

    // order vertices by thresholds
    vector<int> sorted_vectices(n);
    iota(sorted_vectices.begin(), sorted_vectices.end(), 0); //[0, 1, ..., n]
    sort(sorted_vectices.begin(), sorted_vectices.end(), [&](int i, int j){return threshold[i] < threshold[j];} );
    
    // Enumerate over the sorted order
    for(auto u : sorted_vectices)
    {
        if(threshold[u] > neighbors[u].size()) continue; // u will always be in state 0, this should never happen by default
        if(threshold[u] >= obj) continue;
            
        C = GreedyNP_sub(neighbors, threshold, u, obj); // the optimal fixed point that includes u

        if(obj >= C.size())
        {
            obj = C.size();
            A_star = C;
        }

    }

    cout<<"Vertices to be set to state one are:"<<endl;
    cout<<"[";
    for(auto v : A_star) cout<<v<<", ";
    cout<<"]"<<endl;
    cout<<"----------\n";

    // ----------------------------------------------
    // -    Check if it is indeed a fixed point     -
    // ----------------------------------------------

    vector<int> state_one(n, 0);
    for(auto u : A_star) state_one[u] = 1;

    for(int u = 0; u < n; ++u)
    {
        int inp = 0;
        for(auto v : neighbors[u])
        {
            inp += state_one[v];
        }

        if(state_one[u] == 1 && inp < threshold[u])
        {
            cout<<"not a fixed point"<<endl;       
        }

        if(state_one[u] == 0 && inp >= threshold[u])
        {
            cout<<"not a fixed point"<<endl;            
        }
    }

    return A_star; // or we can return A_star
}

     

/* ------------------------------------------ *
 * The subrountine for GreedyThresh *
 * ------------------------------------------ */
vector<int> GreedyThresh_sub(const vector<vector<int>>& neighbors, vector<int> residual_threshold, int u, int current_opt)
{
    int n = neighbors.size(); // the number of nodes
    vector<bool> state_one(n, 0); // state_one[u] = 1 iff u is in state 1
    vector<int> num_unsati_nei(n, 0); // num_unsati_nei[u] : number of neighbors of u who are selected, but unsatisfied
    vector<int> thresh = residual_threshold;

    vector<int> sorted_vectices(n);

    iota(sorted_vectices.begin(), sorted_vectices.end(), 0); //[0, 1, ..., n]

    sort(sorted_vectices.begin(), sorted_vectices.end(), [&](int i, int j){return residual_threshold[i] < residual_threshold[j];} );

    // u must be state 1
    state_one[u] = 1;
    int hamming_weight = 0; // will update later

    // Determine the set of nodes that are passively set to 1
    queue<int> que; // contains vertices that are NEWLY set to state 1
    que.push(u);

    while(!que.empty())
    {
        hamming_weight++;
        int v = que.front();
        que.pop();

        for(auto w : neighbors[v]) // closed neighborhood!
        {
            residual_threshold[w]--;
            if(residual_threshold[w] <= 0 && !state_one[w])
            {
                state_one[w] = 1;
                que.push(w);
            }
        }
    }

    int delta = max(0, residual_threshold[u]); // the additional number of vertices to satisfy u
    
    if((delta + hamming_weight) >= current_opt)
    {
        return vector<int> (n + 1, 0);
    }
    
    vector<int> B(n, 0); // candidate vertices. IMPORTNAT: only state-0 vertices with unsatisfied state-1 neighbors can be cadidates.

    if(residual_threshold[u] > 0)
    {
        for(auto v : neighbors[u]) // only neigbors of u is considered, since all other vertices in A are satisfied
        {
            num_unsati_nei[v]++; // Kinda sloppy here
            if(state_one[v] == 0)
            {
                B[v] = 1;
            }
        }
    }

    // ------------------------------------------------------------------------
    int total = 0;

    while(delta > 0) // each iteration, we select one vertex to set to state 1
    {
        int v_star = -1;

        vector<int> candidates;
        for(int i = 0; i != B.size(); ++i) 
        {
            if(B[i] == 1 && thresh[i] <= neighbors[i].size()) candidates.push_back(i);
        }

        if(candidates.size() == 0) 
        {
            return vector<int> (n+1, 0);
        }
       
        // Select the one with the smallest degree
        sort(candidates.begin(), candidates.end(), [&](int i, int j){return residual_threshold[i] < residual_threshold[j];} );

        v_star = candidates[0];

        if(v_star == -1)
        {
            return vector<int> (n + 1, 0); // no candiate vertices are selected 
        }
        
        // Set v_star to 1
        vector<bool> temp_state_one = state_one;
        temp_state_one[v_star] = 1;
        B[v_star] = 0; // Important
        que.push(v_star); //we can reuse que becasue it is always empty after each loop

        int eps_star = 0;
        
        while(!que.empty())
        {
            hamming_weight++; // each element in que is a new state-1 vertex
            int x = que.front();
            que.pop();

            for(auto y : neighbors[x]) // closed neighborhood!
            {
                if(y != v_star && state_one[y] == 1 && residual_threshold[y] > 0)
                {
                    eps_star++; 
                }

                residual_threshold[y]--;

                if(y != v_star && state_one[y] == 1 && residual_threshold[y] + 1 > 0 && residual_threshold[y] <= 0) // y just got satisfied
                {
                    for(auto z : neighbors[y]) num_unsati_nei[z]--;  
                }

                if(!temp_state_one[y] && residual_threshold[y] <= 0)
                {
                    que.push(y);
                    temp_state_one[y] = 1;
                    B[y] = 0;
                }
            }
        }

        state_one = temp_state_one;

        int sigma_star = max(0, residual_threshold[v_star]); // the number of additional vertices that need to be selected if we selected v
        
        delta = delta + sigma_star - eps_star; 
        
        if(hamming_weight + delta >= current_opt)
        {
            return vector<int> (n + 1, 0);
        }
       
        // -------------------------------------------------------------------
        // - Remove candiate vertices in B whose neighbors are all satisfied -
        // -------------------------------------------------------------------
        
        for(int w = 0; w < B.size(); ++w)
        {
            if(B[w] == 0) continue;

            if(num_unsati_nei[w] <= 0)
            {
                B[w] = 0;
            }
        }

   // -------------------------------------------------------

        if(residual_threshold[v_star] > 0)
        {
            for(auto x : neighbors[v_star])
            {
                if(!state_one[x]) B[x] = 1;
                num_unsati_nei[x]++;
            }
        }

    } //end while(delta != 0)
    
    // Construct A
    vector<int>A(hamming_weight);

    int ind = 0;
    for(int w = 0; w < n; ++w)
    {
        if(state_one[w] == 1)
        {
            A[ind] = w;
            ind++;
        }
    }
    return A;
} // end function


/* --------------- *
 *   GreedyThresh  * 
 * --------------- */
vector<int> GreedyThresh(const vector<vector<int>>& neighbors, const vector<int>& threshold)
{

    cout<<"-----------------------"<<endl;
    cout<<"-    GreedyThresh     -"<<endl;
    cout<<"-----------------------"<<endl;

    cout<<"Running ..."<<endl;

    int n = neighbors.size(); // the numbe of nodes
    int obj = n + 1; // the optimal hamming weight
    vector<int> A_star;
    vector<int> C;
   
    vector<int> sorted_vectices(n);

    iota(sorted_vectices.begin(), sorted_vectices.end(), 0); //[0, 1, ..., n]

    sort(sorted_vectices.begin(), sorted_vectices.end(), [&](int i, int j){return threshold[i] < threshold[j];} );
    
    for(auto u : sorted_vectices)
    {
        if(threshold[u] > neighbors[u].size()) continue; // u will always be in state 0, this should never happen by default

        if(threshold[u] >= obj) continue;
         
        C = GreedyThresh_sub(neighbors, threshold, u, obj); // the optimal fixed point that includes u

        if(obj > C.size())
        {
            obj = C.size();
            A_star = C;
        }
    }

    cout<<"Vertices to be set to state one are:"<<endl;
    cout<<"[";
    for(auto v : A_star) cout<<v<<", ";
    cout<<"]"<<endl;
    cout<<"----------\n";

    // ----------------------------------------------
    // -    Check if it is indeed a fixed point     -
    // ----------------------------------------------
    vector<int> state_one(n, 0);
    for(auto u : A_star) state_one[u] = 1;

    for(int u = 0; u < n; ++u)
    {
        int inp = 0;
        for(auto v : neighbors[u])
        {
            inp += state_one[v];
        }
        if(state_one[u] == 1 && inp < threshold[u])
        {
            cout<<"not a fixed point"<<endl;       
        }
        if(state_one[u] == 0 && inp >= threshold[u])
        {
            cout<<"not a fixed point"<<endl;            
        }
    }

    return A_star; 
}

 
/* ------------------------------ *
 * The subrountine for GreedySub  *
 * ------------------------------ */
vector<int> GreedySub_sub(const vector<vector<int>>& neighbors, vector<int> residual_threshold, int u, int current_opt)
{
    int n = neighbors.size(); // the number of nodes
    vector<bool> state_one(n, 0); // state_one[u] = 1 iff u is in state 1
    vector<int> thresh = residual_threshold;
    vector<int> num_unsati_nei(n, 0); // num_unsati_nei[u] : number of neighbors of u who are selected, but unsatisfied

    // u must be state 1
    state_one[u] = 1;
    int hamming_weight = 0; // will update later

    // Determine the set of nodes that are passively set to 1
    queue<int> que; // contains vertices that are NEWLY set to state 1
    que.push(u);

    while(!que.empty())
    {
        hamming_weight++;
        int v = que.front();
        que.pop();

        for(auto w : neighbors[v]) // closed neighborhood!
        {
            residual_threshold[w]--;
            if(residual_threshold[w] <= 0 && !state_one[w])
            {
                state_one[w] = 1;
                que.push(w);
            }
        }
    }

    int delta = max(0, residual_threshold[u]); // the additional number of vertices to satisfy u
    
    if((delta + hamming_weight) >= current_opt)
    {
        return vector<int> (n + 1, 0);
    }
    
    if(residual_threshold[u] > 0)
    {
        for(auto v : neighbors[u]) // only neigbors of u is considered, since all other vertices in A are satisfied
        {
            num_unsati_nei[v]++; // Kinda sloppy here
        }
    }

    // ------------------------------------------------------------------------
    
    while(delta > 0) // each iteration, we select one vertex to set to state 1
    {
        int v_star = -1; // the next vertex to set to state 1
        int lowest_score = 2 * n;

        for(int v = 0; v < n; ++v)
        {
            if(state_one[v] == 1) continue;

            if(num_unsati_nei[v] <= 0) continue;
            
            if(residual_threshold[v] + hamming_weight >= current_opt) continue;

            if(thresh[v] > neighbors[v].size()) continue;

            int rho_v = 0;
            int eps_v = 0;

            vector<bool> temp_state_one = state_one;
            vector<int> temp_residual_threshold = residual_threshold;

            que.push(v); //we can reuse que becasue it is always empty after each loop

            temp_state_one[v] = 1;

            while(!que.empty())
            {
                int x = que.front();
                que.pop();

                for(auto y : neighbors[x]) // closed neighborhood!
                {
                    if(y != v && state_one[y] == 1 && temp_residual_threshold[y] > 0) eps_v++;
                    
                    temp_residual_threshold[y]--;

                    if(temp_state_one[y] == 0 && temp_residual_threshold[y] <= 0)
                    {
                        que.push(y);
                        temp_state_one[y] = 1;
                        rho_v++;
                    }
                }
            }

            int sigma_v = max(0, temp_residual_threshold[v]); // the number of additional vertices that need to be selected if we selected v 

            if(lowest_score > (rho_v + sigma_v - eps_v))
            {
                v_star = v;
                lowest_score = rho_v + sigma_v - eps_v;
            }
               
        } //end for(auto v : B)

        if(v_star == -1) return vector<int> (n + 1, 0); // no candiate vertices are selected 
        
        // Set v_star to 1
        vector<bool> temp_state_one = state_one;
        temp_state_one[v_star] = 1;
        que.push(v_star); //we can reuse que becasue it is always empty after each loop

        int eps_star = 0;
        
        while(!que.empty())
        {
            hamming_weight++; // each element in que is a new state-1 vertex
            int x = que.front();
            que.pop();

            for(auto y : neighbors[x]) // closed neighborhood!
            {
                if(y != v_star && state_one[y] == 1 && residual_threshold[y] > 0)
                {
                    eps_star++; 
                }

                residual_threshold[y]--;

                if(y != v_star && state_one[y] == 1 && residual_threshold[y] + 1 > 0 && residual_threshold[y]<=0) // y just got satisfied
                {
                    for(auto z : neighbors[y]) num_unsati_nei[z]--;  
                }

                if(!temp_state_one[y] && residual_threshold[y] <= 0)
                {
                    que.push(y);
                    temp_state_one[y] = 1;
                }
            }
        }

        state_one = temp_state_one;

        int sigma_star = max(0, residual_threshold[v_star]); // the number of additional vertices that need to be selected if we selected v

        delta = delta + sigma_star - eps_star; 

        if(hamming_weight + delta >= current_opt)
        {
            return vector<int> (n + 1, 0);
        }
       
   // -------------------------------------------------------

        if(residual_threshold[v_star] > 0)
        {
            for(auto x : neighbors[v_star])
            {
                num_unsati_nei[x]++;
            }
        }
    } //end while(delta != 0)

    // Construct A
    vector<int>A(hamming_weight);

    int ind = 0;
    for(int w = 0; w < n; ++w)
    {
        if(state_one[w] == 1)
        {
            A[ind] = w;
            ind++;
        }
    }

    return A;
} // end function

/* ------------- *
 *   GreedySub   * 
 * ------------- */
vector<int> GreedySub(const vector<vector<int>>& neighbors, const vector<int>& threshold)
{

    cout<<"--------------------"<<endl;
    cout<<"-    GreedySub     -"<<endl;
    cout<<"--------------------"<<endl;

    cout<<"Running ..."<<endl;

    int n = neighbors.size(); // the numbe of nodes
    int obj = n + 1; // the optimal hamming weight
    vector<int> A_star;
    vector<int> C;
    vector<int> selected(n, 0);
    
    for(int u = 0; u < n; ++u)
    {
        if(selected[u] == 1) continue;

        if(threshold[u] > neighbors[u].size()) continue; // u will always be in state 0, this should never happen by default

        if(threshold[u] >= obj) continue;
            
        C = GreedySub_sub(neighbors, threshold, u, obj); // the optimal fixed point that includes u

        for(auto q : C) selected[q] = 1;

        if(obj > C.size())
        {
            obj = C.size();
            A_star = C;
        }

    }

    cout<<"Vertices to be set to state one are:"<<endl;
    cout<<"[";
    for(auto v : A_star) cout<<v<<", ";
    cout<<"]"<<endl;
    cout<<"----------\n";

    // ----------------------------------------------
    // -    Check if it is indeed a fixed point     -
    // ----------------------------------------------
    vector<int> state_one(n, 0);
    for(auto u : A_star) state_one[u] = 1;

    for(int u = 0; u < n; ++u)
    {
        int inp = 0;
        for(auto v : neighbors[u])
        {
            inp += state_one[v];
        }
        if(state_one[u] == 1 && inp < threshold[u])
        {
            cout<<"not a fixed point"<<endl;       
        }
        if(state_one[u] == 0 && inp >= threshold[u])
        {
            cout<<"not a fixed point"<<endl;            
        }
    }

    return A_star; 
}
     

