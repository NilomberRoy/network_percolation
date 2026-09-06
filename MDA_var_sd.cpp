#include <bits/stdc++.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r(), r(), r(), r(), r()};
std::mt19937 gen(seed);

int m0 = 0;
int N = 0;
int filenum = 0;
int M = 0;
int ensemble_count = 0;
int m = 0;

vector<int> ptr;
vector<int> created_bonds;
vector<int> shuffled_bonds;

vector<vector<int>> Nodes_Neigh;
vector<int> NodeA;
vector<int> NodeB;

// Output file
string var_file;
string sd_file;
string mclus_file;
string sdmean_file;

vector<long double> total_sum_S;
vector<long double> total_sum_S2;
vector<long double> total_sum_S3;

map<int, int> size_count;  


void seed_generator()
{
    auto now = system_clock::now();
    auto duration = now.time_since_epoch();
    auto seed_value = duration.count();
    gen.seed(seed_value);
}

void seed_initialization()
{
    Nodes_Neigh.resize(N);

    for (int i = 0; i < m0; i++)
    {
        for (int j = i + 1; j < m0; j++)
        {
            Nodes_Neigh[i].push_back(j);
            Nodes_Neigh[j].push_back(i);
        }
    }

    for (int i = 0; i < m0 - 1; i++)
    {
        for (int j = i + 1; j < m0; j++)
        {
            NodeA.push_back(i);
            NodeB.push_back(j);
        }
    }
}

void network_build()
{
    for (int i = m0; i < N; i++)
    {
        seed_generator();
        uniform_int_distribution<int> distribution(0, i - 1);
        int selected_node = distribution(gen);

        vector<int> neighbours = Nodes_Neigh[selected_node];

        for (int j = 0; j < m; j++)
        {
            uniform_int_distribution<int> dist(j, neighbours.size() - 1);
            int random_node = dist(gen);

            Nodes_Neigh[i].push_back(neighbours[random_node]);
            Nodes_Neigh[neighbours[random_node]].push_back(i);

            NodeA.push_back(i);
            NodeB.push_back(neighbours[random_node]);

            swap(neighbours[j], neighbours[random_node]);
        }

        vector<int>().swap(neighbours);
    }
}

int findroot(int i)
{
    if (ptr[i] < 0)
        return i;
    return ptr[i] = findroot(ptr[i]);
}

void update_root(int father, int child)
{
    ptr[father] += ptr[child];
    ptr[child] = father;
}

void percolation_initialization()
{
    total_sum_S.resize(N + 1, 0.0);
    total_sum_S2.resize(N + 1, 0.0);
    total_sum_S3.resize(N + 1, 0.0);

    // Initialize incremental tracking
    //size_count.clear();
    size_count[1] = N;  // N clusters of size 1
    
}

void ensemble_initialization()
{
    created_bonds.clear();
    for (int a = 0; a < (int)NodeA.size(); a++)
    {
        created_bonds.push_back(a);
    }

    shuffled_bonds.clear();
    for (int h = 0; h < (int)created_bonds.size(); h++)
    {
        shuffled_bonds.push_back(created_bonds[h]);
    }

    shuffle(shuffled_bonds.begin(), shuffled_bonds.end(), gen);
    vector<int>().swap(created_bonds);
}

int bondselection(int A)
{
    vector<pair<long double, int>> product_and_index;
    product_and_index.reserve(M);

    for (int i = 0; i < M; i++)
    {
        int e1 = shuffled_bonds[A + i];
        int s1 = NodeA[e1];
        int s2 = NodeB[e1];

        int s1_root = findroot(s1);
        int s2_root = findroot(s2);

        int s1_size = -ptr[s1_root];
        int s2_size = -ptr[s2_root];

        long double product =
            (long double)s1_size * (long double)s2_size;

        product_and_index.emplace_back(product, A + i);
    }

    sort(product_and_index.begin(), product_and_index.end());

    int original_selection =
        shuffled_bonds[product_and_index[0].second];

    swap(shuffled_bonds[A],
         shuffled_bonds[product_and_index[0].second]);

    uniform_int_distribution<int> dist(A + 1, shuffled_bonds.size() - 1);

    for (int k = 1; k < M; k++)
    {
        int f = dist(gen);
        swap(shuffled_bonds[A + k], shuffled_bonds[f]);
    }

    product_and_index.clear();
    return original_selection;
}

void percolation()
{
    ofstream file1(var_file);
    ofstream file2(sd_file);
    ofstream file3(mclus_file);
    ofstream file4(sdmean_file);

    for (int ensemble = 0; ensemble < ensemble_count; ensemble++)
    {   
        auto ens_start = high_resolution_clock::now(); 
        cout << "Ensemble " << ensemble + 1 << "/" << ensemble_count << " starting..." << endl;
        
        Nodes_Neigh.clear();
        vector<int>().swap(NodeA);
        vector<int>().swap(NodeB); 

        seed_initialization();
        network_build();
        percolation_initialization();
        ensemble_initialization();

        for (int i=0; i<N; i++){
		    ptr.push_back(-1);
	    }

       
        total_sum_S[0] += (long double)N;
        total_sum_S2[0] += (long double)N;
        total_sum_S3[0] += (long double)N;

        int largest_cluster = 1;

        
        for (int a = 0; a < N; a++)
        {   
            int a1;
            if (a <= (int)NodeA.size() - M)
                a1 = bondselection(a);
            else
                a1 = shuffled_bonds[a];

            int x = NodeA[a1];
            int y = NodeB[a1];
            int x1 = findroot(x);
            int y1 = findroot(y);

            

            if (x1 != y1)  // Only process if they are in different clusters
            {
                int size1 = -ptr[x1];
                int size2 = -ptr[y1];
                int new_size = size1 + size2;
                
                size_count[size1]--;
                if (size_count[size1] == 0)
                    size_count.erase(size1);
              
                size_count[size2]--;
                if (size_count[size2] == 0)
                    size_count.erase(size2);
                
                // Add the new merged cluster
                size_count[new_size]++;
                
                // Update largest cluster if needed
                if (new_size > largest_cluster)
                    largest_cluster = new_size;
                
                // Perform union (merge smaller into larger)
                if (size1 >= size2)
                    update_root(x1, y1);
                else
                    update_root(y1, x1);
            }
            
            
            long double sum_s_step = 0.0;
            long double sum_s2_step = 0.0;
            long double sum_s3_step = 0.0;
            bool skipped = false;
            
            // Only iterate over existing cluster sizes (typically far fewer than N)
            for (auto& entry : size_count)
            {
                int sz = entry.first;
                int cnt = entry.second;
                
                int use_cnt = cnt;
                if (!skipped && sz == largest_cluster)
                {
                    use_cnt = cnt - 1;  // Exclude ONE largest cluster (giant component)
                    skipped = true;
                }
                
                sum_s_step += (long double)sz * use_cnt;
                sum_s2_step += (long double)sz * sz * use_cnt;
                sum_s3_step += (long double)sz * sz * sz * use_cnt;
            }
            //cout<<"sum_s_step: "<<sum_s_step<<endl;
            //cout<<"sum_s2_step: "<<sum_s2_step<<endl;
            //cout<<"sum_s3_step: "<<sum_s3_step<<endl;
            // Store RAW sums (not ratios!) for this step
            
            long double S_mean = (long double) sum_s2_step / sum_s_step;
            long double S2_mean = (long double) sum_s3_step / sum_s_step;

            total_sum_S2[a + 1] += S_mean;
            total_sum_S3[a + 1] += S2_mean;
            
            // Progress indicator
        /*    if ((a + 1) % (N / 10) == 0 || a == N - 1)
            {
                cout << "  Step " << a + 1 << "/" << N << " completed" << endl;
            }*/
           
        }
        
        vector<int>().swap(ptr);
        vector<int>().swap(shuffled_bonds);
        map<int, int>().swap(size_count);
        
        auto ens_end = high_resolution_clock::now();
        auto ens_minutes = duration_cast<seconds>(ens_end - ens_start).count();
        cout << "Ensemble " << ensemble + 1 << " completed in " << ens_minutes << " seconds" << endl;
    }

    for (int i = 1; i < N; i++)
    {
        //if (total_sum_S[i] == 0.0 && i > 0) continue;

        long double avg_S = (long double) total_sum_S2[i] / ( ensemble_count);
        long double avg_S2 = (long double) total_sum_S3[i] / ( ensemble_count);
        
        long double variance = avg_S2 - avg_S * avg_S;
        
        
        long double std_dev = sqrt(variance);
        long double t = (long double) i / N;
        
        file1 << fixed << setprecision(10) << t << "  " << variance << endl;
        file2 << fixed << setprecision(10) << t << "  " << std_dev << endl;
        file3 << fixed << setprecision(10) << t << "  " << avg_S << endl;
        if (avg_S > 0) 
        {
            file4 << fixed << setprecision(10) << t << "  " << (long double) std_dev / avg_S << endl;
        }
        
        // Console output for final results
       /* if (i % (N/10) == 0 || i == N)
        {
            cout << "t = " << t << ", <S> = " << avg_S << ", Var = " << variance << endl;
        }*/
    }
    
    file1.close();
    file2.close();
    file3.close();
    file4.close();

    vector<long double>().swap(total_sum_S);
    vector<long double>().swap(total_sum_S2);
    vector<long double>().swap(total_sum_S3);
}

int main(int argc, char *argv[])
{
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    m0 = atoi(argv[3]);
    ensemble_count = atoi(argv[4]);
    m = atoi(argv[5]);
    filenum = atoi(argv[6]);

    var_file = "var_M" + to_string(M) +
               "_N" + to_string(N) +
               "_m0" + to_string(m0) +
               "_en" + to_string(ensemble_count) +
               "_m" + to_string(m) +
               "_fnum" + to_string(filenum) + ".dat";
    sd_file = "sd_M" + to_string(M) +
              "_N" + to_string(N) +
              "_m0" + to_string(m0) +
              "_en" + to_string(ensemble_count) +
              "_m" + to_string(m) +
              "_fnum" + to_string(filenum) + ".dat";
    mclus_file = "mclus_M" + to_string(M) +
                 "_N" + to_string(N) +
                 "_m0" + to_string(m0) +
                 "_en" + to_string(ensemble_count) +
                 "_m" + to_string(m) +
                 "_fnum" + to_string(filenum) + ".dat";
    sdmean_file = "sdmean_M" + to_string(M) +
                  "_N" + to_string(N) +
                  "_m0" + to_string(m0) +
                  "_en" + to_string(ensemble_count) +
                  "_m" + to_string(m) +
                  "_fnum" + to_string(filenum) + ".dat";

    percolation();
    
 
    return 0;
}