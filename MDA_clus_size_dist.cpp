#include <bits/stdc++.h>
#include <chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r(), r(), r(), r(), r()};
std::mt19937 gen(seed);

int m0 = 0;
int N  = 0;
int filenum = 0;
int M  = 0;
int ensemble_count = 0;
int m  = 0;

double T = 0;

vector<int> ptr;
vector<int> created_bonds;
vector<int> shuffled_bonds;

vector<vector<int>> Nodes_Neigh;
vector<int> NodeA;
vector<int> NodeB;

// Cluster size data
vector<int> cluster_sizes; //vector for counting each cluster sizes

// Output file
string dist_file;

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
    if (ptr[i] < 0) return i;
    return ptr[i] = findroot(ptr[i]);
}

void update_root(int father, int child)
{
    ptr[father] += ptr[child];
    ptr[child] = father;
}

void ensemble_initialization()
{
    for (int a = 0; a < NodeA.size(); a++)
    {
        created_bonds.push_back(a);
    }

    for (int h = 0; h < created_bonds.size(); h++)
    {
        shuffled_bonds.push_back(created_bonds[h]);
    }

    shuffle(shuffled_bonds.begin(), shuffled_bonds.end(), gen);
    vector<int>().swap(created_bonds);
}

int bondselection(int A)
{
    vector<pair<long double, int>> product_and_index;

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
    for (int ensemble = 0; ensemble < ensemble_count; ensemble++)
    {
        Nodes_Neigh.clear();
        vector<int>().swap(NodeA);
        vector<int>().swap(NodeB);

        seed_initialization();
        network_build();
        ensemble_initialization();

        for (int i = 0; i < N; i++)
        {
            ptr.push_back(-1);
        }

        double t = T * N;

        for (int a = 0; a < (int)t; a++)
        {
            int a1;

            if (a <= NodeA.size() - M)
                a1 = bondselection(a);
            else
                a1 = shuffled_bonds[a];

            int x = NodeA[a1];
            int y = NodeB[a1];

            int x1 = findroot(x);
            int y1 = findroot(y);

            if (ptr[x1] == -1 && ptr[y1] == -1)
            {
                update_root(x1, y1);
            }
            else if (x1 != y1)
            {
                if (-ptr[x1] >= -ptr[y1])
                    update_root(x1, y1);
                else
                    update_root(y1, x1);
            }
        }

        cluster_sizes.resize(N + 1);

        for (int a = 0; a < N; a++)
        {
            if (ptr[a] < 0)
            {
                int size = -ptr[a]; 
                cluster_sizes[size] += 1;
            }
        }

        vector<int>().swap(shuffled_bonds);
        vector<int>().swap(ptr);

        cout << "Ensemble no. " << ensemble << " is done" << endl;
    }

    ofstream file05(dist_file);

    for (int i = 0; i < N + 1; i++)
    {
        if (cluster_sizes[i] > 0)
        {
            double avg_clus_sizes =
                (long double)cluster_sizes[i] / ensemble_count;

            file05 << fixed << setprecision(6)
                   << i << "  "
                   << setprecision(10)
                   << (long double) avg_clus_sizes << endl;

            cout << fixed << setprecision(6)
                 << i << "  "
                 << setprecision(10)
                 << (long double) avg_clus_sizes << endl;
        }
    }

    file05.close();
    vector<int>().swap(cluster_sizes);
}

int main(int argc, char* argv[])
{
    M = atoi(argv[1]); //no. of bonds we are comparing
    N = atoi(argv[2]); //total no. of nodes
    m0 = atoi(argv[3]); //initial nodes
    ensemble_count = atoi(argv[4]); //total ensemble size
    m = atoi(argv[5]); //no. of incoming edges
    T = atof(argv[6]); // fraction of bonds we are adding while percolation (time)
    filenum = atoi(argv[7]); //file number

    dist_file =
        "cluster_sizes_M" + to_string(M) +
        "_N" + to_string(N) +
        "_m0" + to_string(m0) +
        "_en" + to_string(ensemble_count) +
        "_m" + to_string(m) +
        "_t" + to_string(T) +
        "_fnum" + to_string(filenum) + ".dat";

    percolation();
    return 0;
}

