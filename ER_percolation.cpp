#include<bits/stdc++.h>
#include<chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r(), r(), r(), r(), r()};
std::mt19937 gen(seed);

int N = 0; // Total number of nodes
int filenum = 0;
int M = 0; // Number of bonds to compare
int ensemble_count = 0;
double alpha = 0.0;

vector<int> ptr; // Tracks the root of each node
vector<long double> big; // Order parameter
vector<long int> largest_cluster;
vector<long double> X2;
vector<long double> X1;
vector<long double> X3;
vector<long double> X4;
vector<long double> tot_entropy; // Shannon entropy
vector<long double> tempentropy;
vector<long double> U; // Binder cumulant

string ent;
string largestclus;
string specific_heat;
string susceptibility;
string u;
string powder;
string t_Nalpha;

void seed_generator() {
    auto now = system_clock::now();
    auto duration = now.time_since_epoch();
    auto seed_value = duration.count();
    gen.seed(seed_value);
}

int findroot(int i) {
    if (ptr[i] < 0) return i;
    else return ptr[i] = findroot(ptr[i]);
}

long double entropy(int i) {
    int clustsize = -ptr[i];
    long double meu = (long double)clustsize / (long double)N;
    long double term = -1 * meu * logl(meu);
    return term;
}

void update_root_and_entropy(int bond, int father, int child) {
    tempentropy[bond] = tempentropy[bond - 1] - entropy(father) - entropy(child);
    ptr[father] += ptr[child];
    ptr[child] = father;
    tempentropy[bond] += entropy(father);
}

void percolation_initialization() {
    tot_entropy.resize(N + 1, 0.0);
    tempentropy.resize(N + 1, 0.0);
    big.resize(N + 1, 0.0);
    X1.resize(N + 1, 0.0);
    X2.resize(N + 1, 0.0);
    X3.resize(N + 1, 0.0);
    X4.resize(N + 1, 0.0);
    U.resize(N + 1, 0);
    largest_cluster.resize(N + 1, 0);
    big[0] += (long double)1 / N;
    largest_cluster[0] = 1;
    tot_entropy[0] += (long double)logl(N);
    tempentropy[0] = (long double) logl(N);
    
    
}

pair<int, int> pairselection(set<pair<int, int>>& used_pairs) {
    vector<pair<long double, pair<int, int>>> product_and_pair;

    for (int i = 0; i < M; i++) {
        // Randomly select two distinct nodes
        uniform_int_distribution<int> dist(0, N - 1);
        int s1 = dist(gen);
        int s2 = dist(gen);
        while (s1 == s2 || used_pairs.count({s1, s2}) || used_pairs.count({s2, s1})) { // Ensure s1 and s2 are distinct and pair is not used
            s1 = dist(gen);
            s2 = dist(gen);
        }

        int s1_root = findroot(s1);
        int s2_root = findroot(s2);
        int s1_size = -ptr[s1_root];
        int s2_size = -ptr[s2_root];
        long double product = (long double)s1_size * (long double)s2_size;
        product_and_pair.emplace_back(product, make_pair(s1, s2)); // Store product and pair
    }

    // Sort pairs based on product of cluster sizes
    sort(product_and_pair.begin(), product_and_pair.end());

    // Select the pair with the smallest product
    pair<int, int> selected_pair = product_and_pair[0].second;

    // Mark the selected pair as used
    used_pairs.insert(selected_pair);
    used_pairs.insert({selected_pair.second, selected_pair.first}); // Mark reverse pair as used

    product_and_pair.clear();
    return selected_pair;
}

void percolation() {
    ofstream file1(ent);
    ofstream file2(largestclus);
    ofstream file3(specific_heat);
    ofstream file4(susceptibility);
    ofstream file5(u);
    ofstream file6(powder);
    ofstream file7(t_Nalpha);

    // Debug: Check if files are opened successfully
    if (!file1.is_open()) cerr << "Error opening file: " << ent << endl;
    if (!file2.is_open()) cerr << "Error opening file: " << largestclus << endl;
    if (!file3.is_open()) cerr << "Error opening file: " << specific_heat << endl;
    if (!file4.is_open()) cerr << "Error opening file: " << susceptibility << endl;
    if (!file5.is_open()) cerr << "Error opening file: " << u << endl;
    if (!file6.is_open()) cerr << "Error opening file: " << powder << endl;
    if (!file7.is_open()) cerr << "Error opening file: " << t_Nalpha << endl;

    int explosivity_time_count = 0;
    int time_needed_to_reach_Nalpha = 0;

    for (int ensemble = 0; ensemble < ensemble_count; ensemble++) {
        cout << "Ensemble: " << ensemble << endl; // Debug: Print ensemble number
        percolation_initialization();
        ptr.resize(N , -1); // Initialize each node as its own root

        // Track used pairs to ensure no pair is selected twice
        set<pair<int, int>> used_pairs;

        double N_alpha = pow(N, alpha);
        int current_big = 1;

        for (int a = 0; a < N; a++) {
            int x, y;

            // Select a pair using the product rule
            pair<int, int> bond = pairselection(used_pairs);
            x = bond.first; // NodeA
            y = bond.second; // NodeB

            //cout << "Selected pair: (" << x << ", " << y << ")" << endl; // Debug: Print selected pair

            int x1 = findroot(x);
            int y1 = findroot(y);

            if (ptr[x1] == -1 && ptr[y1] == -1) {
                update_root_and_entropy(a + 1, x1, y1);
            } else if (x1 != y1) {
                if (-ptr[x1] >= -ptr[y1]) {
                    update_root_and_entropy(a + 1, x1, y1);
                } else {
                    update_root_and_entropy(a + 1, y1, x1);
                }
            } else {
                tempentropy[a + 1] = tempentropy[a];
            }

            if (tempentropy[a + 1] < 0) {
                tempentropy[a + 1] = 0;
            }

            if (-ptr[findroot(x)] > current_big) {
                current_big = -ptr[findroot(x)];
            }

            big[a + 1] += (long double)current_big / N;
            long double X = (long double)current_big / N;
            X1[0] += (long double)pow(1/N, 4);
            X2[0] += (long double)pow(1/N, 2);
            X2[a + 1] += (long double)pow(X, 2);
            X1[a + 1] += (long double)pow(X, 4);

            tot_entropy[a + 1] += (long double)tempentropy[a + 1];
            largest_cluster[a + 1] = current_big;

            if (largest_cluster[a + 1] <= (int)N_alpha) {
                time_needed_to_reach_Nalpha += 1;
            }

            if (largest_cluster[a + 1] > (int)N_alpha && largest_cluster[a + 1] < (int)(N / 2)) {
                explosivity_time_count += 1;
            }
        }

        vector<long int>().swap(largest_cluster);
        vector<int>().swap(ptr);
        vector<long double>().swap(tempentropy);
        used_pairs.clear();
    }

    for (int i = 0; i <= N; i++) {
        X3[i] = (long double)X1[i] / ensemble_count;
        X4[i] = (long double)X2[i] / ensemble_count;
    }

    for (int i = 0; i <= N; i++) {
        double avg_entropy = (long double)(tot_entropy[i]) / (ensemble_count);
        double avg_cluster_size = (long double)(big[i]) / (ensemble_count);

        U[i] = (long double)(1.0 - (X3[i] / (3.0 * pow(X4[i], 2))));
        double avg_specific_heat = (long double)((1 - i / N) * (tot_entropy[i] - tot_entropy[i + 1]) * N) / (ensemble_count);
        double avg_susceptibility = (long double)((big[i + 1] - big[i]) * N) / (ensemble_count);

        file4 << fixed << setprecision(10) << (long double)(i) / N << "  " << setprecision(10) << avg_susceptibility << endl;
        file3 << fixed << setprecision(10) << (long double)(i) / N << "  " << setprecision(10) << avg_specific_heat << endl;
        file5 << fixed << setprecision(10) << (long double)(i) / N << "  " << setprecision(10) << U[i] << endl;
        file2 << fixed << setprecision(10) << (long double)(i) / N << "  " << setprecision(10) << avg_cluster_size << endl;
        file1 << fixed << setprecision(10) << (long double)(i) / N << "  " << setprecision(10) << avg_entropy << endl;

        // Debug: Print output data
        //cout << "Step: " << i << ", Avg entropy: " << avg_entropy << ", Avg cluster size: " << avg_cluster_size << endl;
    }

    double time_needed_to_reach_Nalpha_avg = (long double)(time_needed_to_reach_Nalpha) / (ensemble_count);
    double explosivity_time_count_avg = (long double)(explosivity_time_count) / (ensemble_count);
    file6 << N << "  " << setprecision(10) << explosivity_time_count_avg << endl;
    cout << N << "  " << setprecision(10) << explosivity_time_count_avg << endl;
    file7 << N << "  " << setprecision(10) << time_needed_to_reach_Nalpha_avg << endl;
    cout << N << "  " << setprecision(10) << time_needed_to_reach_Nalpha_avg << endl;

    vector<long double>().swap(tot_entropy);
    vector<long double>().swap(big);
    vector<long double>().swap(X1);
    vector<long double>().swap(X2);
    vector<long double>().swap(X3);
    vector<long double>().swap(X4);
    vector<long double>().swap(U);
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " M N ensemble_count alpha filenum" << endl;
        return 1;
    }

    M = atoi(argv[1]);
    N = atoi(argv[2]);
    ensemble_count = atoi(argv[3]);
    alpha = atof(argv[4]);
    filenum = atoi(argv[5]);

    // Debug: Print input arguments
    cout << "M: " << M << ", N: " << N << ", ensemble_count: " << ensemble_count << ", alpha: " << alpha << ", filenum: " << filenum << endl;

    t_Nalpha = "ER_t_Nalpha_M" + to_string(M) + "_N" + to_string(N) + "_en" + to_string(ensemble_count) + "_alpha" + to_string(alpha) + "_fnum" + to_string(filenum) + ".dat";
    powder = "ER_powder_M" + to_string(M) + "_N" + to_string(N) + "_en" + to_string(ensemble_count) + "_alpha" + to_string(alpha) + "_fnum" + to_string(filenum) + ".dat";
    largestclus = "ER_op_M" + to_string(M) + "_N" + to_string(N) + "_en" + to_string(ensemble_count) + "_alpha" + to_string(alpha) + "_fnum" + to_string(filenum) + ".dat";
    ent = "ER_ent_M" + to_string(M) + "_N" + to_string(N) + "_en" + to_string(ensemble_count) + "_alpha" + to_string(alpha) + "_fnum" + to_string(filenum) + ".dat";
    specific_heat = "ER_spe_M" + to_string(M) + "_N" + to_string(N) + "_en" + to_string(ensemble_count) + "_alpha" + to_string(alpha) + "_fnum" + to_string(filenum) + ".dat";
    susceptibility = "ER_sus_M" + to_string(M) + "_N" + to_string(N) + "_en" + to_string(ensemble_count) + "_alpha" + to_string(alpha) + "_fnum" + to_string(filenum) + ".dat";
    u = "ER_U_M" + to_string(M) + "_N" + to_string(N) + "_en" + to_string(ensemble_count) + "_alpha" + to_string(alpha) + "_fnum" + to_string(filenum) + ".dat";

    percolation();

    return 0;
}