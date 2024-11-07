#include<bits/stdc++.h>
#include<chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r(), r(), r(), r(), r()};
std::mt19937 gen(seed);

int m0 = 0; // initial number of seeds, all connected to each other.
int N = 0; // total number of nodes
int filenum = 0;
int M = 0; // number of bonds to compare
int ensemble_count = 0;
int m = 0; // number of randomly connected links
double alpha = 0; 
double beta = 0;
double t_Nalpha = 0;


vector<int> ptr; 
vector<int> created_bonds; 
vector<int> shuffled_bonds;

vector<int>nodes_of_cluster_sizes; //nodes_of_cluster_sizes[i] represents the size of the cluster.
//vector<int> largest_cluster;
//vector<long double> tempentropy;
vector<vector<int>> Nodes_Neigh; // collection of nodes with its neighbors
vector<int> NodeA; // node at one side of a bond
vector<int> NodeB; // node at another side of a bond

string FvsN;

void seed_generator() {
    auto now = system_clock::now();
    auto duration = now.time_since_epoch();
    auto seed_value = duration.count();
    gen.seed(seed_value);
}

void seed_initialization()
{
	
	Nodes_Neigh.resize(N);
	
	// Create seed nodes and connect them with each other
    for (int i=0; i<m0; i++) 
	{
        for (int j= i+1; j<m0; j++) 
		{
           Nodes_Neigh[i].push_back(j);
           Nodes_Neigh[j].push_back(i);
        } 
    }
	//seed bonds initialization
	for(int i=0; i<m0-1; i++)
	{
		for(int j=i+1; j<m0; j++)
		{
			NodeA.push_back(i);
			NodeB.push_back(j);
		}
	}
}

void network_build() //mda
{
	
	// Randomly select a node from the existing nodes and connect it to m neighbours
    for (int i = m0; i <N; i++) 
	{   seed_generator();
        // Select a random node from the existing nodes
        uniform_int_distribution<int> distribution(0, i-1);
        int selected_node = distribution(gen);
        
        // Select neighbours of the selected node
        vector<int> neighbours = Nodes_Neigh[selected_node];
        
        
        // Connect the new node to the selected m neighbours
         for (int j = 0; j <m ; j++) {
        	uniform_int_distribution<int> dist(0+j,neighbours.size()-1);
        	double random_node = dist(gen);
        	Nodes_Neigh[i].push_back(neighbours[random_node]);
        	Nodes_Neigh[neighbours[random_node]].push_back(i);
        	NodeA.push_back(i);
        	NodeB.push_back(neighbours[random_node]);
        	 
        	 
        	swap(neighbours[j], neighbours[random_node]);//for not conneting it in next iteration
        }
        vector<int>().swap(neighbours);
    }
}

int findroot(int i) {
    if (ptr[i] < 0) return i;
    else return ptr[i] = findroot(ptr[i]);
}

void merge_clusters(int father, int child) { 
    //cout << "Merging clusters: " << father << " and " << child << endl;
    //cout << "Before merge: father size = " << nodes_of_cluster_sizes[father] << ", child size = " << nodes_of_cluster_sizes[child] << endl;
    
    // Merge sizes at the current step
    nodes_of_cluster_sizes[father] += nodes_of_cluster_sizes[child];
    nodes_of_cluster_sizes[child] = 0; // Set absorbed cluster size to 0
    ptr[father] += ptr[child];
    ptr[child] = father;
    //cout << "After merge: father size = " << nodes_of_cluster_sizes[father] << endl;
    //cout<<step<<" no. iteration done"<<endl;
}



void percolation_initialization() {
    
    //largest_cluster.resize(N + 1, 0);
    nodes_of_cluster_sizes.resize(N + 1, 1);  
    //largest_cluster[0] = 1;
}

void ensemble_initialization() {
    for (int a = 0; a < NodeA.size(); a++) {
        created_bonds.push_back(a);
    }

    for (int h = 0; h < created_bonds.size(); h++) {
        shuffled_bonds.push_back(created_bonds[h]);
    }
    shuffle(shuffled_bonds.begin(), shuffled_bonds.end(), gen);

    /*for (int i = 0; i < N + 1; i++) {
        tempentropy.push_back(0.0);
    }
    tempentropy[0] = (long double) logl(N);*/

    vector<int>().swap(created_bonds);
}

int bondselection(int A) {
    vector<pair<long double, int>> product_and_index;

    // Calculating the product of cluster sizes for each bond to compare
    for (int i = 0; i < M; i++) {
        int e1 = shuffled_bonds[A + i];
        int s1 = NodeA[e1];
        int s2 = NodeB[e1];
        int s1_root = findroot(s1);
        int s2_root = findroot(s2);
        int s1_size = -ptr[s1_root];
        int s2_size = -ptr[s2_root];
        long double product = (long double) s1_size * (long double) s2_size;
        product_and_index.emplace_back(product, A + i);
    }

    // Sorting the bonds 
    sort(product_and_index.begin(), product_and_index.end());

    // Storing the original selection
    int original_selection = shuffled_bonds[product_and_index[0].second];

    // Swapping the selected bond 
    swap(shuffled_bonds[A], shuffled_bonds[product_and_index[0].second]);


    uniform_int_distribution<int> dist(A + 1, shuffled_bonds.size() - 1);

    
    for (int k = 1; k < M; k++) {
        int f = dist(gen);
        swap(shuffled_bonds[A + k], shuffled_bonds[f]);
    }
     product_and_index.clear();
     return original_selection;
}



void percolation() {
    ofstream file1(FvsN);
    //cout << "Starting percolation with ensemble_count = " << ensemble_count << endl; // Debug
    
    int total_node_count_in_large_clusters = 0; // To accumulate counts across ensembles
    

    for (int ensemble = 0; ensemble < ensemble_count; ensemble++) {
        //cout << "Processing ensemble no. " << ensemble << std::flush; // Added flushing
        
        Nodes_Neigh.clear();
        vector<int>().swap(NodeA);
        vector<int>().swap(NodeB);
        percolation_initialization();
        //cout<<"ok0"<<endl;
        seed_initialization();
        network_build();
        //cout<<"ok1"<<endl;
        ensemble_initialization();

        for (int i = 0; i < N; i++) {
            ptr.push_back(-1);
        }
            //cout<<"ok2"<<endl;
        int current_big = 1;
        //int edges_added = 0;



        for (int a = 0; a < (int)t_Nalpha; a++) {
           
            //edges_added++;
            
            // Simulating bonds and cluster updates
            int x, y, a1;

            if (a <= NodeA.size() - M) {
                a1 = bondselection(a);
            } else {
                a1 = shuffled_bonds[a];
            }

            //cout<<a<<"    "<<a1<<endl;

            x = NodeA[a1];
            y = NodeB[a1];

            int x1 = findroot(x);
            int y1 = findroot(y);

           
            // Copy cluster sizes from the previous step
            //This line copies the entire 1D vector from step a - 1 to step a. Effectively, 
            //it updates the row of the 2D vector nodes_of_cluster_sizes for step a with the values from the previous step
            //nodes_of_cluster_sizes[a + 1] = nodes_of_cluster_sizes[a];

            if (ptr[x1] == -1 && ptr[y1] == -1) {
                //update_root_and_entropy(a + 1, x1, y1);

                // Update cluster sizes
                merge_clusters(x1, y1);
                 //cout<<"ok3"<<endl;
            } 
            else if (x1 != y1) {

                if (-ptr[x1] >= -ptr[y1]) {

                    //update_root_and_entropy(a + 1, x1, y1);
                    // Update cluster sizes
                    merge_clusters(x1, y1);
                    //cout<<"ok4"<<endl;
                } else {

                    //update_root_and_entropy(a + 1, y1, x1);
                    merge_clusters(y1, x1);
                    //cout<<"ok5"<<endl;
                }
            } else {

                //tempentropy[a + 1] = tempentropy[a];
                //cout << "Merging clusters: " << x1 << " and " << y1 << endl;
                //cout << "Before merge: father size = " << nodes_of_cluster_sizes[a + 1][x1] << ", child size = " << nodes_of_cluster_sizes[a + 1][y1] << endl;
                //cout << "After merge: father size = " << nodes_of_cluster_sizes[a + 1][x1] << endl;
                //cout<<a+1<<" no. iteration done"<<endl;
                //cout<<"ok6"<<endl;
            }

            /*if (tempentropy[a + 1] < 0) {
                tempentropy[a + 1] = 0;
            }*/

            if (-ptr[findroot(x)] > current_big) {
                current_big = -ptr[findroot(x)];
            }

            //cout<<"ok7"<<endl;
            //largest_clustelr[a + 1] += current_big;
        }
        
        double N_beta = pow(N, 1 - beta);
    
        
int node_count_in_large_clusters = 0;

//if(edges_added <= t_Nalpha){
//cout << "N_beta: " << N_beta << endl;
//cout<<"nodes_of_cluster_sizes[(int)t_Nalpha - 1].size: "<<nodes_of_cluster_sizes[(int)t_Nalpha - 1].size()<<endl;
for (int i = 0; i < nodes_of_cluster_sizes.size() - 1; i++) {
            if (nodes_of_cluster_sizes[i] > (int)N_beta) {
                node_count_in_large_clusters += nodes_of_cluster_sizes[i];
                //cout << "Cluster no. of father node: " << i << ",  Size: " << nodes_of_cluster_sizes[i] << endl;

            }
            
        }
//cout<<" node_count_in_large_clusters: "<< node_count_in_large_clusters<<endl;


total_node_count_in_large_clusters += node_count_in_large_clusters;
        //cout<<"total_node_count_in_large_clusters:  "<<total_node_count_in_large_clusters<<endl; 

        
        vector<int>().swap(nodes_of_cluster_sizes);
        //nodes_of_cluster_sizes.clear();
        //nodes_of_cluster_sizes.shrink_to_fit(); // Optionally shrink memory allocation
        //vector<int>().swap(largest_cluster);
        vector<int>().swap(ptr);
        vector<int>().swap(shuffled_bonds);
        //vector<long double>().swap(tempentropy);

        //cout << "Ensemble no." << ensemble << " is done" << endl; // This should print
    }
    //cout << "Total node count in large clusters: " << total_node_count_in_large_clusters << endl;

    // Calculate the final average F value
    long double F = (long double) (total_node_count_in_large_clusters) / (N * ensemble_count);
    cout << N << "     " << F << endl;
    file1 << N << "     " << F << endl;
}


int main(int argc, char* argv[]){
	M=atoi(argv[1]);
	N =atoi(argv[2]);
    m0=atoi(argv[3]);
    ensemble_count=atoi(argv[4]);
    m=atoi(argv[5]);
    t_Nalpha=atof(argv[6]);
    beta=atof(argv[7]);
    filenum=atoi(argv[8]); 

    FvsN="FvsN_M"+to_string(M)+"_N"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_m"+to_string(m)+"_tNalpha"+to_string(t_Nalpha)+"_beta"+to_string(beta)+"_fnum"+to_string(filenum)+".dat";
    
	percolation();
	
    return 0;
}
