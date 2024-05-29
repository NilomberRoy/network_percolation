#include<bits/stdc++.h>
#include<chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r()};
std::mt19937 gen(seed);



// variables for network building
//number of hands for each node
int m0; //initial number of seeds, all connected to each other.
int N0;//final number of nodes
// N =N0-m0
int en_count=0;
 //number of bonds to compare
int filenum=0;
int num_edges;
int m=0;
string deg_dist;


vector < vector <int> > Nodes_vect; //collection of nodes, it will contain N elements, N= number of nodes
vector <int> NodeA; //node at one side of a bond, it will contain M elements, M= number of bonds
vector <int> NodeB; //node at another side of a bond, it will contain M elements, M= number of bonds



void seed_initialization()
{
	
	Nodes_vect.resize(N0);
	
	// Create seed nodes and connect them with each other
    for (int i=0; i<m0; i++) 
	{
        for (int j= i+1; j<m0; j++) 
		{
           Nodes_vect[i].push_back(j);
            Nodes_vect[j].push_back(i);
        } 
    }
	
	//cout<<"seed ok"<<endl;
	
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

void network_build()
{
	
	// Randomly select a node from the existing nodes and connect it to m neighbours
    for (int i = m0; i <N0; i++) 
	{
        // Select a random node from the existing nodes
        uniform_int_distribution<int> distribution(0, i-1);
        int selected_node = distribution(gen);
        
        // Select neighbours of the selected node
        vector<int> neighbours = Nodes_vect[selected_node];
        
        
        // Connect the new node to the selected m neighbours
        
        
        
        for (int j = 0; j <m ; j++) {
        	uniform_int_distribution<int> dist(0+j,neighbours.size()-1);
        	 double random_node = dist(gen);
        	 Nodes_vect[i].push_back(neighbours[random_node]);
        	 Nodes_vect[neighbours[random_node]].push_back(i);
        	 NodeA.push_back(i);
        	 NodeB.push_back(neighbours[random_node]);
        	 
        	 
        	 swap(neighbours[j], neighbours[random_node]);
        	 
		}
	vector<int>().swap(neighbours);
    }
}


//double sum_P=0;

void degree_dist(){
    
    ofstream file(deg_dist);
//cout<<"ok1"<<endl;
vector<int> degree_distribution(N0, 0); // initialize vector to store degrees of all nodes to zero

for(int en=0; en<en_count; en++){
  seed_initialization();
  network_build();
// loop through all nodes
for (int i = 0; i < N0; i++) {
    // count the number of edges for the current node
    //cout<<"ok2"<<endl;
     num_edges = Nodes_vect[i].size();
     //cout<<"ok3"<<endl;
    // increment the corresponding index in the degree distribution vector
    degree_distribution[num_edges]++;
    //cout<<"ok4"<<endl;
}
   //vector< vector<int> >().swap(Nodes_vect);
   Nodes_vect.clear();
   vector<int>().swap(NodeA);
   vector<int>().swap(NodeB);
   cout<<"ensem "<<en<<" ok"<<endl;
}
   
//cout<<"ok2"<<endl;
for(int i=0; i<degree_distribution.size(); i++){
    cout<<fixed<< setprecision(10)<<(i)<<"  "<<setprecision(10)<<(double)log(i)<<"  "<<setprecision(10)<<(double)log((double)(degree_distribution[i])/(N0*en_count))<<endl;
    file<<fixed<< setprecision(10)<<(i)<<"  "<<setprecision(10)<<(double)log(i)<<"  "<<setprecision(10)<<(double)log((double)(degree_distribution[i])/(N0*en_count))<<endl;
    //sum_P+=(double)degree_distribution[i]/(N0*en_count);
}
     //cout<<sum_P;
}

int main(int argc, char* argv[]){
	N0 =atoi(argv[1]);
    m0=atoi(argv[2]);
    en_count=atoi(argv[3]);
    m= atoi(argv[4]);
    filenum=atoi(argv[5]);
    
    
    deg_dist="mda_deg_N0"+to_string(N0)+"_m0"+to_string(m0)+"_en"+to_string(en_count)+"_p"+to_string(m)+"_fnum"+to_string(filenum)+".dat";
	degree_dist();
	
    return 0;
}


