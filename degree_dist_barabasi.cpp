#include<bits/stdc++.h>
#include<chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r()};
std::mt19937 gen(seed);



// variables for network building

int m0; //initial number of seeds, all connected to each other.
int N0;//total number of nodes
int en_count=0;
int filenum=0;
int num_edges;
int m=0;
int total_deg;

string deg_dist;

vector<double>probabilities;
vector<double> intervals;
vector < vector <int> > Nodes_Neigh; //collection of nodes, it will contain N elements, N= number of nodes


void seed_initialization()
{
	
	Nodes_Neigh.resize(N0);
	
	// Create seed nodes and connect them with each other
    for (int i=0; i<m0-1; i++) 
	{
        for (int j= i+1; j<m0; j++) 
		{
           Nodes_Neigh[i].push_back(j);
           Nodes_Neigh[j].push_back(i);
        } 
    }
	
}



void network_build()
{
	
    for (int i = m0; i <N0; i++) 
	{     
	
        for(int s=0; s<m; s++){
        	
		
    
    total_deg = 0;
    for (const vector<int>& row : Nodes_Neigh) {
        //total_deg += accumulate(row.begin(), row.end(), 0);
        total_deg+= row.size();
        
    }
    //cout<<"totdeg "<<total_deg<<endl;
    vector<double>().swap(probabilities);
	
	for (size_t j = 0; j <= i; j++) {
    probabilities.push_back((double)Nodes_Neigh[j].size() / total_deg);
    //cout<<"ok ";
    }

	
	
	double interval_end = 0.0;
    
   for(int i=0; i<probabilities.size()-1; i++){
    	intervals.push_back(interval_end);
    	interval_end += probabilities[i];
	}
	
    intervals.push_back(1.0);
    //for(int i=0; i< intervals.size()  ; i++){
    //	cout<<"intervals= "<<intervals[i]<<endl;
	//}  

    // Select a random node from the existing nodes      
	uniform_real_distribution<long double> distribution(0,1);
    long double random_number= distribution(gen);
        int selected_node;
    for (size_t j = 0; j < intervals.size() - 1; ++j) {
        if (random_number >= intervals[j] && random_number < intervals[j + 1]) {
            selected_node = j;
           // cout<<"selec "<<selected_node<<endl;
            break;
        }
    }
        
    Nodes_Neigh[i].push_back(selected_node);
    Nodes_Neigh[selected_node].push_back(i);
    
    vector<double>().swap(intervals);
      }
    
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
     num_edges = Nodes_Neigh[i].size();
     //cout<<"ok3"<<endl;
    // increment the corresponding index in the degree distribution vector
    degree_distribution[num_edges]++;
    //cout<<"ok4"<<endl;
}
   //vector< vector<int> >().swap(Nodes_Neigh);
   Nodes_Neigh.clear();
   
   cout<<"ensem "<<en<<" ok"<<endl;
}
   
//cout<<"ok2"<<endl;
for(int i=0; i<degree_distribution.size(); i++){
    //cout<<fixed<< setprecision(10)<<(i)<<"  "<<setprecision(10)<<(double)log(i)<<"  "<<setprecision(10)<<(double)log((double)(degree_distribution[i])/(N0*en_count))<<endl;
    file<<fixed<< setprecision(10)<<(i)<<"  "<<setprecision(10)<<(double)log(i)<<"  "
	//<<setprecision(10)<<(double)log(pow(i, 0.5))<<"  "
	//<<setprecision(10)<<(double)log((double)(degree_distribution[i])/(N0*pow(i, 0.5)*en_count))<<"  "
	<<setprecision(10)<<(double)log((double)(degree_distribution[i])/(N0*en_count))<<endl;
    //sum_P+=(double)degree_distribution[i]/(N0*en_count);
}
     //cout<<sum_P;
     vector<int>().swap(degree_distribution);
}

int main(int argc, char* argv[]){
	N0 =atoi(argv[1]);
    m0=atoi(argv[2]);
    m= atoi(argv[3]);
    en_count=atoi(argv[4]);
    
    filenum=atoi(argv[5]);
    
    
    deg_dist="deg_barabasi_N0"+to_string(N0)+"_m0"+to_string(m0)+"_m"+to_string(m)+"_en"+to_string(en_count)+"_fnum"+to_string(filenum)+".dat";
	degree_dist();
	
    return 0;
}


