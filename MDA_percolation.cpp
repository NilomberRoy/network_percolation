#include<bits/stdc++.h>
#include<chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r(), r(), r(), r(), r()};
std::mt19937 gen(seed);

int m0=0; //initial number of seeds, all connected to each other.
int N=0 ;//total number of nodes
int filenum=0;
int M=0;//number of bonds to compare
int ensemble_count=0;
//double p=0.0;
int m=0; //number of randomly connected links
long double X;
double alpha = 0.0;

vector<int> ptr; 
vector<int>created_bonds; //rename of each bond
vector<int>shuffled_bonds; //shuffle of renamed bonds

vector<long double>big; //order parameter
vector<long int>largest_cluster;
vector<long double>X2;
vector<long double>X1;
vector<long double>X3;
vector<long double>X4;
vector<long double>tot_entropy; //shannon entropy
vector<long double>tempentropy;
vector<long double>U; //binder cumulant

vector < vector<int> > Nodes_Neigh; //collection of nodes with its neighbours, it will contain N elements, N= number of nodes
vector <int> NodeA; //node at one side of a bond, it will contain N elements, N= number of bonds
vector <int> NodeB; //node at another side of a bond, it will contain N elements, N= number of bonds
//const t=1/N;

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
/*In this function seed_initialization(), we created a complete graph or network, 
where every single node is connected to each other.
Total initial node number is m0.  */
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
/*In network_build() function we used a loop to build MDA network
by iterative method. Per unit time, new node arrives and it chose 
randomly a existing node called mediator. The new node will connect by bonds
with randomly chosen m number of neighbors of 
the mediator. */
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
/*In this function, we defind a method to find the root of each
cluster. This function arrives from the root finding algorithm.*/
int findroot(int i)
{
    if (ptr[i]<0) return i;
    else
    return ptr[i] = findroot(ptr[i]);
}
/*We defined Shannon entropy in entropy() function by the 
definition given by Prof. kamrul Hassan Mamun.*/
long double entropy(int i)
{
    int clustsize = -ptr[i];
          long double meu = (long double)clustsize/(long double)N;
          long double term = -1*meu*logl(meu);
           return term;
}
/*This function is being called when every bond is adding. The father 
node holds the cluster size. So the root of the the father node will 
be negative of its cluster size. ptr vector of father index will be merged. 
ptr vector of child index will be updated by the index of father node. 
From this method, we can easily tell the father index of a random node, 
which belongs to the same cluster. Tempentropy vector will update the 
entropy of each cluster. Initially entropy of each node is fixed. The 
entropy of the new cluster will be the defined by the size of the new
cluster. When cluster forms, the entropy of each new cluster will be 
calculated by subtracting the entropy of previously defined father node
and child node. Because after the change of cluster size, the individual 
entropy of child node and father node will not needed for the entropy of 
new cluster. */
void update_root_and_entropy(int bond,int father, int child)
{
			tempentropy[bond]= tempentropy[bond-1] - entropy(father) - entropy(child);

            ptr[father] += ptr[child];

            ptr[child]= father;
            
            tempentropy[bond] += entropy(father);
}
/*This function is used to initialize each and every vector's size
and their inital values for percolation. */
void percolation_initialization()
{
   tot_entropy.resize(N+1,0.0);
   big.resize(N+1,0.0);
   X1.resize(N+1,0.0);
   X2.resize(N+1,0.0);
   X3.resize(N+1,0.0);
   X4.resize(N+1,0.0);
   U.resize(N+1,0.0);
   largest_cluster.resize(N+1,0);
   big[0]+= (long double)1/N;
   
   
   largest_cluster[0]= 1;

   tot_entropy[0] += (long double)logl(N);

}
/*Same as percolation_initialization()*/
void ensemble_initialization()
{
	 for(int a=0; a<NodeA.size(); a++)
  {
      created_bonds.push_back(a);
  }
	
    for(int h=0; h<created_bonds.size(); h++)
    {
        shuffled_bonds.push_back(created_bonds[h]);
    }
	shuffle(shuffled_bonds.begin(),shuffled_bonds.end(), gen);

for (int i=0; i<N; i++)
	{

		tempentropy.push_back(0.0);
	}

	tempentropy[0] = (long double) logl(N);
	
vector<int>().swap(created_bonds);//clearing created bonds, cause we don't need them after shuffling.


//cout<<"Ensemble initialization ok"<<endl;
}
/*We're selecting one bond from M number of bonds. 
To select a bond, We're using product rule of Achlioptas Process. */
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
    int original_selection = shuffled_bonds[product_and_index[0].second];//selected bond is the minimum product of s1 and s2

    // Swapping the selected bond 
    swap(shuffled_bonds[A], shuffled_bonds[product_and_index[0].second]);// for not chose a bond twice


    uniform_int_distribution<int> dist(A + 1, shuffled_bonds.size() - 1);

    
    for (int k = 1; k < M; k++) {
        int f = dist(gen);
        swap(shuffled_bonds[A + k], shuffled_bonds[f]);
    }

    
    //shuffle(shuffled_bonds.begin() + A , shuffled_bonds.end(), gen);
     product_and_index.clear();
     return original_selection;
}


void percolation()
{   
   
    ofstream file1(ent);
    ofstream file2(largestclus);
    ofstream file3(specific_heat);
    ofstream file4(susceptibility);
    ofstream file5(u);
    ofstream file6(powder);
    ofstream file7(t_Nalpha);

    int explosivity_time_count = 0; 
    int time_needed_to_reach_Nalpha = 0;
    
    for (int ensemble = 0; ensemble < ensemble_count; ensemble++)
    {   
        Nodes_Neigh.clear();
        vector<int>().swap(NodeA);
        vector<int>().swap(NodeB);    
	    //cout<<"ok1"<<endl;
	    percolation_initialization();
	    
        seed_initialization();
        //cout<<"ok1.5"<<endl;
        network_build();
        //cout<<"ok2"<<endl;
       
        ensemble_initialization();
        //cout<<"ok3"<<endl;
        for (int i=0; i<N; i++)
	    {
		ptr.push_back(-1);
	    }
	    //cout<<"ok4"<<endl;
        /*This ptr vector is inputing initial root -1 for every node. 
        The root will change when nodes will form cluster of size greater than
        their previous size.*/

        double N_alpha = pow(N, alpha); 
        int current_big = 1;//initial largest cluster size
        //cout<<"ok5"<<endl;
        for (int a = 0; a < N; a++)
        {//percolation starts from here. we'll connect N number of bonds so that the loop will generate N number of turns.
            int x, y, a1;

            if (a <= NodeA.size() - M)
            {
                a1 = bondselection(a);//selected bond 
            }
            else
            {
                a1 = shuffled_bonds[a];//selected bond
            }
            
            //cout<<a1<<endl;
            
			x = NodeA[a1];//node id number of one corner of this bond
            y = NodeB[a1];//node id number of another corner of this bond

            int x1 = findroot(x);//root number of node id x  
            int y1 = findroot(y);//root number of node id y
            //cout<<"ok4"<<endl;
            if (ptr[x1] == -1 && ptr[y1] == -1)//if ptr is equal and nodes are isolated
            {
                update_root_and_entropy(a+1,x1,y1);
               // cout<<"ok7"<<endl;
            }
            
            else {
            	
			
			if (x1 != y1)//roots are different
            {   //cout<<"ok5"<<endl;
                if (-ptr[x1] >= -ptr[y1])
                {
                  update_root_and_entropy(a+1,x1,y1);
                  //father node will be x1 because x1 holds larger cluster than y1.
                }
                else
                {
                   
                  update_root_and_entropy(a+1,y1,x1);
                  //father node will be y1 because y1 holds larger cluster than x1.
                }
            }
            
            else {
			//new bond is connecting within the same cluster. roots are same. entropy will not change here.
            tempentropy[a+1]=tempentropy[a];
			}
           }
        
            if(tempentropy[a+1]<0)
		    {tempentropy[a+1]=0;} //not needed actually
             //cout<<"ok6"<<endl;
            
           
             if (-ptr[findroot(x)] > current_big)//comparing current largest cluster with the new cluster after adding bond.
            {
                current_big = -ptr[findroot(x)];
            }
            //if(-ptr[y1] > current_big){
            //	current_big = -ptr[y1];
			//}
            big[a+1] +=(long double)current_big/N;
            X =(long double) current_big/N;
            X1[0] += (long double)pow(1/N, 4);
            X2[0] += (long double)pow(1/N, 2);
            X2[a+1] +=(long double)pow(X,2);
            X1[a+1] +=(long double)pow(X,4);
            
			 //cout<<"ok7"<<endl;
            tot_entropy[a+1]+=(long double)tempentropy[a+1];//total entropy is being merged

            largest_cluster[a+1] = current_big;
            
		    if(largest_cluster[a+1] <= (int)N_alpha){
                time_needed_to_reach_Nalpha+=1;
            }	 
            

            if(largest_cluster[a+1] > (int)N_alpha && largest_cluster[a+1] < (int)(N/2)) //powder keg growth-time count 
            {
                explosivity_time_count+= 1;
                
            }    
        //cout<<"ok8"<<endl;
        }
     // Reset the network and bond configurations for each ensemble
        
        //vector< vector<int> >().swap(Nodes_Neigh);
        
        vector<long int>().swap(largest_cluster);

        vector<int>().swap(ptr);
      
        vector<int>().swap(shuffled_bonds);
        
        vector<long double>().swap(tempentropy);
        
        //cout<<"ensemble no."<< ensemble<<" "<<"is done"<<endl;
	}
	 
	 for(int i=0; i<N; i++){
	 	X3[i] = (long double) X1[i]/ensemble_count;
	 	X4[i] = (long double) X2[i]/ensemble_count;
	 	
	    }
   
    
    // Calculate average cluster sizes and print the results
    for (int i = 0; i < N; i++)
    {    
         double avg_entropy = (long double)(tot_entropy[i]) / (ensemble_count);
         double avg_cluster_size = (long double)(big[i]) / (ensemble_count) ;

         U[i] = (long double) (1.0 - (X3[i]/(3.0*pow(X4[i],2))));
         double avg_specific_heat = (long double) ((1- i/N)*(tot_entropy[i] - tot_entropy[i+1])*N)/(ensemble_count);;
         double avg_susceptibility = (long double)((big[i+1] - big[i])*N)/(ensemble_count);
         
        file4 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << avg_susceptibility<<endl;
    	file3 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << avg_specific_heat<<endl;
    	file5 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << U[i]<<endl;
        file2 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << avg_cluster_size<<endl;
        file1 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << avg_entropy << endl;
    }

    double time_needed_to_reach_Nalpha_avg = (long double) (time_needed_to_reach_Nalpha)/(ensemble_count);///for powder keg
    double explosivity_time_count_avg = (long double) (explosivity_time_count)/(ensemble_count); //explosivity time count
    file6 << N << "  " <<setprecision(10) << explosivity_time_count_avg << endl;
    cout << N << "  " <<setprecision(10) << explosivity_time_count_avg << endl;
    file7 << N << "  " <<setprecision(10) << time_needed_to_reach_Nalpha_avg << endl;
    cout << N << "  " <<setprecision(10) << time_needed_to_reach_Nalpha_avg << endl;
   
   vector<long double>().swap(tot_entropy);
   vector<long double>().swap(big);
   vector<long double>().swap(X1);
   vector<long double>().swap(X2);
   vector<long double>().swap(X3);
   vector<long double>().swap(X4);
   vector<long double>().swap(U);
}

int main(int argc, char* argv[]){
	M=atoi(argv[1]);
	N =atoi(argv[2]);
    m0=atoi(argv[3]);
    ensemble_count=atoi(argv[4]);
    m=atoi(argv[5]);
    alpha=atof(argv[6]);
    filenum=atoi(argv[7]);
    
    t_Nalpha="t_Nalpha_M"+to_string(M)+"_N"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_m"+to_string(m)+"_alpha"+to_string(alpha)+"_fnum"+to_string(filenum)+".dat";
    powder="powder_M"+to_string(M)+"_N"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_m"+to_string(m)+"_alpha"+to_string(alpha)+"_fnum"+to_string(filenum)+".dat";
    largestclus="op_M"+to_string(M)+"_N"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_m"+to_string(m)+"_alpha"+to_string(alpha)+"_fnum"+to_string(filenum)+".dat";
    ent="ent_M"+to_string(M)+"_N"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_m"+to_string(m)+"_alpha"+to_string(alpha)+"_fnum"+to_string(filenum)+".dat";
    specific_heat="spe_M"+to_string(M)+"_N"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_m"+to_string(m)+"_alpha"+to_string(alpha)+"_fnum"+to_string(filenum)+".dat";
    susceptibility="sus_M"+to_string(M)+"_N"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_m"+to_string(m)+"_alpha"+to_string(alpha)+"_fnum"+to_string(filenum)+".dat";
    u="U_M"+to_string(M)+"_N"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_m"+to_string(m)+"_alpha"+to_string(alpha)+"_fnum"+to_string(filenum)+".dat";
	percolation();
	
    return 0;
}


