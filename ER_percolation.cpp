#include<bits/stdc++.h>
#include<chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r(), r(), r(), r(), r()};
std::mt19937 gen(seed);


int N=0 ;//final number of nodes
int filenum=0;
int M=0;//number of bonds to compare
int ensemble_count=0;
int a;
int b;

vector<int> ptr;
vector<int>bond_id;
vector<int>sbond_id; //shuffled bonds after numbering them

vector<long double>big;
vector<long double>X2;
vector<long double>X1;
vector<long double>X3;
vector<long double>X4;
vector<long double>tot_entropy;
vector<long double>tempentropy;
vector<long double>U; //binder cumulant


vector <int> NodeA; //node at one side of a bond, it will contain N elements, N= number of bonds
vector <int> NodeB; //node at another side of a bond, it will contain N elements, N= number of bonds

vector <int> nodes; //initializing nodes with putting identity
set<pair<int, int>> unique_bonds; //to store unique sets of pairs of nodes

string ent;
string largestclus;
string specific_heat;
string susceptibility;
string u;

void seed_generator() {
    auto now = system_clock::now();
    auto duration = now.time_since_epoch();
    auto seed_value = duration.count();
    gen.seed(seed_value);
}

void network_build()
{   
    seed_generator();

	/*NodeA.reserve(N);
	NodeB.reserve(N);
	node1.reserve(N);
	node2.reserve(N);*/
	
	for(int i=0; i<N; i++){
		nodes.push_back(i);
	}
	
	shuffle(nodes.begin(),nodes.end(), gen);
	
	while (unique_bonds.size() < N) {
        uniform_int_distribution<int> dist1(0, N - 1);
        uniform_int_distribution<int> dist2(0, N - 1);

         a = dist1(gen);
         b = dist2(gen);

        if (nodes[a] != nodes[b]) {
            unique_bonds.insert({min(nodes[a], nodes[b]), max(nodes[a], nodes[b])});
        }
    }

    for (const auto& bond : unique_bonds) {
        NodeA.push_back(bond.first);
        NodeB.push_back(bond.second);
    }
   
}


int findroot(int i)
{
    if (ptr[i]<0) return i;
    else
    return ptr[i] = findroot(ptr[i]);
}

long double entropy(int i)
{
    int clustsize = -ptr[i];
          long double meu = (long double)clustsize/(long double)N;
          long double term = -1*meu*logl(meu);
           return term;
}

void update_root_and_entropy(int bond,int father, int child)
{
			tempentropy[bond]= tempentropy[bond-1] - entropy(father) - entropy(child);

            ptr[father] += ptr[child];

            ptr[child]= father;
            
            tempentropy[bond] += entropy(father);
}

void percolation_initialization()
{  
   tot_entropy.resize(N,0.0);
   big.resize(N,0.0);	
   X1.resize(N,0.0);
   X2.resize(N,0.0);
   X3.resize(N,0.0);
   X4.resize(N,0.0);
   U.resize(N,0);
   
   big[0]+= (long double)1/N;
   tot_entropy[0] += (long double)logl(N);

}

void ensemble_initialization()
{
	 for(int a=0; a<NodeA.size(); a++)
  {
      bond_id.push_back(a);
  }
	
    for(int h=0; h<bond_id.size(); h++)
    {
        sbond_id.push_back(bond_id[h]);
    }
	shuffle(sbond_id.begin(),sbond_id.end(), gen);

for (int i=0; i<N; i++)
	{

		tempentropy.push_back(0.0);
	}

	tempentropy[0] = (long double) logl(N);
	
//cout<<"Ensemble initialization ok"<<endl;
}

//vector<long double>product_of_clustersize;
/*
int bondselection(int A)
{


   for(int i=0; i<M; i++)
   {
      int e1 = sbond_id[A+ i];
      int s1 = sNodeA[e1];
      int s2 = sNodeB[e1];
      int s1_root = findroot(s1);
      int s2_root = findroot(s2);
	  int s1_size= -ptr[s1_root];
      int s2_size= -ptr[s2_root];
	  long double product = (long double) s1_size* (long double) s2_size;
	 
	  
	  
      product_of_clustersize.push_back(product);
   }

  
  long double selected_bond_size = product_of_clustersize[0];   
   int selected_bond = A;
   for(int j=1; j<M; j++)
   {
       if(product_of_clustersize[j] < selected_bond_size)
       {
          selected_bond_size= product_of_clustersize[j];
          selected_bond= A+j;

       }
   }


    int original_selection = sbond_id[selected_bond];

    swap(sbond_id[A],sbond_id[selected_bond]);
    uniform_int_distribution<int> dist(A+1,sbond_id.size()-1);

    for(int k=1; k<M; k++)
        {
            int f = dist(gen);
            swap(sbond_id[A+k], sbond_id[f]);
        }

    vector<long double>().swap(product_of_clustersize);
	
	

    return original_selection;

}*/

int bondselection(int A) {
	
    vector<pair<long double, int>> product_and_index;

    // Calculate the product of cluster sizes for each bond to compare
    for (int i = 0; i < M; i++) {
        int e1 = sbond_id[A + i];
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
    int original_selection = sbond_id[product_and_index[0].second];

    // Swapping the selected bond 
    swap(sbond_id[A], sbond_id[product_and_index[0].second]);


    uniform_int_distribution<int> dist(A + 1, sbond_id.size() - 1);

    
    for (int k = 1; k < M; k++) {
        int f = dist(gen);
        swap(sbond_id[A + k], sbond_id[f]);
    }
     
     vector<pair<long double, int>>().swap(product_and_index);
     
     return original_selection;
}


void percolation()
{   
   
    ofstream file1(ent);
    ofstream file2(largestclus);
    ofstream file3(specific_heat);
    ofstream file4(susceptibility);
    ofstream file5(u);
    
    for (int ensemble = 0; ensemble < ensemble_count; ensemble++)
    {   
           
	    
	    percolation_initialization();
	    //cout<<"ok1"<<endl;
        
        network_build();
        //cout<<"ok2"<<endl;
       
        ensemble_initialization();
        //cout<<"ok3"<<endl;
        for (int i=0; i<N; i++)
	    {
		ptr.push_back(-1);
	    }
	    //cout<<"ok4"<<endl;
        
         
        int current_big = 1;
       // cout<<"ok5"<<endl;
        for (int a = 0; a < N; a++)
        {
            int x, y, a1;

            if (a <= NodeA.size() - M)
             {
                a1 = bondselection(a);
            }
            else
            {
                a1 = sbond_id[a];
            }
            
            //cout<<a1<<endl;
            
			x = NodeA[a1];
            y = NodeB[a1];

            int x1 = findroot(x);
            int y1 = findroot(y);
            //cout<<"ok4"<<endl;
            if (ptr[x1] == -1 && ptr[y1] == -1)
            {
                update_root_and_entropy(a+1,x1,y1);
               // cout<<"ok7"<<endl;
            }
            
            else {
            	
			
			if (x1 != y1)
            {   //cout<<"ok5"<<endl;
                if (-ptr[x1] >= -ptr[y1])
                {
                  update_root_and_entropy(a+1,x1,y1);
                }
                else
                {
                   
                  update_root_and_entropy(a+1,y1,x1);
                }
            }
            
            else {
			
            tempentropy[a+1]=tempentropy[a];
			}
        }
        
        if(tempentropy[a+1]<0)
		{tempentropy[a+1]=0;} //not needed actually
          // cout<<"ok6"<<endl;
            
           
             if (-ptr[findroot(x)] > current_big)
            {
                current_big = -ptr[findroot(x)];
            }
            //if(-ptr[y1] > current_big){
            //	current_big = -ptr[y1];
			//}
            big[a+1] +=(long double)current_big/N;
            
            X2[a+1] +=(long double)pow(X,2);
            X1[a+1] +=(long double)pow(X,4);
            
			 
            tot_entropy[a+1]+=(long double)tempentropy[a+1];
        
        }
     // Reset the network and bond configurations for each ensemble
        
       
        vector<int>().swap(ptr);
        
        vector<int>().swap(bond_id);
          
        vector<int>().swap(sbond_id);
        
        vector<long double>().swap(tempentropy);
        
        vector<int>().swap(NodeA);
        
        vector<int>().swap(NodeB);  
        
        vector<int>().swap(nodes); 
		 
        set<pair<int, int>>().swap(unique_bonds);
        
        cout<<"ensemble no."<<ensemble<<" "<<"is done"<<endl;
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
         //double avg_clus_pow2=(long double)(big[i]*big[i]*big[i]*big[i])/(ensemble_count*ensemble_count);
         //double avg_clus_pow4= (long double)(big[i]*big[i]*big[i]*big[i])/(ensemble_count);
         U[i] = (long double) (1.0 - (X3[i]/(3.0*pow(X4[i],2))));
         //double avg_specific_heat = (long double) (1- i/N)*((tot_entropy[i+1] - tot_entropy[i])*(-N));;
         //double avg_susceptibility = (long double) ((big[i+1] - big[i])* N);
         
        //file4 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << avg_susceptibility<<endl;
    	//file3 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << avg_specific_heat<<endl;
    	file5 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << U[i]<<endl;
        file2 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << avg_cluster_size<<endl;
        file1 <<fixed<<setprecision(10) << (long double)(i) / N << "  " <<setprecision(10) << avg_entropy << endl;
    }
    
    
   
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
    ensemble_count=atoi(argv[3]);
    filenum=atoi(argv[4]);
    
    largestclus="op_M"+to_string(M)+"_N"+to_string(N)+"_en"+to_string(ensemble_count)+"_fnum"+to_string(filenum)+".dat";
    ent="ent_M"+to_string(M)+"_N"+to_string(N)+"_en"+to_string(ensemble_count)+"_fnum"+to_string(filenum)+".dat";
    //specific_heat="spe_mm"+to_string(mm)+"_N0"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_prob"+to_string(m)+"_fnum"+to_string(filenum)+".dat";
    //susceptibility="sus_mm"+to_string(mm)+"_N0"+to_string(N)+"_m0"+to_string(m0)+"_en"+to_string(ensemble_count)+"_prob"+to_string(m)+"_fnum"+to_string(filenum)+".dat";
    u="U_M"+to_string(M)+"_N"+to_string(N)+"_en"+to_string(ensemble_count)+"_fnum"+to_string(filenum)+".dat";
	percolation();
	
    return 0;
}

