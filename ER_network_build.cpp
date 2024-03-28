#include<bits/stdc++.h>
#include<chrono>

using namespace std;
using namespace std::chrono;

std::random_device r;
std::seed_seq seed{r(), r(), r(), r(), r()};
std::mt19937 gen(seed);


int N0;//final number of nodes

int a;
int b;


vector <int> NodeA; //node at one side of a bond, it will contain M elements, M= number of bonds
vector <int> NodeB; //node at another side of a bond, it will contain M elements, M= number of bonds

vector <int> node1;
//vector <int> node2;
set<pair<int, int>> unique_bonds;


void seed_generator() {
    auto now = system_clock::now();
    auto duration = now.time_since_epoch();
    auto seed_value = duration.count();
    gen.seed(seed_value);
}

void network_build()
{   seed_generator();
	/*NodeA.reserve(N0);
	NodeB.reserve(N0);
	node1.reserve(N0);
	node2.reserve(N0);*/
	
	for(int i=0; i<N0; i++){
		node1.push_back(i);
		//node2.push_back(i);
	}
	shuffle(node1.begin(),node1.end(), gen);
	//shuffle(node2.begin(),node2.end(), gen);
	
	while (unique_bonds.size() < N0) {
        uniform_int_distribution<int> dist1(0, N0 - 1);
        uniform_int_distribution<int> dist2(0, N0 - 1);

         a = dist1(gen);
         b = dist2(gen);

        if (node1[a] != node1[b]) {
            unique_bonds.insert({min(node1[a], node1[b]), max(node1[a], node1[b])});
        }
    }

    for (const auto& bond : unique_bonds) {
        NodeA.push_back(bond.first);
        NodeB.push_back(bond.second);
    }
   
}

int main() {
    N0 = 20; // You can set your desired value for N0
    network_build();

    // Print the generated network for verification
    for (int i = 0; i < N0; i++) {
        cout << NodeA[i] << " " << NodeB[i] << endl;
    }

    return 0;
}
