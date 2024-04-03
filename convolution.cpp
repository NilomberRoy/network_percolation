// compile with : g++ -fopenmp convolute_omp.cpp -o convolute
#include<bits/stdc++.h>
#include<omp.h>

using namespace std;

int convolutenum=0,number_of_threads=2,networksize=319201,m=100,ensemblesize=8000;

double threshold=1e-21;

string file_identity,input_file_string,output_file_string;

void define_input_output_file()

{

	//input_file_string=file_identity+"MDA_"+to_string(networksize)+"node_"+to_string(m)+"m_"+"1.9M_"+to_string(ensemblesize)+"ens_"+to_string(convolutenum)+"convoluted_PR"+".dat";

	//output_file_string=file_identity+"MDA_"+to_string(networksize)+"node_"+to_string(m)+"m_"+"1.9M_"+to_string(ensemblesize)+"ens_"+to_string(convolutenum+1)+"convoluted_PR"+".dat";

//iep_specheat_300k_1convoluted.dat


//input_file_string=file_identity+"_"+to_string(networksize)+"node_"+to_string(convolutenum)+"convoluted.dat";

//output_file_string=file_identity+"_"+to_string(networksize)+"node_"+to_string(convolutenum+1)+"convoluted.dat";

input_file_string = file_identity+".dat";
output_file_string= file_identity+ "conv.dat";

}

void convolute()
{
	ifstream fin;
	ofstream fout;

	/*string identity;

	if (file_identity==1)

	identity="delmaxMDA";

	if( file_identity==2)

	identity="entropyMDA";

	if(file_identity==3)

	identity="maxclusMDA";*/




	int count = 0;
	double a=0,b=0;
	// input data file
	/*string input_file_string1 =to_string(networksize)+"node_"+to_string(m)+"m_"+to_string(M)+"M_"+to_string(ensemblesize)+"ens"+".dat";; //change this

	string input_file_string= file_identity+ input_file_string1;*/


	//int Node = 800000; //change this its connection number

	int Node=networksize; //Network Size=400000
	int node1 =Node;	 //500000;//	2000000;//1500000;//1000000;	//750000;	//change this its node number

	//string output_file_string = "convolute_omp_"+input_file_string;

	/*string output_file_string1= to_string(networksize)+"node_"+to_string(m)+"m_"+to_string(M)+"M_"+to_string(ensemblesize)+"ens"+to_string(convolutenum)+".dat";*/


	/*temporary_file_string=input_file_string;

	cout<<"temporary_file:"<<temporary_file_string<<endl;

	temporary_file_string.erase(temporary_file_string.begin()+39,temporary_file_string.end());

	cout<<"temporary_file:"<<temporary_file_string<<endl;

	output_file_string=temporary_file_string+to_string(convolutenum+1)+"convoluted"+".dat";*/

	cout<<"input_file:"<<input_file_string<<endl;

	cout<<"output_file:"<<output_file_string<<endl;


//sheraj has changed this.

	fin.open(input_file_string.c_str());

	while(fin) //to know the size of array to initialize
	{
		fin>>a>>b;
		count++;

		if(count==(networksize+1))

		{

			cout<<"count:"<<count<<endl;

			break;

		}

		//cout<<count<<endl;
	}



	fin.close();
	cout<<"convolution for total connection : "<<count-1<<endl;

// std::vector can be used here
	double* c = NULL;
	c = new double[count+1];
	double* d = NULL;
	d = new double[count+1];

	for(int i=0;i<count+1;i++)
	{
		c[i] = 0;
		d[i] = 0;
	}

	fin.open(input_file_string.c_str());

	int n = 1;
	cout << "Reading from file ...";
	while(fin) //takes in data
	{
		fin>>a>>b;
		c[n] = a;
		d[n] = b;
		n++;
		if(n==(networksize+2))

		{

			cout<<"n:"<<n<<endl;

			break;

		}

	}

	fin.close();

	c[networksize+1]=0;
	d[networksize+1]=0;

	c[networksize+2]=0;
	d[networksize+2]=0;


	cout << "completed" << endl;
	/*double* binom = NULL;
	binom = new double[count];
	*/
	double* summ = NULL;
	summ = new double[count+1];

	for(int i=0;i<count+1;i++) summ[i] = 0;

	//for(int i=0;i<count;i++) binom[i] = 0;

	//fout.open(output_file_string.c_str());

	/*for(int j=1;j<count;++j) //main convolution of data takes place
	{
		int x,y;
		double prob = j*1.0/(count-1);
		binom[j] = 1;

		for(int i=j+1;i<count;++i)
		binom[i] = binom[i-1]*((count-1)-i+1)*1.0/i*prob/(1-prob);

		for(int i=j-1;i>=0;--i)
		binom[i] = binom[i+1]*(i+1)*1.0/((count-1)-i)*(1-prob)/prob;

		double sum = 0;

		for(int i=0;i<count;++i) sum += binom[i];
		for(int i=0;i<count;++i) binom[i] /= sum;

		sum = 0;

		for(int i=1;i<count;++i) sum += d[i]*binom[i];

		summ[j] = sum;

		//fout<<c[j]<<'\t'<<sum<<endl;
		if(j%1000==0)
		cout<<count-1<<'\t'<<j<<endl;
		//printf("%12d%10.5f%18.8e%18.8e", j, prob, d[j], sum);
		//cout<<endl;
	}*/

	double *factor1 = new double[count];
	double *factor2 = new double[count];

 	cout << "Generating factors...";
	for (int i=0;i<count;++i)
	{
	    factor1[i] = (double) (Node-i+1) / i;
	    factor2[i] = (double) (i+1) / (Node-i);
	}
	cout << "completed" << endl;
	//fout.open(output_file_string.c_str());
	//here the main convolution of data takes place
	cout << "Entering parallel region ..." << endl;
	size_t temp00 = Node / 1000;

	//int number_of_threads=omp_get_num_threads() ;//Sheraj has changed it;in the main code in the RHS it was thread & no semi colon after the line.

	//int number_of_thread = omp_get_num_thread(); //Sheraj has changed it;
	cout << "Total " << omp_get_num_threads() << " threads in sequential region" << endl;

omp_set_dynamic(0);     // Explicitly disable dynamic teams
omp_set_num_threads(number_of_threads);

#pragma omp parallel for



//cout<<"Total " << number_of_threads << " threads_in_parallel_region" << endl;

	for (int j=1;j<=Node;++j) //start from 1
	{

		if(j==1)

		cout<<"Total "<<omp_get_num_threads()<<"threads in parallel region"<<endl;

		if(j % temp00 == 0){
			//int tid = omp_get_thread_num();//Shahnoor's vai's line

			int tid = omp_get_num_threads();
			//cout<<"Total " << tid << " threads in parallel region" << endl;
			//int tid = omp_get_num_threads();//Sheraj have changed it
			//cout << "iteration " << j << " of " << Node << " thread " << tid << endl;
		}
	    double prob     = (double) j / Node;
	    double factor   = 0;
	    double binom    = 0;
	    double prev     = 0;
	    double binomTot = 1;
	    double sum      = d[j];

	    factor = prob / (1-prob);
	    prev   = 1;

	    for (int i=j+1;i<=Node;++i)
		 {
	        binom     = prev * factor1[i] * factor;
	        binomTot += binom;
	        sum      += d[i] * binom;
	        prev      = binom;
	        if(binom < threshold){break;}
	    }

	    factor = (1-prob)/prob;
	    prev   = 1;

	    for (int i=j-1;i>0;--i)
		 {
	        binom     = prev * factor2[i] * factor;
	        binomTot += binom;
	        sum      += d[i] * binom;
	        prev      = binom;
	        if(binom < threshold){break;}
	    }

	    summ[j] = sum / binomTot;
		 //fout<<c[j]<<'\t'<<summ[j]<<endl;
	   // if(j%1000==0)// just to make the process faster and check progress
	  //  cout<<count-1<<'\t'<<j<<endl;
	}


	//cout<<omp_get_num_threads();

	//int tid = omp_get_num_threads();

	//cout<<"Total " << tid << " threads_in_parallel_region" << endl;

	cout << "completed" << endl;
	fout.open(output_file_string.c_str());

	cout<<"writing to file"<<endl;

	for(int i=1;i<=Node;i++){
		fout<<fixed<<setprecision(20)<<(double)i/networksize<<'\t'<<setprecision(40)<<summ[i]<<endl;
	}

	//input_file_string=output_file_string;

	//cout<<"input_file"<<input_file_string<<endl;

	cout<<"done!"<<endl;

	fout.close();
	delete [] factor1;
	delete [] factor2;

	//cout<<"writing to file "<<endl;

	//fout.open(output_file_string.c_str());

	//for(int i=1;i<count;i++) fout<<c[i]<<'\t'<<summ[i]<<endl;



	delete [] c;
	c = NULL;
	delete [] d;
	d = NULL;

	/*delete [] binom;
	binom = NULL;
	*/
	delete [] summ;
	summ = NULL;




}

int main(int argc,char* argv[])
{

	//string input_file_string="maxclusMDA_600000node_50m_2M_30000ens.dat";

	//int networksize=600000;

	/*for(int i=0;i<argc;i++)

	cout<<argv[i]<<endl;

	cout<<"argc"<<argc<<endl;*/

	//string file_identity=argv[3]

	file_identity=argv[1];

	number_of_threads=atoi(argv[2]);

	convolutenum=atoi(argv[3]);

	networksize=atoi(argv[4]);

/*
	m=atoi(argv[5]);

	M=atoi(argv[6]);

	ensemblesize=atoi(argv[7]);

*/

	/*cout<<"number_of_threads"<<number_of_threads<<endl;

	cout<<"networksize"<<networksize<<endl;

	cout<<"convolutenum"<<convolutenum<<endl;

	cout<<"ok"<<endl;*/

	/*int number_of_threads=2;

	int m=50;

	int M=2;

	int ensemblesize=30000;*/
	//int convolutenum=1;

	define_input_output_file();

	convolute();

	return 0;
}
