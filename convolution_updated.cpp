#include<bits/stdc++.h>
using namespace std;

int networksize = 500000;
double threshold = 1e-21;

int main()
{
	ifstream fin ("filename.txt");
	ofstream fout ("filename_conv.txt");

	int i;
	double a1=0, a2=0, a3=0, a4=0, a5=0, a6=0, a7=0, a8=0, a9=0;
	vector <double> x(networksize+2, 0), y(networksize+2, 0), summ(networksize+2, 0);


	for (i=0; i<networksize; i++)
    {
        fin >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >> a8 >> a9;
        x[i] = a2;
        y[i] = a7;
    }


// main loop for convolution
	for (int j=1; j<=networksize; j++) //start from 1
	{
	    double prob = (double) j / networksize;
	    //double prob = x[j];
	    double binom = 0, prev=0, binomTot = 1, sum = y[j];

	    prev = 1;
	    for (int i=j+1; i<=networksize; i++) //loop for the values on the right of the selected cell
		 {
	        binom     = prev * ( (double)(networksize-i+1) / i ) * ( prob / (1-prob) );
	        binomTot += binom;
	        sum      += y[i] * binom;
	        prev      = binom;
	        if (binom < threshold)
                {   break;  }
	    }

	    prev = 1;
	    for (int i=j-1; i>0; i--)  //loop for the values on the left of the selected cell
		 {
	        binom     = prev * ((double) (i+1) / (networksize-i) ) * ( (1-prob)/prob );
	        binomTot += binom;
	        sum      += y[i] * binom;
	        prev      = binom;
	        if (binom < threshold)
                {   break;  }
	    }

	    summ[j] = sum / binomTot;
	}


// output loop
    for (i=0; i<=networksize; i++)
    {
        fout << (double)i/networksize  << "    " << x[i] << "   "  << y[i] << "   " << summ[i]<< endl;
       // fout << c[i] << "   "  << d[i] << "   " << summ[i]<< endl;
    }

	return 0;

}
