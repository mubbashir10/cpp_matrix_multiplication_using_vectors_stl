#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "matrixLib.h"


using namespace std;

int main (int argc, char* argv[]) {

	//defining variables
	int n = 0;

	//asking for size
	cout<< "Enter size of matrix in power of 2 (e.g. 4 for 4x4): ";
	cin>>n;


    vector<int> inner (n);
    vector< vector<int> > A(n, inner), B(n, inner), C(n, inner);
    vector<vector<int> > cTraditional(n, vector<int>(n));

    fill (A, B, n);
    strassen(A, B, C, n);
    printMatrix(C, n);
    matrixTraditional(A,B,cTraditional,n);
    return 0;
}