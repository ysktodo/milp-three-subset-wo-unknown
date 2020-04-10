#include"gurobi_c++.h"
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<string.h>
#include<bitset>
#include<algorithm>
#include<map>
#include<iomanip>
#include<cmath>

using namespace std;

void practicalTestTrivium(void);
int trivium(int evalNumRounds, int threadNumber);

void practicalTestGrain128a(void);
int grain128a(int evalNumRounds, int threadNumber);
int grain128aSub(int evalNumRounds, int threadNumber);
