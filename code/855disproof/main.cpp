#include"gurobi_c++.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<string.h>
#include<bitset>
#include<algorithm>
#include<cmath>
#ifdef WIN32
#include <windows.h>
#else
#ifdef WIN64
#include<windows.h>
#else
#include <unistd.h>
#endif
#endif



using namespace std;
/*
The core function of Trivium:
x[i5]<-x[i3]*x[i4]+x[i2]+x[i1]+x[i5]
x[*]<-x[*] where * is in [0,288)/{i1,i2,i3,i4,i5}
@Para:
model: the MILP model describing the 3-subset division property  
x: the current k
i1...i5: the indices involved
*/
void triviumCoreThree(GRBModel& model, vector<GRBVar>& x, int i1, int i2, int i3, int i4, int i5) {
  GRBVar y1 = model.addVar(0, 1, 0, GRB_BINARY);//<-
  GRBVar y2 = model.addVar(0, 1, 0, GRB_BINARY);
  GRBVar y3 = model.addVar(0, 1, 0, GRB_BINARY);
  GRBVar y4 = model.addVar(0, 1, 0, GRB_BINARY);
  GRBVar y5 = model.addVar(0, 1, 0, GRB_BINARY);

  GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
  GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);

	// remove z3 and z4 because a = z3 = z4
  GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);

  model.addConstr(y1 <= x[i1]);
  model.addConstr(z1 <= x[i1]);
  model.addConstr(y1 + z1 >= x[i1]);

  model.addConstr(y2 <= x[i2]);
  model.addConstr(z2 <= x[i2]);
  model.addConstr(y2 + z2 >= x[i2]);

  model.addConstr(y3 <= x[i3]);
  model.addConstr(a <= x[i3]);
  model.addConstr(y3 + a >= x[i3]);

  model.addConstr(y4 <= x[i4]);
  model.addConstr(a <= x[i4]);
  model.addConstr(y4 + a >= x[i4]);

  model.addConstr(y5 == x[i5] + a + z1 + z2);

  x[i1] = y1;
  x[i2] = y2;
  x[i3] = y3;
  x[i4] = y4;
  x[i5] = y5;
}
/************************************************************
Given 1 monomial, it return the number of division trails, if the number is ODD, 
the monomial exist in the superpoly, otherwise, it do not exist. 
Since z=\sum s[66,93,162,177,243,288], we can evaluate the superpoly corresponding to each s[i] separately.
@Para
I: cube indices
J: monomial indices
evalNumRounds: the number of initialization rounds
target: 1-6 corresponding to s[66,93,162,177,243,288] resepectively. 
***************************************************************/
int veryfi855(vector<int> I, vector<int> J, int evalNumRounds, int target, int threadNumber) {
  ofstream outputfile;
  outputfile.open("log.txt", ios::app);
  outputfile << "++++++++++++++++++++++++++++" << endl;

  //gurobi
  try {
    // Create the environment
    GRBEnv env = GRBEnv();

    // close standard output
    env.set(GRB_IntParam_LogToConsole, 0);
    env.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_BESTBOUND);
    env.set(GRB_IntParam_Threads, threadNumber);
		env.set(GRB_StringParam_LogFile, "log.txt");
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, 2000000000);
		env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);

    // Create the model
    GRBModel model = GRBModel(env);

    // Create variables
    vector<vector<GRBVar>> s(evalNumRounds + 1, vector<GRBVar>(288));
    for (int i = 0; i < 288; i++) 
      s[0][i] = model.addVar(0, 1, 0, GRB_BINARY);
		GRBVar o = model.addVar(0, 1, 0, GRB_BINARY);

    // IV constraint
    for (int i = 0; i < I.size(); i++)
      model.addConstr(s[0][93 + I[i]] == 1);
    
    // Key constraint
    for (int i = 0; i < J.size(); i++)
      model.addConstr(s[0][J[i]] == 1);
    
    // Round function
    for (int r = 0; r < evalNumRounds; r++) {
      vector<GRBVar> tmp = s[r];
      triviumCoreThree(model, tmp, 65, 170, 90, 91, 92);
      triviumCoreThree(model, tmp, 161, 263, 174, 175, 176);
      triviumCoreThree(model, tmp, 242, 68, 285, 286, 287);

      for (int i = 0; i < 288; i++)
        s[r + 1][(i + 1) % 288] = tmp[i];
      
      if ((r + 1) == 210) {
				GRBVar p = model.addVar(0, 1, 0, GRB_BINARY);
				GRBVar newvar = model.addVar(0, 0, 0, GRB_BINARY);
				GRBVar tmp[2] = { p,newvar };
				model.addGenConstrOr(s[210][93], tmp, 2);
				s[210][93] = newvar;

				GRBVar q = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(q == o + p);
				model.addConstr(q == 1);
      }
    }

    // Output constraint
    if (target == 0) {
      GRBLinExpr ks = 0;
      for (int i = 0; i < 288; i++) {
        if ((i == 65) || (i == 92) || (i == 161) || (i == 176) || (i == 242) || (i == 287))
          ks += s[evalNumRounds][i];
        else 
          model.addConstr(s[evalNumRounds][i] == 0);
      }
      model.addConstr(ks == 1);
    }
    else if (target == 1) {
      for (int i = 0; i < 288; i++) {
        if ((i == 65)) 
          model.addConstr(s[evalNumRounds][i] == 1);
        else 
          model.addConstr(s[evalNumRounds][i] == 0);
      }
    }
    else if (target == 2) {
      for (int i = 0; i < 288; i++) {
        if ((i == 92)) 
          model.addConstr(s[evalNumRounds][i] == 1);
        else 
          model.addConstr(s[evalNumRounds][i] == 0);
      }
    }
    else if (target == 3) {
      for (int i = 0; i < 288; i++) {
        if ((i == 161)) 
          model.addConstr(s[evalNumRounds][i] == 1);
        else 
          model.addConstr(s[evalNumRounds][i] == 0);
      }
    }
    else if (target == 4) {
      for (int i = 0; i < 288; i++) {
        if ((i == 176)) 
          model.addConstr(s[evalNumRounds][i] == 1);
        else 
          model.addConstr(s[evalNumRounds][i] == 0);
      }
    }
    else if (target == 5) {
      for (int i = 0; i < 288; i++) {
        if ((i == 242)) 
          model.addConstr(s[evalNumRounds][i] == 1);
        else 
          model.addConstr(s[evalNumRounds][i] == 0);
      }
    }
    else if (target == 6) {
      for (int i = 0; i < 288; i++) {
        if ((i == 287))
          model.addConstr(s[evalNumRounds][i] == 1);
        else 
          model.addConstr(s[evalNumRounds][i] == 0);
      }
    }

    // dummy objective function
    GRBLinExpr sumMiddle = 0;
    for (int i = 0; i < 288; i++)
      sumMiddle += s[evalNumRounds / 2][i];
    model.setObjective(sumMiddle, GRB_MAXIMIZE);

    // Solve
    model.update();
    model.optimize();

    //
		int solCount = model.get(GRB_IntAttr_SolCount);
    double dulation = model.get(GRB_DoubleAttr_Runtime);
    cout << "Running time = " << dulation << " sec." << endl;
    cout << "There are " << solCount << " solutions." << endl;

		//
		ofstream outputSolution;
		outputSolution.open("solution.txt", ios::app);
		for (int i = 0; i < solCount; i++) {
			outputSolution << "++++++++++++++++++++++++++++" << endl;
			model.set(GRB_IntParam_SolutionNumber, i);
			for (int r = 0; r < evalNumRounds; r++) {
				for (int j = 0; j < 288; j++) {
					if (round(s[r][j].get(GRB_DoubleAttr_Xn)) == 1) outputSolution << "1";
					else outputSolution << "0";

					if ((j == 92) || (j == 176))
						outputSolution << " ";
				}
				outputSolution << endl;
			}
		}

    // return
    if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
      outputfile << endl;
      return -1;
    }
    else if ((model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)) {
      int upperBound = round(model.get(GRB_DoubleAttr_ObjVal));
      outputfile << endl;
      return upperBound;
    }
    else {
      cout << model.get(GRB_IntAttr_Status) << endl;
      outputfile << endl;
      return -2;
    }
  }
  catch (GRBException e) {
    cerr << "Error code = " << e.getErrorCode() << endl;
    cerr << e.getMessage() << endl;
  }
  catch (...) {
    cerr << "Exception during optimization" << endl;
  }

  return -1;
}

int main(int argc, char const* argv[]){

  int evalNumRounds = 0;
  int threadNumber = 2;
  ofstream outputfile;
  outputfile.open("log.txt");
  ofstream outputSolution;
  outputSolution.open("solution.txt");

  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-r")) evalNumRounds = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "-t")) threadNumber = atoi(argv[i + 1]);
  }

  if (evalNumRounds == 0){
    cerr << "Please set option about number of rounds as '-r [number of rounds]]'" << endl;
    return 0;
  }

  cerr << threadNumber << " cores are used. " << endl;
  cerr << "if you want to change the core number, please set option as '-t [core number]'" << endl;

  // Fu's cube (v74, v60, v75, v30, and v48 are fixed to 0, and other bits are active)
  cout << "index of cube I" << endl;
  vector<int> I;
	for (int i = 0; i < 80; i++) {
		if ((i == 74) || (i == 60) || (i == 75) || (i == 30) || (i == 48)) {
      printf("   , ");
		}
		else {
			I.push_back(i);
      printf("v%02d, ", (i + 1));
		}
    if (i % 10 == 9)
      cout << endl;
	}
  cout << endl;

  //
  cout << "target monomial" << endl;
  vector<int> J = { 40,41,42,53,54,55,56,57,58,61,62,63,65,66,67,68,69,70,71,72,73,74,75,76,78,79 };
  for (int i = 0; i < J.size(); i++) {
    printf("x%02d * ", (J[i] + 1));
  }
  for (int i = 0; i < I.size() - 1; i++) {
    printf("v%02d * ", (I[i] + 1));
  }
  printf("v%02d\n", (I[I.size() - 1] + 1));
  cout << endl;

	//
  cout << "Target s288" << endl;
  veryfi855(I, J, evalNumRounds, 6, threadNumber);
  cout << endl;

  cout << "Target s177" << endl;
  veryfi855(I, J, evalNumRounds, 4, threadNumber);
  cout << endl;

  cout << "Target s93" << endl;
  veryfi855(I, J, evalNumRounds, 2, threadNumber);
  cout << endl;
	
  cout << "Target s243" << endl;
  veryfi855(I, J, evalNumRounds, 5, threadNumber);
  cout << endl;
  
  cout << "Target s162" << endl;
  veryfi855(I, J, evalNumRounds, 3, threadNumber);
  cout << endl;

  cout << "Target s66" << endl;
  veryfi855(I, J, evalNumRounds, 1, threadNumber);
  cout << endl;
  
  return 0;
}
