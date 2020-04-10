#include"main.h"

//Regard two 288-bit vectors a and b as ordinary integers and compare them to see whether a<b
struct cmpBitset288 {
	bool operator()(const bitset<288>& a, const bitset<288>& b) const {
		for (int i = 0; i < 288; i++) {
			if (a[i] < b[i])
				return true;
			else if (a[i] > b[i])
				return false;
		}
		return false;
	}
};
/*
The 3-subset division property requires to evaluate 2 parameters: the monomial u and the number of division trails corresponding to u, denoted as J[u].
If J[u] is EVEN, monomial u is cancelled and cannot appear in the superpoly;
If J[u] is ODD, monomial u is guaranteed to appear in the superpoly.
Both u and J[u] are to be determined by constructing a MILP model and solve with the Gurobi solvers. 
For R-round Trivium, the MILP model for evaluating J[u] can be quite complicated to be solved. 
In order to compute J[u], we cut the full R-round Trivium into 2 stages:
0->R/2->R
u->k'->1
If the number of division trail u->k is J[u], the number of division trail from u->k' is J[u->k'] and k'->1 is J'[k'->1], we have
J[u]=\sum_{k'} (J[u->k']*J[k'->k])
Therefore, the 1st stage is to compute J[u->k'] and the 2nd stage is to find J[k'->k]
The twoTwoStage structure is to handle such a process. 
@Para:
userTwoStage: true if we use such a two stage strategy
divRound: usually R/2. Of course, other number is also OK. 
hint: the respresentation of k'
*/
struct twoStage {
	bool useTwoStage;
	int divRound;
	vector<bitset<288>> hint;
};

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

	GRBVar y1 = model.addVar(0, 1, 0, GRB_BINARY);
	GRBVar y2 = model.addVar(0, 1, 0, GRB_BINARY);
	GRBVar y3 = model.addVar(0, 1, 0, GRB_BINARY);
	GRBVar y4 = model.addVar(0, 1, 0, GRB_BINARY);
	GRBVar y5 = model.addVar(0, 1, 0, GRB_BINARY);

	GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
	GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);
	//GRBVar z3 = model.addVar(0, 1, 0, GRB_BINARY);
	//GRBVar z4 = model.addVar(0, 1, 0, GRB_BINARY);

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
/*
cube: cube indices
flag: the situation of the 288 bits at current round: secret 3, active 2, const 1, const 0
evalNumRounds: current round
countingBox: map each monomial u(represented as a 288-bit vector) to an integer J[u] representing the number of available division trails 
dulation: time for solving the model
threadNumber: the number of threads used for solving the model
target: 0: evaluate directly the exact output z=\sum ss[66,93,162,177,243,288]; 1-6 corresponding to s[66,93,162,177,243,288] resepectively to save the solving time.  
opt: tell the solver to construct and solve the model corresponding to the 1st or 2nd stage
*/
int triviumThreeEnumuration(vector<int> cube, vector<int> flag, int evalNumRounds, map<bitset<288>, int, cmpBitset288>& countingBox, double& dulation, int threadNumber, int target = 0, struct twoStage opt = { false, 0, });
class threeEnumuration : public GRBCallback
{
public:
	vector<int> cube;
	vector<int> flag;
	vector<vector<GRBVar>> s;
	int target;
	map<bitset<288>, int, cmpBitset288>* countingBox;
	int threadNumber;
	ofstream* outputfile;
	threeEnumuration(vector<int> xcube, vector<int> xflag, vector<vector<GRBVar>> xs, int xtarget, map<bitset<288>, int, cmpBitset288>* xcountingBox, int xthreadNumber, ofstream* xoutputfile) {
		cube = xcube;
		flag = xflag;
		s = xs;
		target = xtarget;
		countingBox = xcountingBox;
		threadNumber = xthreadNumber;
		outputfile = xoutputfile;
	}
protected:
	void callback() {
		try {
			if (where == GRB_CB_MIPSOL) {

				int evalNumRounds = s.size() - 1;
				int divRound = evalNumRounds / 2;

				*outputfile << "found \t divide in " << divRound << "\t" << getDoubleInfo(GRB_CB_RUNTIME) << "sec" << endl;

				// store found solution into trail
				vector<bitset<288>> trail(evalNumRounds + 1);
				for (int r = 0; r <= evalNumRounds; r++) {
					for (int i = 0; i < 288; i++) {
						if (round(getSolution(s[r][i])) == 1) trail[r][i] = 1;
						else trail[r][i] = 0;
					}
				}

				// 2nd stage
				double dulation = 0;
				int solCnt = triviumThreeEnumuration(cube, flag, evalNumRounds, *countingBox, dulation, threadNumber, target, { true, divRound, trail });

				//
				int solTotal = 0;
				auto it = (*countingBox).begin();
				while (it != (*countingBox).end()) {
					solTotal += (*it).second;
					it++;
				}
				(*outputfile) << "\t" << solCnt << "( total : " << solTotal << ")" << endl;
				(*outputfile) << "\t" << (*countingBox).size() << " monomials are involved" << endl;

				// remove
				GRBLinExpr addCon = 0;
				for (int i = 0; i < 288; i++) {
					if (round(getSolution(s[divRound][i])) == 1) {
						addCon += (1 - s[divRound][i]);
					}
					else {
						addCon += s[divRound][i];
					}
				}
				addLazy(addCon >= 1);

			}
			else if (where == GRB_CB_MESSAGE) {
				// Message callback
				string msg = getStringInfo(GRB_CB_MSG_STRING);
				*outputfile << msg << flush;
			}
		}
		catch (GRBException e) {
			cerr << "Error number: " << e.getErrorCode() << endl;
			cerr << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}
	}
};
int triviumThreeEnumuration(vector<int> cube, vector<int> flag, int evalNumRounds, map<bitset<288>, int, cmpBitset288>& countingBox, double& dulation, int threadNumber, int target, struct twoStage opt) {

	//
	ofstream outputfile;
	if ((opt.useTwoStage == true) && (opt.hint.size() == 0)) {
		outputfile.open("log_trivium.txt", ios::app);
	}
	else if ((opt.useTwoStage == true) && (opt.hint.size() > 0)) {
		outputfile.open("log_trivium2.txt", ios::app);
	}
	else {
		outputfile.open("log_trivium.txt", ios::app);
	}


	//
	if (opt.useTwoStage == true) {
		if (opt.hint.size() == 0) {
			outputfile << endl;
			outputfile << "++++++++++++++++++++++++++++" << endl;
			outputfile << "1st stage" << endl;
		}
		else {
			outputfile << "---" << endl;
			outputfile << "2nd stage" << endl;
		}
	}

	//gurobi
	try {
		// Create the environment
		GRBEnv env = GRBEnv();

		// close standard output
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_Threads, threadNumber);
		env.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_BESTBOUND);

		if ((opt.useTwoStage == true) && (opt.hint.size() == 0)) {
			env.set(GRB_IntParam_LazyConstraints, 1);
		}else if ((opt.useTwoStage == true) && (opt.hint.size() > 0)) {
			env.set(GRB_StringParam_LogFile, "log_trivium2.txt");
			env.set(GRB_IntParam_PoolSearchMode, 2);
			env.set(GRB_IntParam_PoolSolutions, 2000000000);
			env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);
		}
		else {
			env.set(GRB_StringParam_LogFile, "log_trivium.txt");
			env.set(GRB_IntParam_PoolSearchMode, 2);
			env.set(GRB_IntParam_PoolSolutions, 2000000000);
			env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);
		}

		
		// Create the model
		GRBModel model = GRBModel(env);

		// Create variables
		vector<vector<GRBVar>> s(evalNumRounds + 1, vector<GRBVar>(288));
		for (int i = 0; i < 288; i++) {
			s[0][i] = model.addVar(0, 1, 0, GRB_BINARY);
		}

		// IV constraint
		for (int i = 0; i < 80; i++) {
			if (cube[i] == 1)
				model.addConstr(s[0][93 + i] == 1);
		}

		// Const 0, constraint
		for (int i = 0; i < 288; i++) {
			if (flag[i] == 0)
				model.addConstr(s[0][i] == 0);
		}

		// Round function
		for (int r = 0; r < evalNumRounds; r++) {
			vector<GRBVar> tmp = s[r];
			triviumCoreThree(model, tmp, 65, 170, 90, 91, 92);
			triviumCoreThree(model, tmp, 161, 263, 174, 175, 176);
			triviumCoreThree(model, tmp, 242, 68, 285, 286, 287);
			
			for (int i = 0; i < 288; i++) {
				s[r + 1][(i + 1) % 288] = tmp[i];
			}
		}

		// Output constraint
		if (target == 0) {
			GRBLinExpr ks = 0;
			for (int i = 0; i < 288; i++) {
				if ((i == 65) || (i == 92) || (i == 161) || (i == 176) || (i == 242) || (i == 287)) {
					ks += s[evalNumRounds][i];
				}
				else {
					model.addConstr(s[evalNumRounds][i] == 0);
				}
			}
			model.addConstr(ks == 1);
		}
		else if (target == 1) {
			for (int i = 0; i < 288; i++) {
				if ((i == 65)) {
					model.addConstr(s[evalNumRounds][i] == 1);
				}
				else {
					model.addConstr(s[evalNumRounds][i] == 0);
				}
			}
		}
		else if (target == 2) {
			for (int i = 0; i < 288; i++) {
				if ((i == 92)) {
					model.addConstr(s[evalNumRounds][i] == 1);
				}
				else {
					model.addConstr(s[evalNumRounds][i] == 0);
				}
			}
		}
		else if (target == 3) {
			for (int i = 0; i < 288; i++) {
				if ((i == 161)) {
					model.addConstr(s[evalNumRounds][i] == 1);
				}
				else {
					model.addConstr(s[evalNumRounds][i] == 0);
				}
			}
		}
		else if (target == 4) {
			for (int i = 0; i < 288; i++) {
				if ((i == 176)) {
					model.addConstr(s[evalNumRounds][i] == 1);
				}
				else {
					model.addConstr(s[evalNumRounds][i] == 0);
				}
			}
		}
		else if (target == 5) {
			for (int i = 0; i < 288; i++) {
				if ((i == 242)) {
					model.addConstr(s[evalNumRounds][i] == 1);
				}
				else {
					model.addConstr(s[evalNumRounds][i] == 0);
				}
			}

		}
		else if (target == 6) {
			for (int i = 0; i < 288; i++) {
				if ((i == 287)) {
					model.addConstr(s[evalNumRounds][i] == 1);
				}
				else {
					model.addConstr(s[evalNumRounds][i] == 0);
				}
			}
		}

		//
		GRBLinExpr sumKey = 0;
		for (int i = 0; i < 80; i++) {
			sumKey += s[0][i];
		}
		model.setObjective(sumKey, GRB_MAXIMIZE);
	
		//
		if (opt.useTwoStage == true) {
			if (opt.hint.size() > 0) {

				// fix
				for (int i = 0; i < 288; i++) {
					if (opt.hint[opt.divRound][i] == 1)
						model.addConstr(s[opt.divRound][i] == 1);
					else
						model.addConstr(s[opt.divRound][i] == 0);
				}

				// hint
				for (int r = 0; r < evalNumRounds; r++) {
					for (int i = 0; i < 288; i++) {
						if (opt.hint[r][i] == 1)
							s[r][i].set(GRB_DoubleAttr_Start, 1);
						else
							s[r][i].set(GRB_DoubleAttr_Start, 0);
					}
				}


			}
		}

		// Solve
		model.update();
		if ((opt.useTwoStage == true) && (opt.hint.size() == 0) ) {
			threeEnumuration cb = threeEnumuration(cube, flag, s, target, &countingBox, threadNumber, &outputfile);
			model.setCallback(&cb);
			model.optimize();
		}
		else {
			model.optimize();
		}


		//
		int solCount = model.get(GRB_IntAttr_SolCount);
		dulation = model.get(GRB_DoubleAttr_Runtime);


		//
		if (opt.useTwoStage == true) {
			if (opt.hint.size() > 0) {
				// check solution limit
				if (solCount >= 2000000000) {
					cerr << "Number of solutions is too large" << endl;
					exit(0);
				}

				// store the information about solutions
				for (int i = 0; i < solCount; i++) {
					model.set(GRB_IntParam_SolutionNumber, i);
					bitset<288> tmp;
					for (int j = 0; j < 288; j++) {
						if (round(s[0][j].get(GRB_DoubleAttr_Xn)) == 1) tmp[j] = 1;
						else tmp[j] = 0;
					}
					countingBox[tmp]++;
				}

				return solCount;
			}
		}
		else {
			// check solution limit
			if (solCount >= 2000000000) {
				cerr << "Number of solutions is too large" << endl;
				exit(0);
			}

			// store the information about solutions
			for (int i = 0; i < solCount; i++) {
				model.set(GRB_IntParam_SolutionNumber, i);
				bitset<288> tmp;
				for (int j = 0; j < 288; j++) {
					if (round(s[0][j].get(GRB_DoubleAttr_Xn)) == 1) tmp[j] = 1;
					else tmp[j] = 0;
				}
				countingBox[tmp]++;
			}
		}



		// display result
		auto it = countingBox.begin();
		while (it != countingBox.end()) {

			cout << ((*it).second % 2) << " | " << (*it).second << "\t";
			bitset<288> tmp = (*it).first;
			for (int i = 0; i < 80; i++) {
				if ((tmp[i] == 1)) {
					cout << "k" << (i + 1) << " ";
				}
			}
			for (int i = 0; i < 80; i++) {
				if ((cube[i] == 0) && (tmp[93 + i] == 1)) {
					cout << "v" << (i + 1) << " ";
				}
			}
			for (int i = 285; i < 288; i++) {
				if ((tmp[i] == 1)) {
					cout << "c" << (i + 1) << " ";
				}
			}
			cout << endl;

			it++;
		}
		cout << dulation << "sec" << endl;



		//result
		if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
			return -1;
		}
		else if ((model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)) {
			int upperBound = round(model.get(GRB_DoubleAttr_ObjVal));
			return upperBound;
		}
		else {
			cout << model.get(GRB_IntAttr_Status) << endl;
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
int trivium(int evalNumRounds, int threadNumber) {

  //
  ofstream outputfile, outputfile2;
  outputfile.open("log_trivium.txt");
  outputfile2.open("log_trivium2.txt");

  // active bits
  vector<int> cube(80, 0);
  for (int i = 0; i < 80; i++) {

		if (evalNumRounds == 840) {
			if ((i == 33) || (i == 46)) {

			}
			else {
				cube[i] = 1;
			}
		}else if (evalNumRounds == 841) {
			if ((i == 8) || (i == 78)) {

			}
			else {
				cube[i] = 1;
			}
		}
		else {
			cube[i] = 1;
		}
		
  }
  cout << "the index of cube" << endl;
  for (int i = 0; i < 80; i++) {
    if(cube[i] == 1)
      cout << "iv" << (i + 1) << ", ";
  }
  cout << endl;

	// constant (secret 3, active 2, const 1, const 0)
  vector<int> flag(288, 0);
  for (int i = 0; i < 80; i++) {
    flag[i] = 3;
  }
  for (int i = 0; i < 80; i++) {
    if (cube[i] == 1) {
      flag[93 + i] = 2;
    }
  }
  flag[285] = 1;
  flag[286] = 1;
  flag[287] = 1;
	



  // 
  map<bitset<288>, int, cmpBitset288> countingBox;
  double dulation = 0;

  //
  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Target s288" << endl;
	triviumThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 6, { true,0, });
  if (countingBox.size() == 0) {
    cout << "zero sum\t" << dulation << "sec" << endl;
  }

  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Target s177" << endl;
  triviumThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 4, { true,0, });
  if (countingBox.size() == 0) {
    cout << "zero sum\t" << dulation << "sec" << endl;
  }

  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Target s93" << endl;
  triviumThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 2, { true,0, });
  if (countingBox.size() == 0) {
    cout << "zero sum\t" << dulation << "sec" << endl;
  }

  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Target s243" << endl;
  triviumThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 5, { true,0, });
  if (countingBox.size() == 0) {
    cout << "zero sum\t" << dulation << "sec" << endl;
  }

  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Target s162" << endl;
  triviumThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 3, { true,0, });
  if (countingBox.size() == 0) {
    cout << "zero sum\t" << dulation << "sec" << endl;
  }

  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Target s66" << endl;
  triviumThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 1, { true,0, });
  if (countingBox.size() == 0) {
    cout << "zero sum\t" << dulation << "sec" << endl;
  }



	cout << "*****************************" << endl;
	cout << "Final solution" << endl;
	cout << countingBox.size() << " solutions are found" << endl;

	map<bitset<288>, int, cmpBitset288> countingBox2;
	auto it = countingBox.begin();
	while (it != countingBox.end()) {
		bitset<288> tmp = (*it).first;
		for (int i = 0; i < 80; i++) {
			if(cube[i] == 1)
				tmp[93 + i] = 0;
		}
		tmp[285] = 0;
		tmp[286] = 0;
		tmp[287] = 0;
		countingBox2[tmp] += (*it).second;
		it++;
	}

	auto it2 = countingBox2.begin();
	while (it2 != countingBox2.end()) {
		if (((*it2).second % 2) == 1) {
			cout << ((*it2).second % 2) << " | " << (*it2).second << "\t";
			bitset<288> tmp = (*it2).first;
			for (int i = 0; i < 80; i++) {
				if ((tmp[i] == 1)) {
					cout << "k" << (i + 1) << " ";
				}
			}
			cout << endl;
		}
		it2++;
	}


	it2 = countingBox2.begin();
	while (it2 != countingBox2.end()) {
		if (((*it2).second % 2) == 0) {
			cout << ((*it2).second % 2) << " | " << (*it2).second << "\t";
			bitset<288> tmp = (*it2).first;
			for (int i = 0; i < 80; i++) {
				if ((tmp[i] == 1)) {
					cout << "k" << (i + 1) << " ";
				}
			}
			cout << endl;
		}
		it2++;
	}


}

// for the practical verification
int roundFuncTrivium(bitset<288>& s) {

  bitset<288> o = s;

  int x1 = o[92] ^ o[65];
  int x2 = o[176] ^ o[161];
  int x3 = o[287] ^ o[242];
  int z = x1 ^ x2 ^ x3;

  o[92] = x1 ^ (o[91] & o[90]) ^ o[170];
  o[176] = x2 ^ (o[174] & o[175]) ^ o[263];
  o[287] = x3 ^ (o[286] & o[285]) ^ o[68];


  s = (o << 1) ^ (o >> 287);

  return z;
}
int encryptionSum(int evalNumRounds, vector<int> cube, vector<int> iv, vector<int> key) {

	bitset<288> s;
	for (int i = 0; i < 80; i++) {
		s[i] = key[i];
	}
	for (int i = 0; i < 80; i++) {
		s[93 + i] = iv[i];
	}
	s[285] = 1;
	s[286] = 1;
	s[287] = 1;

	int DATA_SIZE = 0;
	vector<int> map;
	for (int i = 0; i < cube.size(); i++) {
		DATA_SIZE += cube[i];
		if (cube[i] == 1) {
			map.push_back(i);
		}
	}

	int sum = 0;
	for (int in = 0; in < (1 << DATA_SIZE); in++) {

		bitset<288> tmp = s;
		for (int i = 0; i < DATA_SIZE; i++) {
			tmp[93 + map[i]] = ((in >> i) & 1);
		}

		int z = 0;
		for (int r = 0; r <= evalNumRounds; r++) {
			z = roundFuncTrivium(tmp);
		}
		sum ^= z;
	}

	return sum;
}
int theoreticalSum(map<bitset<288>, int, cmpBitset288> countingBox, vector<int> cube, vector<int> iv, vector<int> key) {

	int sum = 0;
	auto it = countingBox.begin();
	while (it != countingBox.end()) {

		if (((*it).second % 2) == 1) {

			int var = 1;
			for (int i = 0; i < 288; i++) {
				if ((*it).first[i] == 1) {
					if (i < 80) {
						var *= key[i];
					}
					else if( (93 <= i) && (i < 93 + 80) ) {
						if (cube[i - 93] == 0) {
							var *= iv[i - 93];
						}
					}
				}
			}
			sum ^= var;

		}


		it++;
	}
	return sum;

}
void practicalTestTrivium(void) {


	ofstream outputfile, outputfile2;
	outputfile.open("log_trivium.txt");
	outputfile2.open("log_trivium2.txt");

  srand(time(NULL));

  // create cube index at random
  int numActBits = 1 + rand() % 4;
  cout << numActBits << " active bits" << endl;

  vector<int> cube(80,0);
  for (int i = 0; i < numActBits; i++) {
    int index;
    do {
      index = rand() % 80;
    } while (cube[index] == 1);
    cube[index] = 1;
  }

  // C_0 = null
	// constant (secret 3, active 2, const/undermine 1, const 0)
	vector<int> flag(288, 0);
  for (int i = 0; i < 80; i++) {
    flag[i] = 3;
    if (cube[i] == 1) {
      flag[93 + i] = 2;
		}
		else {
			flag[93 + i] = 1;
		}
  }
  flag[285] = 1;
  flag[286] = 1;
  flag[287] = 1;
  
	//
	cout << "the index of cube" << endl;
	for (int i = 0; i < 80; i++) {
		if (cube[i] == 1)
			cout << "iv" << (i + 1) << ", ";
	}
	cout << endl;

  //
  for (int r = 300; r < 600; r++) {
    cout << "##############################" << endl;
    cout << r << " rounds" << endl;

    double dulation;
    map<bitset<288>, int, cmpBitset288> countingBox;
		triviumThreeEnumuration(cube, flag, r, countingBox, dulation, 2);
		
    //triviumThreeEnumuration(cube, flag, r, countingBox, dulation, 2, 1);
		//triviumThreeEnumuration(cube, flag, r, countingBox, dulation, 2, 2);
		//triviumThreeEnumuration(cube, flag, r, countingBox, dulation, 2, 3);
		//triviumThreeEnumuration(cube, flag, r, countingBox, dulation, 2, 4);
		//triviumThreeEnumuration(cube, flag, r, countingBox, dulation, 2, 5);
		//triviumThreeEnumuration(cube, flag, r, countingBox, dulation, 2, 6);
		
    if (countingBox.size() == 0) {
      cout << "zero sum" << endl;
    }
    else {

			cout << "         key        \t";
			cout << "          iv          ";
			cout << "expe.   ";
			cout << "theo.   ";
			cout << endl;

      for (int trial = 0; trial < 100; trial++) {

        vector<int> key(80);
        for (int i = 0; i < 80; i++)
          key[i] = rand() % 2;


				vector<int> iv(80);
				for (int i = 0; i < 80; i++)
					iv[i] = rand() % 2;
        


        int sum1 = encryptionSum(r, cube, iv, key);
        int sum2 = theoreticalSum(countingBox, cube, iv, key);

				for (int i = 9; i >= 0; i--) {
					int hexvar = 0;
					for (int j = 7; j >= 0; j--) {
						hexvar ^= (key[8 * i + j] << j);
					}
					printf("%02x", hexvar);
				}
				cout << "\t";

				for (int i = 9; i >= 0; i--) {
					int hexvar = 0;
					for (int j = 7; j >= 0; j--) {
						hexvar ^= (iv[8 * i + j] << j);
					}
					printf("%02x", hexvar);
				}
				cout << "\t";

				cout << sum1 << "\t" << sum2 << "\t";

				if (sum1 == sum2) {
					cout << "OK" << endl;
				}
				else {
					cout << endl;
					cout << "error" << endl;
					cerr << "error" << endl;
				}

      }

    }


    cout << endl << endl;
  }



}