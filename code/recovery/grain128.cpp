#include"main.h"

//Regard two 256-bit vectors a and b as ordinary integers and compare them to see whether a<b
struct cmpBitset256 {
  bool operator()(const bitset<256>& a, const bitset<256>& b) const {
    for (int i = 0; i < 256; i++) {
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
struct twoStageGrain {
	bool useTwoStage;
	int divRound;
	vector<bitset<256>> hint;
};


/***************************************
 * MILP model for each components
 ***************************************/
/*
The tap operation of NLFSRs. It corresponds to the COPY operation of division property:
A variable x is first copied to x->(y,z). x is replaced with z and y will be involved in other operations (such as AND, XOR etc.)
@Para
model: the MILP model
x: the 3-subset division property of the variable to be tapped.
*/
GRBVar tap(GRBModel& model, GRBVar& x) {

  GRBVar y = model.addVar(0, 1, 0, GRB_BINARY);
  GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
  GRBVar tmp[2] = { y,z };
  model.addGenConstrOr(x, tmp, 2);
  x = z;
  return y;

}



/***************************************
 * MILP model for the H function
 ***************************************/
/*
The h function of the Grain-128a stream cihper--the non-linear part of the output function z=h+LinearPart
h= b12s8 + s13s20 + b95s42 + s60s79 + b12b95s94
An additional parameter "target" is introduced enabling to evaluate the division property of each term b12s8, s12s20 ... separately.
For example, if we set target=1, we regard the the output bit as z'=b12s8 so the corresponding superpoly is only related to b12s8. 
target=2,...5 corresponds to the rest terms s13s20,..., b12b95s94
target=6 corresponds to the linear part of the output function
target=-1 is used during the initialization phase.
@Para
model: the MILP model
b: the three-subset division property for the NFSR
s: the three-subset division property for the LFSR
target: the evaluation target: -1->initialization round, 6->ignore the non-linear part and only regard the linear parts as the output, 1-5->corresponds to the evaluation of b12s8,...,b12b94s94 separately
*/
GRBVar funcH(GRBModel& model, vector<GRBVar>& b, vector<GRBVar>& s, int target = -1) {

  GRBVar b12x = tap(model, b[12]);
  GRBVar s8 = tap(model, s[8]);
  model.addConstr(b12x == s8);

  GRBVar s13 = tap(model, s[13]);
  GRBVar s20 = tap(model, s[20]);
  model.addConstr(s13 == s20);

  GRBVar b95x = tap(model, b[95]);
  GRBVar s42 = tap(model, s[42]);
  model.addConstr(b95x == s42);

  GRBVar s60 = tap(model, s[60]);
  GRBVar s79 = tap(model, s[79]);
  model.addConstr(s60 == s79);

  GRBVar b12y = tap(model, b[12]);
  GRBVar b95y = tap(model, b[95]);
  GRBVar s94 = tap(model, s[94]);
  model.addConstr(b12y == b95y);
  model.addConstr(b12y == s94);

  GRBVar y = model.addVar(0, 1, 0, GRB_BINARY);
  model.addConstr(y == b12x + s13 + b95x + s60 + b12y);

	if (target == 1)
		model.addConstr(b12y == 1);
	else if (target == 2)
		model.addConstr(s60 == 1);
	else if (target == 3)
		model.addConstr(b95x == 1);
	else if (target == 4)
		model.addConstr(s13 == 1);
	else if (target == 5)
		model.addConstr(b12x == 1);
	else if (target == 6)
		model.addConstr(y == 0);

  return y;

}

/***************************************
 * MILP model for the output function 
 ***************************************/
/*
The linear part of the output function
LinearPart= s93 + b2 + b15 + b36 + b45 + b64 + b73 + b89
The parameter "target" corresponds to the "target" in "funcH" so 
-1 corresponds to the initialization phase, 
1-5 corresponds to the 5 non-linear terms in h function so there is an additional y=0 constraint in the fundO
6 corresponds the situation that regards the linear part as the output so there is and additional y=1 constraint in funcO
@Para
model: the MILP model
b: the three-subset division property for the NFSR
s: the three-subset division property for the LFSR
target: the evaluation target: -1->initialization round, 6->ignore the non-linear part and only regard the linear parts as the output, 1-5->corresponds to the evaluation of b12s8,...,b12b94s94 separately
*/
GRBVar funcO(GRBModel& model, vector<GRBVar>& b, vector<GRBVar>& s, int target = -1) {
  GRBVar s93 = tap(model, s[93]);
  GRBVar b2 = tap(model, b[2]);
  GRBVar b15 = tap(model, b[15]);
  GRBVar b36 = tap(model, b[36]);
  GRBVar b45 = tap(model, b[45]);
  GRBVar b64 = tap(model, b[64]);
  GRBVar b73 = tap(model, b[73]);
  GRBVar b89 = tap(model, b[89]);

  GRBVar y = model.addVar(0, 1, 0, GRB_BINARY);
  model.addConstr(y == s93 + b2 + b15 + b36 + b45 + b64 + b73 + b89);  

	if (target == 6)
		model.addConstr(y == 1);
	else if (target >= 1)
		model.addConstr(y == 0);

  return y;
}

/***************************************
 * MILP model for the F function
 ***************************************/
/*
The updating function of the LFSR
f=s0 + s7 + s38 + s70 + s81 + s96
@Para
model: the MILP model
s: the three-subset division property for the LFSR
*/
GRBVar funcF(GRBModel& model, vector<GRBVar>& s) {
  GRBVar s0 = tap(model, s[0]);
  GRBVar s7 = tap(model, s[7]);
  GRBVar s38 = tap(model, s[38]);
  GRBVar s70 = tap(model, s[70]);
  GRBVar s81 = tap(model, s[81]);
  GRBVar s96 = tap(model, s[96]);

  GRBVar f = model.addVar(0, 1, 0, GRB_BINARY);
  model.addConstr(f == s0 + s7 + s38 + s70 + s81 + s96);

  return f;
}

/***************************************
 * MILP model for the G function
 ***************************************/
/*
The updating function of the LFSR
g=b0 + b26 + b56 + b91 + b96 + b3b67 + b11b13 + b17b18 + b27b59 + b40b48 + b61b65 + b68b84 + b88b92b93b95 + b22b24b25 + b70b78b82 
@Para
model: the MILP model
b: the three-subset division property for the NFSR
*/
GRBVar funcG(GRBModel& model, vector<GRBVar>& b) {
  // nonlinear
  GRBVar b26 = tap(model, b[26]);
  GRBVar b56 = tap(model, b[56]);
  GRBVar b91 = tap(model, b[91]);
  GRBVar b96 = tap(model, b[96]);

  GRBVar b3 = tap(model, b[3]);
  GRBVar b67 =  tap(model, b[67]);
  model.addConstr(b3 == b67);

  GRBVar b11 = tap(model, b[11]);
  GRBVar b13 = tap(model, b[13]);
  model.addConstr(b11 == b13);

  GRBVar b17 = tap(model, b[17]);
  GRBVar b18 = tap(model, b[18]);
  model.addConstr(b17 == b18);

  GRBVar b27=  tap(model, b[27]);
  GRBVar b59= tap(model, b[59]);
  model.addConstr(b27 == b59);

  GRBVar b40=  tap(model, b[40]);
  GRBVar b48=  tap(model, b[48]);
  model.addConstr(b40 == b48);

  GRBVar b61=  tap(model, b[61]);
  GRBVar b65=  tap(model, b[65]);
  model.addConstr(b61 == b65);

  GRBVar b68=  tap(model, b[68]);
  GRBVar b84=  tap(model, b[84]);
  model.addConstr(b68 == b84);

  GRBVar b88=  tap(model, b[88]);
  GRBVar b92=  tap(model, b[92]);
  GRBVar b93=  tap(model, b[93]);
  GRBVar b95=  tap(model, b[95]);
  model.addConstr(b88 == b92);
  model.addConstr(b88 == b93);
  model.addConstr(b88 == b95);

  GRBVar b22=  tap(model, b[22]);
  GRBVar b24=  tap(model, b[24]);
  GRBVar b25=  tap(model, b[25]);
  model.addConstr(b22 == b24);
  model.addConstr(b22 == b25);

  GRBVar b70=  tap(model, b[70]);
  GRBVar b78=  tap(model, b[78]);
  GRBVar b82=  tap(model, b[82]);
  model.addConstr(b70 == b78);
  model.addConstr(b70 == b82);

  // nonlinear feed back
  GRBVar g = model.addVar(0, 1, 0, GRB_BINARY);
  model.addConstr(g == b[0] + b26 + b56 + b91 + b96 + b3 + b11 + b17 + b27 + b40 + b61 + b68 + b88 + b22 + b70);

  return g;
}





/*
The main function for the three-subset division property attack on Grain128a
@Para
cube: the cube indices
flag: the flag values of each state bit: secret 3, active 2, const 1, const 0
evalNumRounds: the number of initialization rounds to be attacked
countingBox: store the vector u and the corresponding trail number J[u]
dulation: the solution time
threadNumber: the number of threads used for solving the MILP model
target: same with the "target" in funcH and funcO
opt: the parameters used in the two-stage strategy
*/
int grainThreeEnumuration(vector<int> cube, vector<int> flag, int evalNumRounds, map<bitset<256>, int, cmpBitset256>& countingBox, double& dulation, int threadNumber, int target = -1, struct twoStageGrain opt = { false, 0, });
/*
The class defining the callback strategy for enumerating the trails to acquire J[u]
*/
class threeEnumurationGrain : public GRBCallback
{
public:
	vector<int> cube;
	vector<int> flag;
	vector<vector<GRBVar>> s;
	vector<vector<GRBVar>> b;
	map<bitset<256>, int, cmpBitset256>* countingBox;
	int threadNumber;
	ofstream* outputfile;
	int target;
	threeEnumurationGrain(vector<int> xcube, vector<int> xflag, vector<vector<GRBVar>> xs, vector<vector<GRBVar>> xb, int xtarget, map<bitset<256>, int, cmpBitset256>* xcountingBox, int xthreadNumber, ofstream* xoutputfile) {
		cube = xcube;
		flag = xflag;
		s = xs;
		b = xb;
		countingBox = xcountingBox;
		threadNumber = xthreadNumber;
		outputfile = xoutputfile;
		target = xtarget;
	}
protected:
	void callback() {
		try {
			if (where == GRB_CB_MIPSOL) {

				int evalNumRounds = s.size() - 1;
				int divRound = evalNumRounds / 2;

				(*outputfile) << "\tfound \t divide in " << divRound << "\t" << getDoubleInfo(GRB_CB_RUNTIME) << "sec" << endl;


				// store found solution into trail
				vector<bitset<256>> trail(evalNumRounds + 1);
				for (int r = 0; r <= evalNumRounds; r++) {
					for (int i = 0; i < 128; i++) {
						if (round(getSolution(b[r][i])) == 1) trail[r][i] = 1;
						else trail[r][i] = 0;
					}
					for (int i = 0; i < 128; i++) {
						if (round(getSolution(s[r][i])) == 1) trail[r][128 + i] = 1;
						else trail[r][128 + i] = 0;
					}
				}

				//
				double dulation;
				int solCnt = grainThreeEnumuration(cube, flag, evalNumRounds, *countingBox, dulation, threadNumber, target, { true, divRound, trail });

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
				for (int i = 0; i < 128; i++) {
					if (trail[divRound][i] == 1) {
						addCon += (1 - b[divRound][i]);
					}
					else {
						addCon += b[divRound][i];
					}
				}
				for (int i = 0; i < 128; i++) {
					if (trail[divRound][128 + i] == 1) {
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
				(*outputfile) << msg << flush;
			}
		}
		catch (GRBException e) {
			cerr << "Error number: " << e.getErrorCode() << endl;
			cerr << e.getMessage() << endl;
		}
		catch (...) {
			cerr << "Error during callback" << endl;
		}
	}
};
int grainThreeEnumuration(vector<int> cube, vector<int> flag, int evalNumRounds, map<bitset<256>, int, cmpBitset256>& countingBox, double& dulation, int threadNumber, int target, struct twoStageGrain opt) {


	//
	ofstream outputfile;
	if ((opt.useTwoStage == true) && (opt.hint.size() == 0)) {
		outputfile.open("log_grain128a.txt", ios::app);
	}
	else if ((opt.useTwoStage == true) && (opt.hint.size() > 0)) {
		outputfile.open("log_grain128a2.txt", ios::app);
	}
	else {
		outputfile.open("log_grain128a.txt", ios::app);
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
		//env.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_BESTBOUND);

		if ((opt.useTwoStage == true) && (opt.hint.size() == 0)) {
			env.set(GRB_IntParam_LazyConstraints, 1);
		}
		else if ((opt.useTwoStage == true) && (opt.hint.size() > 0)) {
			env.set(GRB_StringParam_LogFile, "log_grain128a2.txt");
			env.set(GRB_IntParam_PoolSearchMode, 2);
			env.set(GRB_IntParam_PoolSolutions, 2000000000);
			env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);
		}
		else {
			env.set(GRB_StringParam_LogFile, "log_grain128a.txt");
			env.set(GRB_IntParam_PoolSearchMode, 2);
			env.set(GRB_IntParam_PoolSolutions, 2000000000);
			env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);
		}

		// Create the model
		GRBModel model = GRBModel(env);

		// Create variables
		vector<vector<GRBVar>> s(evalNumRounds + 1, vector<GRBVar>(128));
		vector<vector<GRBVar>> b(evalNumRounds + 1, vector<GRBVar>(128));
		for (int i = 0; i < 128; i++) {
			s[0][i] = model.addVar(0, 1, 0, GRB_BINARY);
			b[0][i] = model.addVar(0, 1, 0, GRB_BINARY);
		}

		// IV constraint
		for (int i = 0; i < 96; i++) {
			if (cube[i] == 1)
				model.addConstr(s[0][i] == 1);
			else if (flag[128 + i] == 0)
				model.addConstr(s[0][i] == 0);
		}
		model.addConstr(s[0][127] == 0);



		// Round function
		for (int r = 0; r <= evalNumRounds; r++) {
			vector<GRBVar> tmpb = b[r];
			vector<GRBVar> tmps = s[r];

			if (r < evalNumRounds) {
				GRBVar h = funcH(model, tmpb, tmps);
				GRBVar o = funcO(model, tmpb, tmps);

				GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(z == h + o);

				GRBVar z1 = model.addVar(0, 1, 0, GRB_BINARY);
				GRBVar z2 = model.addVar(0, 1, 0, GRB_BINARY);
				GRBVar tmpVar[2] = { z1,z2 };
				model.addGenConstrOr(z, tmpVar, 2);

				GRBVar f = funcF(model, tmps);
				GRBVar g = funcG(model, tmpb);

				GRBVar news = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(news == z1 + f);

				GRBVar newb = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(newb == z2 + g + tmps[0]);

				for (int i = 0; i < 127; i++) {
					b[r + 1][i] = tmpb[i + 1];
					s[r + 1][i] = tmps[i + 1];
				}
				b[r + 1][127] = newb;
				s[r + 1][127] = news;

				// remove (s[r][0], z[r]) = (1,1) 
				model.addConstr((1 - s[r][0]) + (1 - z) >= 1);
			}
			else {
				GRBVar h = funcH(model, tmpb, tmps, target);
				GRBVar o = funcO(model, tmpb, tmps, target);

				GRBVar z = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(z == h + o);

				model.addConstr(z == 1);

				for (int i = 0; i < 128; i++) {
					model.addConstr(tmpb[i] == 0);
					model.addConstr(tmps[i] == 0);
				}
			}
		}


		//
		GRBLinExpr sumKey = 0;
		for (int i = 0; i < 128; i++) {
			sumKey += b[0][i];
		}
		model.setObjective(sumKey, GRB_MAXIMIZE);

		//
		if (opt.useTwoStage == true) {
			if (opt.hint.size() > 0) {

				// fix
				for (int i = 0; i < 128; i++) {
					if (opt.hint[opt.divRound][i] == 1)
						model.addConstr(b[opt.divRound][i] == 1);
					else
						model.addConstr(b[opt.divRound][i] == 0);

					if (opt.hint[opt.divRound][128 + i] == 1)
						model.addConstr(s[opt.divRound][i] == 1);
					else
						model.addConstr(s[opt.divRound][i] == 0);
				}

				// hint
				for (int r = 0; r < evalNumRounds; r++) {
					for (int i = 0; i < 128; i++) {
						if (opt.hint[r][i] == 1)
							b[r][i].set(GRB_DoubleAttr_Start, 1);
						else
							b[r][i].set(GRB_DoubleAttr_Start, 0);

						if (opt.hint[r][128 + i] == 1)
							s[r][i].set(GRB_DoubleAttr_Start, 1);
						else
							s[r][i].set(GRB_DoubleAttr_Start, 0);
					}
				}
			}
		}

		// Solve
		model.update();
		if ((opt.useTwoStage == true) && (opt.hint.size() == 0)) {
			threeEnumurationGrain cb = threeEnumurationGrain(cube, flag, s, b, target, &countingBox, threadNumber, &outputfile);
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
					bitset<256> tmp;
					for (int j = 0; j < 128; j++) {
						if (round(b[0][j].get(GRB_DoubleAttr_Xn)) == 1) tmp[j] = 1;
						else tmp[j] = 0;

						if (round(s[0][j].get(GRB_DoubleAttr_Xn)) == 1) tmp[128 + j] = 1;
						else tmp[128 + j] = 0;
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
				bitset<256> tmp;
				for (int j = 0; j < 128; j++) {
					if (round(b[0][j].get(GRB_DoubleAttr_Xn)) == 1) tmp[j] = 1;
					else tmp[j] = 0;

					if (round(s[0][j].get(GRB_DoubleAttr_Xn)) == 1) tmp[128 + j] = 1;
					else tmp[128 + j] = 0;
				}
				countingBox[tmp]++;
			}
		}



		// disp
		auto it = countingBox.begin();
		while (it != countingBox.end()) {

			cout << ((*it).second % 2) << " | " << (*it).second << "\t";

			bitset<256> tmp = (*it).first;
			for (int i = 0; i < 128; i++) {
				if ((tmp[i] == 1)) {
					cout << "k" << (i + 1) << " ";
				}
			}
			for (int i = 0; i < 96; i++) {
				if ((tmp[128 + i] == 1) && (cube[i] == 0)) {
					cout << "v" << (i + 1) << " ";
				}
			}
			for (int i = 96; i < 128; i++) {
				if (tmp[128 + i] == 1) {
					cout << "v" << (i + 1) << " ";
				}
			}
			cout << endl;

			it++;
		}
		cout << dulation << "sec" << endl;
		cout << endl;

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



/***************************************
 * Computing Superpolies with 96 active bits
 ***************************************/
/*
evalNumRounds: the number of initialization rounds
threadNumber: the number of threads
*/
int grain128a(int evalNumRounds, int threadNumber) {

	//
	ofstream outputfile, outputfile2, outputfile3;
	outputfile.open("log_grain128a.txt");
	outputfile2.open("log_grain128a2.txt");

	// active bits
	vector<int> cube(96, 0);
	for (int i = 0; i < 96; i++) {
		cube[i] = 1;
	}
	//cube[46] = 0;
	cout << "the index of cube" << endl;
	for (int i = 0; i < 96; i++) {
		if (cube[i] == 1)
			cout << "iv" << (i + 1) << ", ";
	}
	cout << endl;


	// constant (secret 3, active 2, const 1, const 0)
	vector<int> flag(256, 0);
	for (int i = 0; i < 128; i++) {
		flag[i] = 3;
	}
	for (int i = 0; i < 96; i++) {
		if (cube[i] == 1) {
			flag[128 + i] = 2;
		}
	}
	for (int i = 96; i < 127; i++) {
		flag[128 + i] = 1;
	}

	//
	map<bitset<256>, int, cmpBitset256> countingBox;
	double dulation = 0;

	//Seperately evaluate the non-linear terms and the linear part of the output bit
	cout << "++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Target s93 + b2 + b15 + b36 + b45 + b64 + b73 + b89" << endl;
	grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 6, { true, 0, });
	if (countingBox.size() == 0) {
		cout << "zero sum\t" << dulation << "sec" << endl;
	}

	cout << "++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Target b12 * s8" << endl;
	grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 5, { true, 0, });
	if (countingBox.size() == 0) {
		cout << "zero sum\t" << dulation << "sec" << endl;
	}

	cout << "++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Target s13 * s20" << endl;
	grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 4, { true, 0, });
	if (countingBox.size() == 0) {
		cout << "zero sum\t" << dulation << "sec" << endl;
	}

	cout << "++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Target b95 * s42" << endl;
	grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 3, { true, 0, });
	if (countingBox.size() == 0) {
		cout << "zero sum\t" << dulation << "sec" << endl;
	}

	cout << "++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Target s60 * s79" << endl;
	grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 2, { true, 0, });
	if (countingBox.size() == 0) {
		cout << "zero sum\t" << dulation << "sec" << endl;
	}

	cout << "++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Target b12 * b95 * s94" << endl;
	grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 1, { true, 0, });
	if (countingBox.size() == 0) {
		cout << "zero sum\t" << dulation << "sec" << endl;
	}

	

	cout << "*****************************" << endl;
	cout << "Final solution" << endl;
	cout << countingBox.size() << " solutions are found" << endl;

	map<bitset<256>, int, cmpBitset256> countingBox2;
	auto it = countingBox.begin();
	while (it != countingBox.end()) {
		bitset<256> tmp = (*it).first;
		for (int i = 0; i < 96; i++) {
			if (cube[i] == 1) {
				tmp[128 + i] = 0;
			}
		}
		for (int i = 96; i < 128; i++) {
			tmp[128 + i] = 0;
		}
		countingBox2[tmp] += (*it).second;
		it++;
	}

	// 
	cout << "odd list" << endl;
	auto it2 = countingBox2.begin();
	while (it2 != countingBox2.end()) {
		if (((*it2).second % 2) == 1) {
			cout << ((*it2).second % 2) << " | " << (*it2).second << "\t";
			bitset<256> tmp = (*it2).first;
			for (int i = 0; i < 128; i++) {
				if ((tmp[i] == 1)) {
					cout << "k" << (i + 1) << " ";
				}
			}
			cout << endl;
		}
		it2++;
	}
	cout << endl;

	cout << "even list" << endl;
	it2 = countingBox2.begin();
	while (it2 != countingBox2.end()) {
		if (((*it2).second % 2) == 0) {
			cout << ((*it2).second % 2) << " | " << (*it2).second << "\t";
			bitset<256> tmp = (*it2).first;
			for (int i = 0; i < 128; i++) {
				if ((tmp[i] == 1)) {
					cout << "k" << (i + 1) << " ";
				}
			}
			cout << endl;
		}
		it2++;
	}



}


/***************************************
 * Computing 15 Superpolies with 95 active bits
 ***************************************/
int grain128aSub(int evalNumRounds, int threadNumber) {

	//
	ofstream outputfile, outputfile2;
	outputfile.open("log_grain128a.txt");
	outputfile2.open("log_grain128a2.txt");

	vector<int> cons_pos_vector = { 26,29,30,31,33,40,43,44,45,47,57,58,63,69,71 };

	for (int id = 0; id < cons_pos_vector.size(); id++) {

		int cons_pos = cons_pos_vector[id];

		cout << "///////////////////////////////////////////" << endl;
		cout << "        CONSTANT IV[" << cons_pos + 1 << "]" << endl;
		cout << "///////////////////////////////////////////" << endl;

		// active bits
		vector<int> cube(96, 0);
		for (int i = 0; i < 96; i++) {
			cube[i] = 1;
		}
		cube[cons_pos] = 0;
		cout << "the index of cube" << endl;
		for (int i = 0; i < 96; i++) {
			if (cube[i] == 1)
				cout << "iv" << (i + 1) << ", ";
		}
		cout << endl;


		// constant (secret 3, active 2, const 1, const 0)
		vector<int> flag(256, 0);
		for (int i = 0; i < 128; i++) {
			flag[i] = 3;
		}
		for (int i = 0; i < 96; i++) {
			if (cube[i] == 1) {
				flag[128 + i] = 2;
			}
		}
		for (int i = 96; i < 127; i++) {
			flag[128 + i] = 1;
		}

		//
		map<bitset<256>, int, cmpBitset256> countingBox;
		double dulation = 0;

		//Seperately evaluate the non-linear terms and the linear part of the output bit
		cout << "++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Target s93 + b2 + b15 + b36 + b45 + b64 + b73 + b89" << endl;
		grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 6, { true, 0, });
		if (countingBox.size() == 0) {
			cout << "zero sum\t" << dulation << "sec" << endl;
		}

		cout << "++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Target b12 * s8" << endl;
		grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 5, { true, 0, });
		if (countingBox.size() == 0) {
			cout << "zero sum\t" << dulation << "sec" << endl;
		}

		cout << "++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Target s13 * s20" << endl;
		grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 4, { true, 0, });
		if (countingBox.size() == 0) {
			cout << "zero sum\t" << dulation << "sec" << endl;
		}

		cout << "++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Target b95 * s42" << endl;
		grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 3, { true, 0, });
		if (countingBox.size() == 0) {
			cout << "zero sum\t" << dulation << "sec" << endl;
		}

		cout << "++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Target s60 * s79" << endl;
		grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 2, { true, 0, });
		if (countingBox.size() == 0) {
			cout << "zero sum\t" << dulation << "sec" << endl;
		}

		cout << "++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Target b12 * b95 * s94" << endl;
		grainThreeEnumuration(cube, flag, evalNumRounds, countingBox, dulation, threadNumber, 1, { true, 0, });
		if (countingBox.size() == 0) {
			cout << "zero sum\t" << dulation << "sec" << endl;
		}


		cout << "*****************************" << endl;
		cout << "Final solution" << endl;
		cout << countingBox.size() << " solutions are found" << endl;

		map<bitset<256>, int, cmpBitset256> countingBox2;
		auto it = countingBox.begin();
		while (it != countingBox.end()) {
			bitset<256> tmp = (*it).first;
			for (int i = 0; i < 128; i++) {
				if (cube[i] == 1) {
					tmp[128 + i] = 0;
				}
			}
			for (int i = 96; i < 128; i++) {
				tmp[128 + i] = 0;
			}
			countingBox2[tmp] += (*it).second;
			it++;
		}

		cout << "odd list" << endl;
		auto it2 = countingBox2.begin();
		while (it2 != countingBox2.end()) {
			if (((*it2).second % 2) == 1) {
				cout << ((*it2).second % 2) << " | " << (*it2).second << "\t";
				bitset<256> tmp = (*it2).first;
				for (int i = 0; i < 128; i++) {
					if ((tmp[i] == 1)) {
						cout << "k" << (i + 1) << " ";
					}
				}
				cout << endl;
			}
			it2++;
		}
		cout << endl;

		cout << "even list" << endl;
		it2 = countingBox2.begin();
		while (it2 != countingBox2.end()) {
			if (((*it2).second % 2) == 0) {
				cout << ((*it2).second % 2) << " | " << (*it2).second << "\t";
				bitset<256> tmp = (*it2).first;
				for (int i = 0; i < 128; i++) {
					if ((tmp[i] == 1)) {
						cout << "k" << (i + 1) << " ";
					}
				}
				cout << endl;
			}
			it2++;
		}

		cout << endl;
		cout << endl;

	}

	return 0;






}













/***************************************
 * Code for the practical verification
 ***************************************/
static int roundFuncGrain128a(bitset<128>& b, bitset<128>& s) {

	int f = s[0] ^ s[7] ^ s[38] ^ s[70] ^ s[81] ^ s[96];
	int g = s[0] ^ b[0] ^ b[26] ^ b[56] ^ b[91] ^ b[96] ^ (b[3] & b[67]) ^ (b[11] & b[13]) ^ (b[17] & b[18]) ^ (b[27] & b[59]) ^ (b[40] & b[48]) ^ (b[61] & b[65]) ^ (b[68] & b[84]);
	g ^= (b[88] & b[92] & b[93] & b[95]) ^ (b[22] & b[24] & b[25]) ^ (b[70] & b[78] & b[82]);
	int h = (b[12] & s[8]) ^ (s[13] & s[20]) ^ (b[95] & s[42]) ^ (s[60] & s[79]) ^ (b[12] & b[95] & s[94]);
	int y = h ^ s[93] ^ b[2] ^ b[15] ^ b[36] ^ b[45] ^ b[64] ^ b[73] ^ b[89];

	b >>= 1;
	b[127] = g ^ y;

	s >>= 1;
	s[127] = f ^ y;

	return y;
}
static int encryptionSum(int evalNumRounds, vector<int> cube, vector<int> iv, vector<int> key) {

	bitset<128> b, s;
	for (int i = 0; i < 128; i++) {
		b[i] = key[i];
	}
	for (int i = 0; i < 128; i++) {
		if (iv[i] == 0) s[i] = 0;
		else if (iv[i] == 1) s[i] = 1;
	}

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

		bitset<128> tmp_b = b;
		bitset<128> tmp_s = s;
		for (int i = 0; i < DATA_SIZE; i++) {
			tmp_s[map[i]] = ((in >> i) & 1);
		}

		int z = 0;
		for (int r = 0; r <= evalNumRounds; r++) {
			z = roundFuncGrain128a(tmp_b, tmp_s);
		}
		sum ^= z;
	}

	return sum;
}
static int theoreticalSum(map<bitset<256>, int, cmpBitset256> countingBox, vector<int> cube, vector<int> iv, vector<int> key) {

	int sum = 0;
	auto it = countingBox.begin();
	while (it != countingBox.end()) {

		if (((*it).second % 2) == 1) {

			int var = 1;
			for (int i = 0; i < 256; i++) {
				if ((*it).first[i] == 1) {
					if (i < 128) {
						var *= key[i];
					}
					else {
						if (cube[i - 128] == 0) {
							var *= iv[i - 128];
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
void practicalTestGrain128a(void) {


	//
	ofstream outputfile, outputfile2;
	outputfile.open("log_grain128a.txt");
	outputfile2.open("log_grain128a2.txt");


	srand(time(NULL));

	// create cube index at random
	int numActBits = 1 + rand() % 4;
	cout << numActBits << " active bits" << endl;

	vector<int> cube(96, 0);
	for (int i = 0; i < numActBits; i++) {
		int index;
		do {
			index = (rand() % 96);
		} while (cube[index] == 1);
		cube[index] = 1;
	}
	cube[0] = 1;
	cube[1] = 1;
	cube[2] = 1;


	// fix non-IV bits at random
	// constant (secret 3, active 2, const/undermined 1, const 0)
	vector<int> flag(256, 0);
	for (int i = 0; i < 128; i++) {
		flag[i] = 3;

		if (i < 96) {
			if (cube[i] == 1) {
				flag[128 + i] = 2;
			}
			else {
				flag[128 + i] = 1;
			}
		}
		else if (i == 127) {
			flag[i] = 0;
		}
		else {
			flag[i] = 1;
		}
	}
	


	// input
	for (int i = 0; i < 96; i++) {
		cout << flag[128 + i];
	}
	cout << endl;

	//
	for (int r = 50; r < 120; r++) {
		cout << "##############################" << endl;
		cout << r << " rounds" << endl;

		double dulation;
		map<bitset<256>, int, cmpBitset256> countingBox;
		grainThreeEnumuration(cube, flag, r, countingBox, dulation, 1);

		if (countingBox.size() == 0) {
			cout << "zero sum" << endl;
		}
		else {
			cout << "               key              \t";
			cout << "            iv                ";
			cout << "expe.   ";
			cout << "theo.   ";
			cout << endl;

			for (int trial = 0; trial < 100; trial++) {

				vector<int> key(128);
				for (int i = 0; i < 128; i++)
					key[i] = rand() % 2;

				vector<int> iv(128);
				for (int i = 0; i < 96; i++)
					iv[i] = rand() % 2;
				for (int i = 96; i < 127; i++)
					iv[i] = 1;

				int sum1 = encryptionSum(r, cube, iv, key);
				int sum2 = theoreticalSum(countingBox, cube, iv, key);

				for (int i = 15; i >= 0; i--) {
					int hexvar = 0;
					for (int j = 7; j >= 0; j--) {
						hexvar ^= (key[8 * i + j] << j);
					}
					printf("%02x", hexvar);
				}
				cout << "\t";

				for (int i = 11; i >= 0; i--) {
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
