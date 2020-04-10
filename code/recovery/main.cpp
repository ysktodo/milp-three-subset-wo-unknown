#include"main.h"

int main(int argc, char const* argv[]){
	
  int target = 0;
  int evalNumRounds = 0;
  int threadNumber = 1;
	int practical = 0;
	int subcube = 0;

  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-r")) evalNumRounds = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "-t")) threadNumber = atoi(argv[i + 1]);

    if (!strcmp(argv[i], "-trivium")) target = 1;
    if (!strcmp(argv[i], "-grain")) target = 2;

		if (!strcmp(argv[i], "-practical")) practical = 1;

		if (!strcmp(argv[i], "-subcube")) subcube = 1;

  }

  cerr << endl;
  if (target == 1) {
		if (practical) {
			cerr << "Practical verification for trivium." << endl;
		}
		else {
			cerr << evalNumRounds << " round trivium." << endl;
		}
  }
  else if (target == 2) {
		if (practical) {
			cerr << "Practical verification for Grain-128AEAD." << endl;
		}
		else {
			cerr << evalNumRounds << " round Grain128a." << endl;
		}
  }
  else {
    cerr << "Please set option " << endl;
    cerr << "  -trivium for Trivium" << endl;
    cerr << "  -grain for Grain128a" << endl;
    return 0;
  }

	if ((subcube == 1) && (target == 1)) {
		cerr << "Sorry, subcube option only works in the application to Grain." << endl;
	}

  if (evalNumRounds == 0) {
		if (practical == 0) {
			cerr << "Please set option about number of rounds as '-r [number of rounds]'" << endl;
			return 0;
		}
  }
	if(practical == 0)
	  cerr << evalNumRounds << " cores are used. to change, set '-r [number of rounds]'" << endl;
  cerr << threadNumber << " cores are used. to change, set '-t [number of threads]'" << endl;

  


  if (target == 1) {

		if (practical) {
			practicalTestTrivium();
		}
		else {
			trivium(evalNumRounds, threadNumber);
		}

  }else if (target == 2) {

		if (practical) {
			practicalTestGrain128a();
		}
		else if (subcube) {
			grain128aSub(evalNumRounds, threadNumber);
		}
		else {
			grain128a(evalNumRounds, threadNumber);
		}

  }

  return 0;
}


