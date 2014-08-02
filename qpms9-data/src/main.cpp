#include <string>
#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <ctime>
#include "RandomGen.h"

using namespace std;

string dnaAlphabet = "ACGT";
string proteinAlphabet = "ACDEFGHIKLMNPQRSTVWY";

void parseIntOrExit(char option, char *argument, int *result) {
	if (1 != sscanf(argument, "%d", result)) {
		cerr << "Option -" << option << " requires integer argument" << endl;
		exit(1);
	}
}

void parseFloatOrExit(char option, char *argument, float *result) {
	if (1 != sscanf(argument, "%f", result)) {
		cerr << "Option -" << option << " requires numerical argument" << endl;
		exit(1);
	}
}

bool hasRepeatedChars(string& s) {
	for (int i = s.size() - 1; i >= 0; --i)
		if (string::npos != s.find(s[i], i + 1))
			return true;
	return false;
}

int main(int argc, char**argv) {
	if (argc == 1) {
		cout
				<< "Arguments: -o <outputFile> -n <nStrings> -m <stringLen> -l <motifLen> -d <motifMutations> -q <quorumPercent> -s <aSize> -a <alphabet> -r <randSeed>"
				<< endl;
		cout << "Where:" << endl;
		cout << "  -o: outputFile [stdout if not specified]" << endl;
		cout << "  -n: number of strings to generate [default 20]" << endl;
		cout << "  -m: length of each string [default 600]" << endl;
		cout << "  -l: length of motif" << endl;
		cout << "  -d: max changes for planted instances of the motif" << endl;
		cout
				<< "  -q: quorum percent - in what percentage of strings to plant motif [default 100]"
				<< endl;
		cout
				<< "  -s: alphabet size - 4 for DNA, 20 for protein, or specify arbitrary alphabet using -alphabet [default 4]"
				<< endl;
		cout << "  -a: alphabet - specify arbitrary alphabet" << endl;
		cout
				<< "  -r: seed for the random number generator; use the same seed for repeatable output"
				<< endl;
		exit(0);
	}

	int n = 20;
	int m = 600;
	int l = -1, d = -1;
	float q = 100;
	int s = 4;
	string alphabet;
	string outputFile;
	int r = time(NULL);

	for (int c; (c = getopt(argc, argv, "o:n:m:l:d:q:s:a:r:")) != -1;)
		switch (c) {
		case 'o':
			outputFile = string(optarg);
			break;
		case 'n':
			parseIntOrExit(c, optarg, &n);
			break;
		case 'm':
			parseIntOrExit(c, optarg, &m);
			break;
		case 'l':
			parseIntOrExit(c, optarg, &l);
			break;
		case 'd':
			parseIntOrExit(c, optarg, &d);
			break;
		case 'q':
			parseFloatOrExit(c, optarg, &q);
			if (q < 0 || q > 100) {
				cerr << "Quorum -q has to be between 0 and 100" << endl;
				exit(1);
			}
			break;
		case 's':
			parseIntOrExit(c, optarg, &s);
			break;
		case 'a':
			alphabet = string(optarg);
			if (hasRepeatedChars(alphabet)) {
				cerr << "Alphabet cannot have repeated characters" << endl;
				exit(1);
			}
			break;
		case 'r':
			parseIntOrExit(c, optarg, &r);
			break;
		case '?':
			cerr << "Unknown option `-" << optopt
					<< "' or missing required argument." << endl;
			exit(1);
		default:
			exit(1);
		}
	if (l < 0 || d < 0) {
		cerr << "Arguments -l and -d must be specified" << endl;
		exit(0);
	}

	if (alphabet.size() > 0) {
		s = alphabet.size();
	} else if (s == 4)
		alphabet = dnaAlphabet;
	else if (s == 20)
		alphabet = proteinAlphabet;
	else {
		cerr << "Alphabet size " << s
				<< " not recognized. Use 4 for DNA, 20 for protein, or specify a custom alphabet using -a."
				<< endl;
		exit(0);
	}

	if (!outputFile.empty()) {
		freopen((const char *) outputFile.c_str(), "w", (FILE*)stdout);
	}

	RandomGen rg(r);

	int toModify = q * n / 100;
	bool shouldPlant[n];
	rg.generateDistinctPositions(toModify, n, shouldPlant);

	string motif = rg.generateString(l, alphabet);
	for (int i = 0; i < n; ++i) {
		string str = rg.generateString(m, alphabet);
		if (shouldPlant[i]) {
			string mutatedMotif = rg.mutate(motif, d, alphabet);
			int p = rg.randomPos(m - l + 1);
			str.replace(p, l, mutatedMotif);
			cout << "> " << i << " Motif " << motif << " planted as "
					<< mutatedMotif << " at position " << p << "." << endl;
		} else {
			cout << "> " << i << " No motif planted." << endl;
		}
		cout << str << endl;
	}

	fclose((FILE*)stdout);
	return 0;
}
