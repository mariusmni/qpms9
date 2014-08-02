/*
 * main.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: marius
 */
#include "MotifWorker.h"
#include <ctime>
#include <cmath>
#include <getopt.h>
#include "utils.h"
#include <cstdio>

#ifndef NOMPI
#include <mpi.h>
#endif

#define MAX_MOTIFS_PER_PROC 100000

void *runWork(void *w) {
	((MotifWorker*) w)->doWork();
	return NULL;
}
#ifndef NOMPI
void createWorkerThread(pthread_t& t, MotifWorker *w) {
	pthread_create(&t, NULL, runWork, (void *) w);
}
#endif

void force_quit() {
#ifndef NOMPI
	MPI::COMM_WORLD.Abort(1);
	MPI::Finalize();
#endif
	exit(1);
}

void printHelp() {
	cout
			<< "Arguments: inputFile -l <motifLength> -d <#Mutations> [-q <quorumPercent>] [-s <stackSz>] [-n nPrime] [-o <outputFile>] [-h]"
			<< endl;
	cout << "Where:" << endl;
	cout << "  -l: length of motif" << endl;
	cout << "  -d: max changes for planted instances of the motif" << endl;
	cout << "  -o: outputFile [stdout if not specified]" << endl;
	cout << "  -q: min percentage of strings that have motif [default 100]"
			<< endl;
	cout
			<< "  -s: (optional) number of lmers for which we generate common neighborhoods"
			<< endl;
	cout
			<< "  -n: (optional) number of sequences for which we initially find motifs"
			<< endl;
	cout << "  -h: print this help" << endl;
}

void parseArgs(int argc, char **argv, MotifConfig& mc, string& inputFile,
		string& outputFile) {
	mc.L = -1;
	mc.d = -1;
	mc.minStackSize = -1;
	mc.nPrime = -1;
	mc.q_percent = 100;
	for (int c; (c = getopt(argc, argv, "o:l:d:q:s:n:h")) != -1;) {
		switch (c) {
		case 'h':
			printHelp();
			force_quit();
			break;
		case 'o':
			outputFile = string(optarg);
			break;
		case 'l':
			parseIntOrExit(c, optarg, &mc.L);
			break;
		case 'd':
			parseIntOrExit(c, optarg, &mc.d);
			break;
		case 'q':
			parseFloatOrExit(c, optarg, &mc.q_percent);
			if (mc.q_percent < 0 || mc.q_percent > 100) {
				cerr << "Quorum percent -q has to be between 0 and 100" << endl;
				exit(1);
			}
			break;
		case 's':
			parseIntOrExit(c, optarg, &mc.minStackSize);
			break;
		case 'n':
			parseIntOrExit(c, optarg, &mc.nPrime);
			break;
		case '?':
			cerr << "Unknown option `-" << optopt
					<< "' or missing required argument." << endl;
			force_quit();
			break;
		default:
			force_quit();
		}
	}

	if (optind < argc)
		inputFile = string(argv[optind]);
	else {
		cerr << "No input file specified!" << endl;
		force_quit();
	}

	if (mc.L < 0 || mc.d < 0) {
		cerr << "Arguments -l and -d must be specified" << endl;
		force_quit();
	}
}

int main(int argc, char **argv) {
	clock_t startTime = clock();

	int myRank = 0;
	int nProcs = 1;
	int totalMsgWords;
	int *b = NULL;

#ifndef NOMPI
	int MotifTransferTag = 2;
	MPI::Status status;
	MPI::Init_thread(argc, argv, MPI::THREAD_FUNNELED);
	myRank = MPI::COMM_WORLD.Get_rank();
	nProcs = MPI::COMM_WORLD.Get_size();
#endif
	if (argc == 1) {
		if (myRank == 0)
		{
			printHelp();
		}
		exit(0);
	}

	string inputFile;
	string outputFile;
	if (myRank == 0) {

		MotifConfig mc;
		parseArgs(argc, argv, mc, inputFile, outputFile);

		b = MotifWorker::readAndEncodeInput(inputFile.c_str(), totalMsgWords, mc);
		if (b == NULL) {
			cerr << "Unable to open file " << inputFile << endl;
			force_quit();
		}
	}

#ifndef NOMPI
	MPI::COMM_WORLD.Bcast(&totalMsgWords, sizeof(int), MPI::CHAR, 0);
	if (myRank != 0) {
		b = new int[totalMsgWords];
	}

	MPI::COMM_WORLD.Bcast(b, totalMsgWords * sizeof(int), MPI::CHAR, 0);
#endif

	MotifWorker worker(myRank, nProcs, startTime, b);
	delete[] b;

	//#ifndef NOMPI
	//	if (myRank != 0) {
	//		cerr << "Processor " << myRank << " received ";
	//		worker.printConfig();
	//	}
	//#endif

	if (myRank == 0 && nProcs > 1) {
#ifndef NOMPI // should have MPI but we avoid pthread also if possible
		pthread_t thread;
		createWorkerThread(thread, &worker);
		worker.schedulerLoop();
		pthread_join(thread, NULL);
		cerr << endl;
#else
		cerr << "This should never happen" << endl;
		force_quit();
#endif
	} else {
		worker.doWork();
	}

	set<MyString> motifs = worker.getMotifs();

	if (myRank == 0) {
#ifndef NOMPI
		int sz;
		int *buf = worker.allocateMotifBuffer(MAX_MOTIFS_PER_PROC, sz);
		for (int i = 1; i < nProcs; ++i) {
			MPI::COMM_WORLD.Recv(buf, sz, MPI::INT, MPI::ANY_SOURCE,
					MotifTransferTag, status);
			worker.decodeMotifs(buf, motifs);
			//			int nm = worker.decodeMotifs(buf, motifs);
			//			cerr << "Processor " << myRank << " received " << nm
			//					<< " motifs from " << status.Get_source() << endl;
		}
		delete[] buf;
#endif
		if (!outputFile.empty())
			freopen((const char *) outputFile.c_str(), "w", (FILE*) stdout);
		worker.printMotifs(motifs);
		if (!outputFile.empty())
			fclose((FILE*) stdout);
		cerr << "Total motifs found: " << motifs.size() << endl;
	} else {
#ifndef NOMPI
		int nMotifs = motifs.size();
		if (nMotifs > MAX_MOTIFS_PER_PROC) {
			cerr << "Processor " << myRank << " found " << motifs.size()
							<< " motifs; Keeping first " << MAX_MOTIFS_PER_PROC
							<< " of them" << endl;
			nMotifs = MAX_MOTIFS_PER_PROC;
		}
		int requiredMem;
		int *buf = worker.encodeMotifs(motifs, nMotifs, requiredMem);
		MPI::COMM_WORLD.Send(buf, requiredMem, MPI::INT, 0, MotifTransferTag);
		//		cerr << "Processor " << myRank << " sent " << nMotifs
		//				<< " motifs to proc 0" << endl;
		delete[] buf;
#endif
	}

	if (myRank == 0) {
		clock_t endTime = clock();
		float seconds = (float) (endTime - startTime) / CLOCKS_PER_SEC;
		cerr << "Time on processor " << myRank << ": " << seconds << "s; file "
				<< inputFile << endl;
	}

#ifndef NOMPI
	MPI::Finalize();
#endif
	return 0;
}

