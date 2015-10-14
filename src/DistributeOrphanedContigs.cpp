/*
 * =====================================================================================
 *
 *       Filename:  DistributeOrphanedContigs.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/08/12 13:03:48
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <fstream>

#include <sys/stat.h>
#include <getopt.h>
#include <omp.h>

#include <google/dense_hash_map>
#include <boost/unordered_map.hpp>

#include <boost/algorithm/string.hpp>

#define MAX_BUFFER_SIZE 100000
#define TRUE 1
using namespace std;

// Function declarations
void printUsage();
string get_working_directory (unsigned int LastDir);

unsigned int parse_graph(char *, vector<string> *);
unsigned int parse_scaffold_file(char *, unsigned int *, vector<unsigned int> *, vector<set <unsigned int> > *, vector<string> *);
unsigned int parse_seq(string , vector<string> *);
vector <string> parse_sam (string, string  , vector<unsigned int> *, unsigned int, bool);
unsigned int parse_contig(char *, unsigned int *, vector<unsigned int> *, vector<set <unsigned int> > *, vector<string> *);
typedef boost::unordered_map<string, string> hashmap;

string create_working_directories (unsigned int LastDir);
unsigned long get_last_read_num(string );

int sub_max_threads = 0;
bool fill_gap_only = 0;
bool output_unmatched = 0;
bool no_lib = 0;

int main(int argc, char* argv[]) {
	// Variables for options.
	char c;
	extern char *optarg;

	char *ScaffoldFile = NULL;
	char *OrphanFile = NULL;
	//char *SeqFile = NULL;
	//char *SamFile = NULL;
	char *GraphFile = NULL;

	// All purpose index variables
	//unsigned int i, j; 
	//unsigned long numLines = 0;
	//int lastSize = 0;

	// Get the different options necessary to run
	while ((c = getopt (argc, argv, "i:c:g:h")) != -1) {
		switch(c) {
			case 'i':
				ScaffoldFile = optarg;
				break;
			case 'c':
				OrphanFile = optarg;
				break;
			case 'g':
				GraphFile = optarg;
				break;
			default:
				printUsage();
				return 0;
		}
	}

	// If a req'd parameter isn't given
	if ((NULL == OrphanFile) || (NULL == ScaffoldFile) || (NULL == GraphFile)) { 
		printUsage();
		return 0;
	}

	vector <string> Orphans;
	vector <unsigned int> ContigsToScaffolds;
	vector <string> OutputStrings;
	vector <unsigned int> UsedOrphans;

	ifstream OrphanFH;
	OrphanFH.open(OrphanFile);
	if (!OrphanFH.is_open()) {
		cerr << "Couldn't open file: " << OrphanFile << endl;
		return 1;
	}
	while (OrphanFH.good()) {
		string Acc, Seq;
		getline(OrphanFH,Acc);
		getline(OrphanFH,Seq);
		if (!OrphanFH.good()) { break; }


		string Acc_string;
		Acc_string.assign(Acc);
		Acc_string.append("\n");
		Acc_string.append(Seq);
		Acc_string.append("\n");
		//cout << Acc << endl << Seq << endl << "Here " <<  Acc_string << endl;
		Acc.erase(0,1);
		unsigned int Acc_num = atoi(Acc.c_str());
		if (Orphans.size() <= Acc_num) {
			Orphans.resize(Acc_num+1,"");
			UsedOrphans.resize(Acc_num+1,0);
		}
		Orphans[Acc_num].assign(Acc_string); 
	}

	for (unsigned int i = 1; i < Orphans.size(); i++) { 
//		cout << i << " " << Orphans[i] << endl;
	}
	ifstream ScaffFH;
	ScaffFH.open(ScaffoldFile);
	if (!ScaffFH.is_open()) {
		cerr << "Couldn't open file: " << ScaffoldFile << endl;
		return 1;
	}
	unsigned int scaffold_number = 1;

	while (ScaffFH.good()) {
		string Line;
		stringstream tmp;
		unsigned int val;
		getline(ScaffFH,Line);
		if (!ScaffFH.good()) { break; }

		tmp << Line;
		while (tmp) { 
			tmp >> val;

			if (ContigsToScaffolds.size() <= val) {
				ContigsToScaffolds.resize(val+1);
			}
			ContigsToScaffolds[val] = scaffold_number;
		}
		scaffold_number++;
	}

	OutputStrings.resize(scaffold_number, "");
	vector <unsigned int> OutputStringFilled;
	OutputStringFilled.resize(scaffold_number,0);

	ifstream GraphFH;
	GraphFH.open(GraphFile);
	if (!GraphFH.is_open()) {
		cerr << "Couldn't open file: " << GraphFile << endl;
		return 1;
	}

	while (GraphFH.good()) {
		string Line;
		stringstream tmp;
		getline(GraphFH,Line);
		if (!GraphFH.good()) { break; }
		string param1, param2, param3;
		tmp << Line;
		tmp >> param1;

		if (param1 == "edge:") {
			tmp >> param2;
			tmp >> param3;
			unsigned int iparam2 = atoi(param2.c_str());
			unsigned int iparam3 = atoi(param3.c_str());

			if ((ContigsToScaffolds[iparam2] > 0) && (ContigsToScaffolds[iparam3] == 0)) { 
				OutputStrings[ContigsToScaffolds[iparam2]].append(Orphans[iparam3]);
				OutputStringFilled[ContigsToScaffolds[iparam2]] = 1;
				UsedOrphans[iparam3] = 1;
				//cout << iparam3 << endl;
			}
			else if ((ContigsToScaffolds[iparam2] == 0) && (ContigsToScaffolds[iparam3] > 0)) { 
				OutputStrings[ContigsToScaffolds[iparam3]].append(Orphans[iparam2]);
				OutputStringFilled[ContigsToScaffolds[iparam3]] = 1;
				UsedOrphans[iparam2] = 1;
				//cout << iparam2 << endl;
			}
		}
	}
//	cerr << "starting outputting data.. " << endl;
//	cerr << "Size: " << OutputStrings.size() << endl;

	for (unsigned int Scaffold = 1; Scaffold < OutputStrings.size(); Scaffold++) { 
		if (OutputStringFilled[Scaffold]) { 
			stringstream ddir;
			string dir;
			unsigned int ScaffoldDir1, ScaffoldDir2, ScaffoldDir3, ScaffoldDir4, ScaffoldDir5;

			ScaffoldDir1 = Scaffold % 10;
			ScaffoldDir2 = Scaffold % 100;
			ScaffoldDir3 = Scaffold % 1000;
			ScaffoldDir4 = Scaffold % 10000;
			ScaffoldDir5 = Scaffold % 100000;
	
			ddir.clear();


			ddir << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2 << "/" << ScaffoldDir3 
			   	  << "/" << ScaffoldDir4 << "/" << ScaffoldDir5 << "/" << Scaffold;
			ddir >> dir;
			dir.append("/contigs.orphan.fasta");
//			cout << dir << endl;
			ofstream Orphan_OUT;
			Orphan_OUT.open(dir.c_str());
			if (Orphan_OUT.is_open()) {
				Orphan_OUT << OutputStrings[Scaffold];
			}
			Orphan_OUT.close();
		}
		//cout << "Scaffold: " << Scaffold << endl << "Dir: " << dir << endl;
	}

	cerr << "finished outputting data.. " << endl;

	for (unsigned int i = 1; i < UsedOrphans.size(); i++) { 
		//cout << i << " " << UsedOrphans[i] << endl;
		if ((!UsedOrphans[i]) && (!Orphans[i].empty())) { 
			cout << Orphans[i];
		}
	}
//		cout << i << endl << OutputStrings[i] << endl << "-------" << endl;

}

void printUsage() {
	cerr << "DistributeOrphanedContigs -i <scaffold file> -c <orphan file> -g <graph file>\n"  << endl;
}
