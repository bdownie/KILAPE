//============================================================================
// Name        : BuildGraphFromSAM.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <fstream>

#include <sys/stat.h>
#include <getopt.h>


#include <boost/unordered_map.hpp>

#include <boost/algorithm/string.hpp>

using namespace std;
void printUsage();
bool check_sam_multihit (string l);

//#typedef boost::unordered_map<string, string> hashmap;

int main(int argc, char* argv[]) {
	char c;
	char *SamFile = NULL;
	char *ContigFile =  NULL;

	unsigned int insert_size = 0;
	unsigned int insert_sd = 0;

	bool soap = 0;

	while ((c = getopt (argc, argv, "s:c:i:d:o")) != -1) {
		switch(c) {
			case 's':
				SamFile = optarg;
				break;
			case 'c':
				ContigFile = optarg;
				break;
			case 'i':
				insert_size = atoi(optarg);
				break;
			case 'd':
				insert_sd = atoi(optarg);
				break;
			case 'o':
				soap = 1;
				break;
			case 'h':
				printUsage();
				return 0;
			default:
				printUsage();
				return 0;
		}
	}
	istream *in;
	ifstream fhSAM, fhContig;
	if (ContigFile == NULL) { 
		printUsage();
		return 0;
	}
		
	if (insert_size && !insert_sd) { 
		printUsage();
		return 0;
	}
//	if (argc != 5) { printUsage(); return 0; }
	
	if (SamFile == NULL) { 
		in = &cin;
	}
	else { 
		fhSAM.open(SamFile);
		if (!fhSAM.is_open()) {
			cerr  << "Couldn't open file: " << SamFile << endl;
			return 0;
		}
		in = &fhSAM;
	}
	fhContig.open(ContigFile);
	//unsigned int insert_size = atoi(argv[2]);

	if (!fhContig.is_open()) {
		cerr  << "Couldn't open file: " << ContigFile << endl;
		return 0;
	}

	// Modified 28.8.2011 BRD
	//unsigned int max_size = insert_size + (insert_sd * 3);
	unsigned int max_size = insert_size + insert_sd;

	//unsigned int max_size = insert_size;
	//cout << max_size << endl; return 0;

	//boost::unordered_map<unsigned int, unsigned int> contig_sizes;
	boost::unordered_map<string, unsigned int> edges;
	boost::unordered_map<string, unsigned int>::iterator it;

	boost::unordered_map<string, unsigned int> edges_dist;
	boost::unordered_map<string, unsigned int>::iterator dist_it;
	//boost::unordered_map<unsigned int, unsigned int>::iterator uit;

	boost::unordered_map<string, unsigned int> edges_dist_count;
	boost::unordered_map<string, unsigned int>::iterator dist_count_it;

	boost::unordered_map<string, unsigned int> orientation;
	boost::unordered_map<string, unsigned int>::iterator or_it;

	boost::unordered_map<string, unsigned int> contig1_read_direction;
	boost::unordered_map<string, unsigned int>::iterator c1direct_it;
	boost::unordered_map<string, unsigned int> contig2_read_direction;
	boost::unordered_map<string, unsigned int>::iterator c2direct_it;
	//boost::unordered_map<unsigned int, unsigned int> edge;
	//boost::unordered_map<unsigned int, unsigned int>::iterator edge_it;

	//boost::unordered_map<string, int> edge_dist;
	//boost::unordered_map<string, int>::iterator dist_it;

	vector<unsigned int> contig_sizes;
	contig_sizes.resize(1500000, 0);

	//unsigned int min_contig_size = (int) insert_size/2;

//	cerr << "Processing SAM file..\n" << endl;

	unsigned long numLines = 0;
	int lastSize = 0;

	while (fhContig.good()) {
		string Line1, Line2;
		stringstream edge;
		vector<string> vLine1, vLine2;
		
		getline (fhContig, Line1);
		getline (fhContig, Line2);
		
		if (!fhContig.good()) { break; }

		Line1.erase(0,1);
		unsigned int contig = atoi(Line1.c_str());
		if (contig_sizes.capacity() < contig) { 
			contig_sizes.resize(contig_sizes.capacity() + 500000, 0);
		}
		contig_sizes[contig] = Line2.size();
	}

	ofstream fhDebug;
	stringstream debug_string_stream;
	string debug_string;
	//debug_string.assign("graph_edges");
	debug_string_stream << "graph_edges.";
	debug_string_stream << insert_size;
	debug_string_stream << ".sam";
	debug_string_stream >> debug_string;


	fhDebug.open(debug_string.c_str());

	while (in->good()) {
		string Line1, Line2;
		string ContigPairs;
		stringstream edge;
		vector<string> vLine1, vLine2;
		
		getline (*in, Line1);
		
		if (Line1[0] == '@') { 
			if ((Line1[1] == 'S') && (Line1[2] == 'Q'))   {
				boost::split(vLine1, Line1, boost::is_any_of(" \t"));
				vLine1[1].erase(0,3);
				vLine1[2].erase(0,3);
				unsigned int contig, size;
				contig = atoi(vLine1[1].c_str());
				size = atoi(vLine1[2].c_str());
				contig_sizes[contig] = size;
			}
			continue;
		}

		getline (*in, Line2);
		numLines +=2;
		if (!in->good()) { break; }

		if ((numLines % 100000) == 0) {
			for (int k = 0; k < lastSize; k++) {
				//cerr << "\b";
			}
			stringstream Output;
			string strOutput;
			Output << numLines;
			Output >> strOutput;
			lastSize = strOutput.size();
			//cerr << numLines;
		}
		stringstream ssLine1, ssLine2;
		string trash;
		ssLine1 << Line1;
		ssLine2 << Line2;
//cout << endl;
		bool direct1, direct2;
		int flag1, flag2;
		string contig1, contig2;
		unsigned int start1, start2;
		//unsigned int Seq1Length, Seq2Length;
		unsigned int acc1, acc2;

		ssLine1 >> acc1;
		ssLine2 >> acc2;
		if (acc2 != (acc1 + 1)) { cerr << "accessions of paired reads don't match!" << endl << Line1 << endl << Line2 << endl; return 0; }

		//boost::split(vLine1, Line1, boost::is_any_of(" \t"));
		//boost::split(vLine2, Line2, boost::is_any_of(" \t"));


		//unsigned int flag1, flag2;	
		ssLine1 >> flag1;
		ssLine2 >> flag2;

		//string contig1, contig2;
		ssLine1 >> contig1;
		ssLine2 >> contig2;


		//unsigned int start1, start2;
		ssLine1 >> start1;
		ssLine2 >> start2;


		string tmp = Line1;
		bool skip_lines = 0;	
		skip_lines = check_sam_multihit(tmp);
		if (skip_lines) { continue; }
		tmp = Line2;
		skip_lines = check_sam_multihit(tmp);
		if (skip_lines) { continue; }

				/*		

		// Field phred score
		ssLine1 >> trash;
		ssLine2 >> trash;

		// Field CIGAR
		ssLine1 >> trash;
		ssLine2 >> trash;

		// Field 6
		ssLine1 >> trash;
		ssLine2 >> trash;

		// Field 7
		ssLine1 >> trash;
		ssLine2 >> trash;

		// Field 8
		ssLine1 >> trash;
		ssLine2 >> trash;

		// Field sequence
		string Seq1, Seq2;
		ssLine1 >> Seq1;
		ssLine2 >> Seq2;
		Seq1Length = Seq1.length();
		Seq2Length = Seq2.length();

		// Field 9
		ssLine1 >> trash;
		ssLine2 >> trash;

		// Field 10
		ssLine1 >> trash;
		ssLine2 >> trash;


		// Field 11
		ssLine1 >> trash;
		ssLine2 >> trash;

		// Field 12
		ssLine1 >> trash;
		ssLine2 >> trash;

		//contig1 = vLine1[2];
		//contig2 = vLine2[2];
		*/
	
		if ((contig1[0] ==  '*') || (contig2[0] == '*')) {
			continue;
		}

		if (contig1 == contig2) { 
			continue;
		}

		//cout << Line1 << endl << Line2 << endl;
		//	cout << endl<< Seq1Length << " " << Seq2Length << endl;
		//	cout << flag1 << " " << flag2 << endl;
		//	cout << contig1 << " " << contig2 << endl;
		//	cout << start1 << " " << start2 << endl;
		unsigned long iContig1, iContig2;
		iContig1 = atol(contig1.c_str());
		iContig2 = atol(contig2.c_str());

	//	if ((contig_sizes[iContig1] < min_contig_size) || (contig_sizes[iContig2] < min_contig_size)) { 
	//		continue;
	//	}
		//if (iContig1 == iContig2) { 
	//		continue;
	//	}
		//start1 = atoi(vLine1[3].c_str());
		//start2 = atoi(vLine2[3].c_str());

		//bool direct1, direct2;
		direct1 = flag1 & 16 ? 0 : 1;
		direct2 = flag2 & 16 ? 0 : 1;

		int distance = 0;

		fhDebug << Line1 << endl << Line2 << endl;
		if (direct1) { 
			distance += contig_sizes[iContig1] - start1;
		//	cout << "direct1: " << contig_sizes[iContig1] << " " << start1 << endl;
		}
		else { distance += start1; 
		//	cout << "not direct1: " << contig_sizes[iContig1] << " " << start1 << endl;
		}
		//cout << endl << iContig1 << " " << contig_sizes[iContig1]  <<  endl<< direct1 << " " << distance << endl;
		
		if (direct2) { 
			distance += contig_sizes[iContig2] - start2;
		//	cout << "direct2: " << contig_sizes[iContig2] << " " << start2 << endl;
		}
		else { distance += start2; 
		//	cout << "not direct2: " << contig_sizes[iContig2] << " " << start2 << endl;
		}
		//cout << distance << endl;
		//cout << iContig2 << " " << contig_sizes[iContig2]  << endl << direct2 << " " << distance << endl;

		bool FirstContigReadDirection;
		bool SecondContigReadDirection;
		//distance += Seq1Length;
		//distance += Seq2Length;
		//cout << distance << endl;


//cout << endl << "Start:" << endl << Line1 << endl << Line2 << endl <<  "(" <<  distance << " " << max_size << ")" << endl;
///cout <<  "(" <<  distance << " " << max_size << ")" << endl;
		if (max_size && (distance <= max_size)) { 
			if (iContig1 < iContig2) { 
				edge << iContig1 << "-" << iContig2;
				FirstContigReadDirection = direct1;
				SecondContigReadDirection = direct2;
			}
			else {
				edge << iContig2 << "-" << iContig1;
				FirstContigReadDirection = direct2;
				SecondContigReadDirection = direct1;
			}
	
	
	
			edge >> ContigPairs;
			ContigPairs[ContigPairs.find('-')] = ' ';
	
			c1direct_it = contig1_read_direction.find(ContigPairs);
			if (c1direct_it == contig1_read_direction.end()) {
				contig1_read_direction.insert(pair<string, unsigned int> (ContigPairs, 0));
				c1direct_it = contig1_read_direction.find(ContigPairs);
			}
			c2direct_it = contig2_read_direction.find(ContigPairs);
			if (c2direct_it == contig2_read_direction.end()) {
				contig2_read_direction.insert(pair<string, unsigned int> (ContigPairs, 0));
				c2direct_it = contig2_read_direction.find(ContigPairs);
			}
			//cerr << endl<< Line1<< endl << Line2 << endl;
			//cerr << ContigPairs<< endl;
		//	cerr << direct1 << "\t" << direct2 << "\t" << endl;
			if (FirstContigReadDirection) { (*c1direct_it).second++; }
			if (SecondContigReadDirection) { (*c2direct_it).second++; }
			//cerr << FirstContigReadDirection << "\t" << SecondContigReadDirection << "\t" << endl;
	
			//or_it = orientation.find(Line1);
			//if (or_it == orientation.end()) {
			//	orientation.insert(pair<string, unsigned int> (Line1, 0));
			//	or_it = orientation.find(Line1);
			//}
			//if (direct1^direct2) {
			//	(*or_it).second++;
			//}
				
			//	if (ContigPairs == "150 345") { cout << distance << " " << max_size << endl; }
			//cout << ContigPairs << endl;
			//cout << ContigPairs << " " << distance << endl;
			it = edges.find(ContigPairs);
			if (it == edges.end()) { 
				edges.insert(pair <string, unsigned int>(ContigPairs, 1));
				edges_dist.insert(pair <string, unsigned int>(ContigPairs, distance));
				edges_dist_count.insert(pair <string, unsigned int>(ContigPairs, 1));
			}
			else {
				unsigned int val;
				val = (*it).second;
				val++;
				(*it).second = val;
				dist_it = edges_dist.find(ContigPairs);
				val = (*dist_it).second;
				val += distance;
			//if (ContigPairs == "150 345") { cout << val  << endl; }
				(*dist_it).second = val;
	
				dist_count_it = edges_dist_count.find(ContigPairs);
				val = (*dist_count_it).second;
				val++;
				(*dist_count_it).second = val;
			}
		}
	}
	for (int k = 0; k < lastSize; k++) {
		//cerr << "\b";
	}
	//cerr << numLines;
	//cerr << endl;

	//cerr << "Finishing graph..\n" << endl;
	//int avg_dist = 0;
	for (it = edges.begin(); it != edges.end(); it++) {
		string val = (*it).first;
		dist_it = edges_dist.find(val);
		dist_count_it = edges_dist_count.find(val);
		//or_it = orientation.find(val);
		c1direct_it = contig1_read_direction.find(val);
		c2direct_it = contig2_read_direction.find(val);
		unsigned int avg_dist = (unsigned int) (*dist_it).second/(*dist_count_it).second;
		// Modified BRD
		//unsigned int maximum_size = insert_size + insert_sd + (2 * insert_sd/(*dist_count_it).second);
		//unsigned int maximum_size = 2* insert_size;
		if (max_size && (avg_dist < max_size)) { 	
	//	if (((*dist_it).second/(*dist_count_it).second)  < max_size) { 
			cout << "edge: " << (*it).first << "\t"<< (*it).second << endl;
			cout << "dist " << avg_dist << endl;
			cout << "1direct " << (*c1direct_it).second << endl;
			cout << "2direct " << (*c2direct_it).second << endl;
		}
		//}
	}


	return 0;
}

void printUsage() { 
	cerr << "Usage: BuildGraphFromSAM -s <sam file> -c <contig file> -i <insert size> -d <insert std dev> [ -o (use soap mode (not implemented)) ]" << endl;
}

bool check_sam_multihit (string l) {
	stringstream ss;
	ss << l;
	string tmp, tmp2;
	stringstream sstmp;
	int itmp = 0;

	while (ss.good()) { 
		ss >> tmp;
		if ((tmp[0] == 'X') && (tmp[1] == '0')) { 
			itmp = 0;
			tmp.erase(0,5);
			sstmp << tmp;
			sstmp >> itmp;
			if (itmp > 1) { return 1; }
		}
		if ((tmp[0] == 'X') && (tmp[1] == '1')) { 
			itmp = 0;
			tmp.erase(0,5);
			sstmp << tmp;
			sstmp >> itmp;
			if (itmp > 0) { return 1; }
		}
	//	cout << l << endl;
	}
	return 0;
}

