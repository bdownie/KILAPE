//============================================================================
// Name        : BuildScaffolds.cpp
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
string rev_comp(string);

int main(int argc, char* argv[]) {
	if (argc != 4) { printUsage(); return 0; }
	//if (argc != 5) { printUsage(); return 0; }
	
	unsigned int max_contigs = 4096000;
	//int avg_insert = atoi(argv[4]);

	ifstream fhGraph, fhContigs, fhScaff;
	fhGraph.open(argv[2]);

	if (!fhGraph.is_open()) {
		cerr  << "Couldn't open file: " << argv[1] << endl;
		return 0;
	}

	vector<unsigned long> contig_sizes;
	contig_sizes.resize(max_contigs, 0);

//	cerr << "Processing graph file..\n" << endl;

//	unsigned long numLines = 0;
//	int lastSize = 0;
//	cerr << "Progress: ";

	fhContigs.open(argv[1]);
	if (!fhContigs.is_open()) {
		cerr  << "Couldn't open file: " << argv[2] << endl;
		return 0;
	}

	map <unsigned int, string> ContigSeqs;
	map <unsigned int, string>::iterator ContigSeqs_it;
	int contig_count = 0;

	while (fhContigs.good()) {
		string Line1, Line2;
		getline (fhContigs, Line1);
		if (!fhContigs.good()) { break; }	
		getline (fhContigs, Line2);

		Line1.erase(0,1);
		unsigned int contig = atoi(Line1.c_str());

		ContigSeqs.insert(pair <unsigned int, string> (contig, Line2));
		contig_count++;
	}
	if (contig_count == 1) { 
		for (ContigSeqs_it = ContigSeqs.begin(); ContigSeqs_it != ContigSeqs.end(); ContigSeqs_it++) { 
			cout << ">" << ContigSeqs_it->first << endl << ContigSeqs_it->second << endl;
		}
		return 0;
	}

	map <string, int> distances;
	map <string, int>::iterator distances_it;
	//map <string, unsigned int> aligned;
	//map <string, unsigned int>::iterator aligned_it;
	//map <string, unsigned int> orientation;
	//map <string, unsigned int>::iterator orientation_it;

	map <string, unsigned int> first_direction;
	map <string, unsigned int>::iterator first_direction_it;
	map <string, unsigned int> second_direction;
	map <string, unsigned int>::iterator second_direction_it;

	int dist, direct1, direct2;
	//bool reversed;
	string edge;

	while(fhGraph.good()) {
		string Line1;
		stringstream linestream, ls2;
		unsigned int contig1, contig2, count;
		
		getline (fhGraph, Line1);
		if (!fhGraph.good()) { break; }	

		string tag;  
		linestream << Line1;
		linestream >> tag;

		switch (tag[0]) {
			case 'e':
				linestream >> contig1;
				linestream >> contig2;
				linestream >> count;

				ls2 << contig1 << "-" << contig2;
				ls2 >> edge;


				break;
			case 'd':
				linestream >> dist;
				break;
			case '1':
				linestream >> direct1;
				break;
			case '2':
				linestream >> direct2;
				if (!edge.empty()) {
					distances.insert(pair <string, int> (edge, dist));	
					first_direction.insert(pair <string, unsigned int> (edge, direct1));	
					second_direction.insert(pair <string, unsigned int> (edge, direct2));	

				//	cout << edge << " " << dist << " " << direct1 << " " << direct2 << endl;
				}
				edge.clear();
				break;
			default:
				cerr << "Flag " << tag << " not recognized." << endl;
				break;
		}
			

	}

		
		//cerr << endl << endl;
	fhScaff.open(argv[3]);
	if (!fhScaff.is_open()) {
		cerr  << "Couldn't open file: " << argv[3] << endl;
		return 0;
	}

	string Line;
	vector<string> vScaffold;
	getline (fhScaff, Line);
	//cerr << "(" << Line << ")" << endl;
	//if (Line[0] == '\0') { 
	if (Line.empty()) { 
		//cerr << "here" << endl;
		//cerr << "Inside" << endl;
		//cerr << ContigSeqs.size() << endl;
		//cerr << (*ContigSeqs[2305]).first << endl;
		//for (ContigSeqs_it == ContigSeqs.begin(); ContigSeqs_it != ContigSeqs.end(); ContigSeqs_it++) {
		//ContigSeqs_it == ContigSeqs.begin();
		//cerr << (*ContigSeqs_it).first << endl;
			//cout << ">" << (*ContigSeqs_it).first << endl << (*ContigSeqs_it).second << endl;
		//}
		return 1;
	}


/*
	for (first_direction_it = first_direction.begin(); first_direction_it != first_direction.end(); first_direction_it++) { 
		cout << first_direction_it->first << " " << first_direction_it->second << endl;
	}
	for (second_direction_it = second_direction.begin(); second_direction_it != second_direction.end(); second_direction_it++) { 
		cout << second_direction_it->first << " " << second_direction_it->second << endl;
	}
	*/
//	return 0;
		
	boost::split(vScaffold, Line, boost::is_any_of(" \t"));

	string FinalSequence;
	//cerr << Line << endl; 
	unsigned int iContig = atoi(vScaffold[0].c_str());
	//unsigned int iSecondContig = 0;

	//if (vScaffold.size() > 2) { iSecondContig = vScaffold[1]; }

	ContigSeqs_it = ContigSeqs.find(iContig);
	FinalSequence = (*ContigSeqs_it).second;



	bool FirstContig = 1;
	bool LastReversed = 0;
	bool reverse_next = 0;
	int	first_direction_val = 0; 
	int second_direction_val = 0; 

	for (unsigned int i = 1; i <= vScaffold.size(); i++) { 
		string Sequence, edge;
		stringstream tmp;
		iContig = atoi(vScaffold[i-1].c_str());
		ContigSeqs_it = ContigSeqs.find(iContig);
		unsigned int contig1, contig2;


		if (i != vScaffold.size()) { 
			contig1 = atoi(vScaffold[i-1].c_str());
			contig2 = atoi(vScaffold[i].c_str());
			if (contig1 < contig2) { 
				tmp << contig1 << "-" << contig2;
				tmp >> edge;
				first_direction_it = first_direction.find(edge);
				//first_direction_val = first_direction_it->second;
				second_direction_it = second_direction.find(edge);
				//second_direction_val = second_direction_it->second;
				if ((first_direction_it == first_direction.end()) || second_direction_it == second_direction.end()) {
					throw "Bad edge";
				}
			}
			else { 
				tmp << contig2 << "-" << contig1;
				tmp >> edge;
				// Switching first and second here is necessary to maintain the ordering
				// information of the contigs
				//cout << edge << endl;
				first_direction_it = second_direction.find(edge);
				//first_direction_val = second_direction_it->second;
				second_direction_it = first_direction.find(edge);
				//second_direction_val = first_direction_it->second;
				//cout << first_direction_it->first << first_direction_it->second << endl; 
				//cout << second_direction_it->first << second_direction_it->second << endl; 
				if ((first_direction_it == first_direction.end()) || second_direction_it == second_direction.end()) {
					throw "Bad edge";
				}
			}

		//	first_direction_it = first_direction.find(edge);
		//	second_direction_it = second_direction.find(edge);
			first_direction_val = first_direction_it->second;
			second_direction_val = second_direction_it->second;
			
		//	cout << edge << endl;
		//	cout << first_direction_val << " " << second_direction_val << endl;
		//	cout << first_direction_it->first << " " << first_direction_it->second << " " << first_direction_val << endl; 
		//	cout << second_direction_it->first << " " << second_direction_it->second << " " << second_direction_val<< endl; 
		//	cout << "Processing: " << iContig << endl;




	
	
			if (FirstContig) {
		//		cout << "First contig: " << contig1 << " " <<  first_direction_val << endl;
				if (first_direction_val == 0) {
					FinalSequence = rev_comp(FinalSequence);
					LastReversed = 1;
		//			cout << "First Contig: Reversing " << iContig << endl;
				}
		//		else { cout << "First Contig: keeping orient " << iContig << endl; }
				//if (second_direction_val != 0) { 
				//	reverse_next = 1;
				//}
				//else { reverse_next = 0; }
			//cerr <<contig1 << "\t" <<  first_direction_val << endl;
				FirstContig = 0;
			}
			else if (first_direction_val == 0) {
				Sequence = rev_comp((*ContigSeqs_it).second);  
		//		cout << "Reversing " << iContig << endl;
			}
			else { 
				Sequence = (*ContigSeqs_it).second;  
		//		cout << "Keeping orient " << iContig << endl;
			}
		//	cout << "Seq Length " << Sequence.size() << endl;
		//	cerr <<contig1 << "\t" <<  first_direction_val << endl;
		//	cerr << contig2 << "\t" << second_direction_val << endl;
		/*
			else { 
				if (reverse_next) { 
					Sequence = rev_comp((*ContigSeqs_it).second); 
					cout << "Reversing " << iContig << endl;
					if (second_direction_val != 0) { 
						reverse_next = 0;
					}
					else { reverse_next = 1; }
				}
				else { 
					Sequence = (*ContigSeqs_it).second;  
					if (second_direction_val != 0) { 
						reverse_next = 1;
					}
					else { reverse_next = 0; }
					cout << "keeping orient " << iContig << endl; 
				}
			}
			*/
			distances_it = distances.find(edge);
		
		//cerr << edge << endl << (*distances_it).second<<endl;
		//int gap = avg_insert - (*distances_it).second;
			int gap = (*distances_it).second;
		//	cout << FinalSequence.size() << " " << edge << " " << gap<< endl;
			//cout << "edge " << edge << " " << gap << endl;
			if (gap > 0) { 
				for (int j = 0; j <= gap; j++) { 
					FinalSequence.push_back('N');
				}
				FinalSequence.append(Sequence);
			}
			else {
				int j = ((0 - gap) / 2) + 2;
				//cerr << j << endl << gap << endl;
				long len = FinalSequence.size();
				//if (len == 0) { FinalSequence = Sequence; }
				if (len < j) { len = j; }
				FinalSequence.erase(len - j, j);
				Sequence.erase(0,j);
				FinalSequence.append("NN");
				FinalSequence.append(Sequence);
			}
		}
		else { 
			Sequence = (*ContigSeqs_it).second;  
			if (second_direction_val != 0) {
				Sequence = rev_comp(Sequence);
				//cout << "Reversing " << iContig << endl;
			}
			//else { cout << "not reversing " << iContig << endl; }
			FinalSequence.append(Sequence);
		}





		//cerr << contig1 << "\t" << contig2 << endl;

		//aligned_it = aligned.find(edge);
		//int aligned_val = (*aligned_it).second;

		//orientation_it = orientation.find(edge);
		//int orientation_val = (*orientation_it).second;

		//if (second_direction_val != 0) { 
	//		FinalSequence = rev_comp(FinalSequence); 
	//		if (LastOrient == 1) { LastOrient = 0; }
	//		else { LastOrient = 1; }
	//	}

		//cerr << LastOrient << "\t" << aligned_val << endl;

	//	if (LastOrient == 1) { 
	//		if (aligned_val == 0) { 
	//			Sequence = rev_comp((*ContigSeqs_it).second);
	//			LastOrient = 0;
	//		}
	//		else { Sequence = (*ContigSeqs_it).second; }
	//	}
	//	// Last Orient == 0
	//	else if (aligned_val == 0) { 
	//		Sequence = (*ContigSeqs_it).second; 
	//		LastOrient = 1; 
	//	}
	//	else { 
	//		Sequence = rev_comp((*ContigSeqs_it).second);
	//	}

	}

//cout << FinalSequence.size() << endl;
	cout << ">Scaffold_" << vScaffold[0] << endl << FinalSequence << endl;



	return 0;
}

void printUsage() { 
	cerr << "Usage: BuildScaffolds <contig file> <graph file> <scaffold file> " << endl;
	//cerr << "Usage: BuildScaffolds <contig file> <graph file> <scaffold file> <avg insert>" << endl;
}

string rev_comp(string strNucleotides) { 
	string strReverseNucleotides(strNucleotides.size(),'N');
	
	for (int i = strNucleotides.size()-1, j = 0; i >= 0; i--,j++) {
		switch (toupper(strNucleotides[i])) { 
			case 'A' :
				strReverseNucleotides[j] = 'T';
				break;
			case 'C' :
				strReverseNucleotides[j] = 'G';
				break;
			case 'G' :
				strReverseNucleotides[j] = 'C';
				break;
			case 'T' :
				strReverseNucleotides[j] = 'A';
				break;
			default :
				strReverseNucleotides[j] = 'N';
				break;
		}
	}
	return strReverseNucleotides;
}

