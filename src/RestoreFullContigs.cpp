/*
 * =====================================================================================
 *
 *       Filename:  RestoreFullContigs.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/07/12 11:38:52
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
#include <set>
#include <vector>
#include <sstream>
#include <fstream>

#include <stdlib.h>

using namespace std;
void printUsage(void);
string process_line(string l, vector <unsigned long> *cs, vector <unsigned long> *efl);

int main(int argc, char* argv[]) {
	// Variables for options.
	char c;
	extern char *optarg;

	//char *SeqFile = NULL;
	//char *SamFile = NULL;
	char *SamFile= NULL;
	char *ContigFile = NULL;

	// All purpose index variables
	
	// Get the different options necessary to run
	while ((c = getopt (argc, argv, "s:c:")) != -1) {
		switch(c) {
			case 's':
				SamFile = optarg;
				break;
			case 'c':
				ContigFile = optarg;
				break;
			default:
				printUsage();
				return 0;
		}
	}

	// If a req'd parameter isn't given
	if ((NULL == ContigFile) || (NULL == SamFile)) { 
		printUsage();
		return -1;
	}

	//map <unsigned int, unsigned long> contig_sizes;
	vector <unsigned long> contig_sizes;
	ifstream contig_fh;
	contig_fh.open(ContigFile);
	if (!contig_fh.is_open()) { 
		cerr << "Couldn't open " << ContigFile << endl;
		return 1;
	}
	
	contig_sizes.push_back(0);
	while (contig_fh.good()) {
		string acc, fasta;
		stringstream tmp;
		unsigned int iacc;
		getline(contig_fh, acc);
		if (!contig_fh.good()) { break; }
		getline(contig_fh, fasta);

		acc.erase(0,1);
		tmp << acc;
		tmp >> iacc;

		if (contig_sizes.size() != iacc) { 
			cerr << "Bad contig accession number: "<< acc << endl;
			return 1;
		}
		contig_sizes.push_back(fasta.size());
		//cout << contig_sizes.size() << " " << iacc << endl;

		//contig_sizes.insert(pair<unsigned int,unsigned long>(iacc,fasta.size()));	
	}

	contig_fh.close();

	ifstream sam_fh;
	sam_fh.open(SamFile);
	if (!sam_fh.is_open()) {
		cerr << "Couldn't open " << SamFile << endl;
		return 2;
	}
	vector <unsigned long> end_frag_lengths;
	end_frag_lengths.resize(contig_sizes.size(), 0);
	while (sam_fh.good()) {
		string line1, line2;
		getline(sam_fh, line1);
		if (!sam_fh.good()) { break; }
		if (line1[0] == '@') { 
			stringstream tmp, tmp2, tmp3;
			string entry1, entry2, entry3;
			unsigned long length;
			unsigned int contig;
			tmp << line1;
			tmp >> entry1;
			if (entry1 == "@SQ") {
				tmp >> entry2;
				tmp >> entry3;
				entry2.erase(0,3);
				entry3.erase(0,3);
				if ((entry2[entry2.size() - 1]) == 'b') {
					entry2.erase(entry2.size() -2 , 2);
					contig =atoi(entry2.c_str());
					length = atol(entry3.c_str());
					end_frag_lengths[contig] = length;
					continue;
				}
				else if ((entry2[entry2.size() - 1]) == 'a') {
					entry2.erase(entry2.size() -2 , 2);
					contig =atoi(entry2.c_str());
				}
				else {
					contig =atoi(entry2.c_str());
				}
				//map <string, unsigned long> cs_it;
				//cs_it = contig_sizes.find(entry2);
				//if (cs_it == contig_sizes.end()) { 
			//		cerr << "Couldn't find contig " << entry2 << " in contig file!" << endl;
		//			return 3;
		//		}
		//		unsigned long len = cs_it->second;
				unsigned long len = contig_sizes[contig];
				if (len == 0) { 
					cerr << "Contig not found for line: " << line1 << endl;
					return 1; 
				}
				cout << "@SQ\tSN:" << contig << "\tLN:" << len << endl;


			}
			else { cout << line1 << endl; }
		}
		else { 
			getline(sam_fh, line2);
			line1 = process_line(line1, &contig_sizes, &end_frag_lengths);
			line2 = process_line(line2, &contig_sizes, &end_frag_lengths);
			//if ((line1[0] != 'X') && (line2[0] != 'X')) {
				cout << line1 << endl << line2 << endl;
			//}
		}
	}
}

string process_line(string l, vector <unsigned long> *cs, vector <unsigned long> *efl) { 
	unsigned int contig;
	unsigned int start_pos;

	size_t i = l.find('.');
	if (i != string::npos) { 
		if (l[i+1] == 'a') { 
			l.erase(i,2);
			return l;
		}
		else if (l[i+1] == 'b') {
			stringstream tmp, tmp2, tmp3, return_ss;
			string entry1, entry2, entry3, entry4;
			tmp << l;
			tmp >> entry1;
			tmp >> entry2;
			tmp >> entry3;
			tmp >> entry4;
			
		//	cout << entry1 << " " << entry2 << " " << entry3 <<  " " << entry4 << endl;

			entry3.erase(entry3.size()-2,2);
			contig = atoi(entry3.c_str());
			start_pos = atoi(entry4.c_str());
		//	cout << entry1 << " " << entry2 << " " << entry3 <<  " " << entry4 << " " << contig << " " << start_pos << endl;
	
	
			unsigned long frag_length  = (*efl)[contig];
			unsigned long contig_length = (*cs)[contig];
			if (!frag_length) { 
				cerr << "bad fragment length for contig: " << contig << endl;
				cerr << "Bad line: " << l << endl;
				throw "bad";
			}
	
			unsigned long new_coord = (contig_length - frag_length) + start_pos;
	
	
			return_ss << entry1 << "\t" << entry2 << "\t" << contig << "\t" << new_coord;
			l = return_ss.str();
		}
	}

	return l;
}

void printUsage(void) { 
	cerr << "Usage: RestoreFullContigs -s <sam file> -c <contig file>" << endl << endl;
}
