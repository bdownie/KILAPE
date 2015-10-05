/*
 * =====================================================================================
 *
 *       Filename:  filterSamOnCount.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/21/13 13:56:34
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <algorithm>
#include <fstream>

using namespace std;

void printUsage(void);

int main(int argc, char* argv[]) {
	char c;
	char *SamFile = NULL;
	char *CountFile = NULL;
	unsigned int quant_val = 99;

	while ((c = getopt (argc, argv, "s:q:c:h")) != -1) {
		switch(c) {
			case 's':
				SamFile = optarg;
				break;
			case 'c':
				CountFile = optarg;
				break;
			case 'q':
				quant_val = atoi(optarg);
				break;
			case 'h':
				printUsage();
				return 0;
			default:
				printUsage();
				return 0;
		}
	}
	if (quant_val > 100) { cerr << "Quantile value must be integer 0-100" << endl; }

	ifstream SAM;
	ifstream COUNT;
	if ((SamFile == NULL) || (CountFile == NULL)) {
		printUsage();
		return 0;
	}
	SAM.open(SamFile);
	COUNT.open(CountFile);
	if (!SAM.is_open()) { 
		cerr << "Couldn't open sam file: " << SamFile << endl;
		return 1;
	}
	if (!COUNT.is_open()) { 
		cerr << "Couldn't open count file: " << CountFile << endl;
		return 1;
	}


	map <unsigned long,unsigned int> loc_density;
	map <unsigned long,unsigned int>::iterator ld_it;
	vector <unsigned int> numbers;
	while (COUNT.good()) { 
		string line;
		getline(COUNT,line);
		if (!COUNT.good()) { break; }
		stringstream tmp;
		unsigned int chr, loc, number;
		tmp << line;
		tmp >> chr;
		tmp >> loc;
		tmp >> number;

		unsigned long unique_site =  chr;
		unique_site = unique_site<<32;
		unique_site = unique_site | loc;
		loc_density.insert(pair<unsigned long, unsigned int> (unique_site,number));
		numbers.push_back(number);
	}
	COUNT.close();
	sort (numbers.begin(),numbers.end());
//	for (unsigned int i = 0; i <= 20; i++) { 
//		cout << numbers[i] << " ";
//	}
//	cout << endl;
	unsigned int threshold;
	unsigned int index_to_check;
	if (numbers.size() == 0) { threshold = 0; }
	else { 
		index_to_check = int(quant_val * numbers.size() / 100);
		threshold = numbers[--index_to_check];
	}
	cerr << "Threshold determined " << threshold << " for quantile " << quant_val << "%" << endl;
//	cerr << "Index: " << index_to_check << " of " << numbers.size() << endl;

	string trash_line_suffix = "\t4\t*\t0";
	while (SAM.good()) {
		string line;
		getline(SAM,line);
		string out_line = "";
		if (line[0] == '@') { cout << line << endl; }
		else { 
			string trash;
			string chr_string;
			string read_num;
			stringstream tmp;
			unsigned int loc;
	
			tmp << line;
			tmp >> read_num;
			tmp >> trash;
			tmp >> chr_string;
			tmp >> loc;
	
			if (chr_string[0] == '*') { 
				out_line = read_num;
				out_line.append(trash_line_suffix);
			}
			else { 
				unsigned int chr;
				chr = atoi(chr_string.c_str());
	
				unsigned long unique_site =  chr;
				unique_site = unique_site<<32;
				unique_site = unique_site | loc;
				ld_it = loc_density.find(unique_site);
				if (ld_it == loc_density.end()) { 
					out_line = line;
				}
				else if (ld_it->second > threshold) { 
					out_line = read_num;
					out_line.append(trash_line_suffix);
				}
				else { 
					out_line = line;
				}
			}
			cout << out_line << endl;
		}
	}
}


void printUsage() { 
	cerr << "Usage: filterSamOnCount -s <sam file> -c <count file> [-q 0-100]" << endl;
}
