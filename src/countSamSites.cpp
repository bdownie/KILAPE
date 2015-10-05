/*
 * =====================================================================================
 *
 *       Filename:  countSamSites.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/20/13 09:25:39
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
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
	char *CountFile = NULL;
	unsigned int quant_val = 99;
	char c;

	while ((c = getopt (argc, argv, "c:")) != -1) {
		switch(c) {
			case 'c':
				CountFile = optarg;
				break;
		}
	}



	map <unsigned long,unsigned int> loc_density;
	map <unsigned long,unsigned int>::iterator ld_it;
	while (cin.good()) { 
		string line;
		stringstream tmp;
		getline(cin,line);
		if (!cin.good()) { break; }
		cout << line << endl;

		if (line[0] == '@') { continue; }

		string trash;
		string chr_string;
		unsigned int chr;
		unsigned int loc;

		tmp << line;
		tmp >> trash;
		tmp >> trash;
		tmp >> chr_string;
		tmp >> loc;

		if (chr_string[0] == '*') { continue; }
		chr = atoi(chr_string.c_str());

		unsigned long unique_site =  chr;
		unique_site = unique_site<<32;
		unique_site = unique_site | loc;
		//cout << chr << " " << loc << " " << unique_site << endl;
		ld_it = loc_density.find(unique_site);
		if (ld_it == loc_density.end()) { 
			loc_density.insert(pair<unsigned long, unsigned int> (unique_site,1));
			ld_it = loc_density.find(unique_site);
		}
		else { 
			ld_it->second++;
		}
	}

	ofstream count_file;
	if (CountFile == NULL) { 
		count_file.open("sites.count");
	}
	else { 
		count_file.open(CountFile);
	}
	if (count_file.is_open()) { 
		for (ld_it = loc_density.begin(); ld_it != loc_density.end(); ld_it++) { 
			unsigned int chr = ld_it->first>>32;
			unsigned int loc = (ld_it->first<<32)>>32;
	
			count_file << chr << " " << loc << " " << ld_it->second << endl;
		}
	}
}



