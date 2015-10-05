//============================================================================
// Name        : K-masker
// Author      : Bryan Downie
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#define MAX_SEQUENCE_LENGTH 100000  // Change this to the maximum sequence length


#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <tr1/array>
#include <sstream>


#include <omp.h>
#include <set>

#include <bitset>
#include <limits.h>

#include <time.h>
#include <math.h>

#include <getopt.h>

#include <boost/unordered_set.hpp>


//#include <boost/filesystem/operations.hpp>
//#include <boost/thread.hpp>
//#include <boost/array.hpp>

#if ID_TYPE==1
#define PROGRAM_NAME "Kmasker_reverse"
#else
#define PROGRAM_NAME "Kmasker"
#endif
#define MAX_BUFFER_SIZE 100000

using namespace std;

int DEBUG = 0;

// Function declarations here
void printUsage();
string rev_comp(string strNucleotides);
void encode_kmer(string strKmer, unsigned long *ulValue);
string print_time(time_t *timer);

class CL_nucleotide {
		char value;
		vector <char> fwd_lookup;
		vector <char> rev_lookup;
	public:
		CL_nucleotide(char);

		void set_value (char c) { value = c; };
		char get_value(void) { return value; };

		int nuc_to_int(char c) { return fwd_lookup[c]; } ;
		char int_to_nuc(int i) { return rev_lookup[i]; } ;

		char complement(void);
};

char CL_nucleotide::complement (void) { 
	if (value == 'N') { return value; }
	else { 
		//cerr << abs(nuc_to_int(value) - 3);
		char c = (char) int_to_nuc(abs(nuc_to_int(value) - 3));
		return c;
	}
}

CL_nucleotide::CL_nucleotide (char c) {
	value = c;
	fwd_lookup.resize(255, 255);
	rev_lookup.resize(4);

	fwd_lookup['A'] = 0;
	fwd_lookup['a'] = 0;
	fwd_lookup['C'] = 1;
	fwd_lookup['c'] = 1;
	fwd_lookup['G'] = 2;
	fwd_lookup['g'] = 2;
	fwd_lookup['T'] = 3;
	fwd_lookup['t'] = 3;
	fwd_lookup['N'] = 9;
	fwd_lookup['n'] = 9;

	rev_lookup[0] = 'A';
	rev_lookup[1] = 'C';
	rev_lookup[2] = 'G';
	rev_lookup[3] = 'T';
	rev_lookup[9] = 'N';
}


class CL_kmer {
		string strKmer;
		vector <char> fwd_lookup;
		vector <char> rev_lookup;
	public: 
		CL_kmer(void);
		string reverse_complement(void);
		unsigned long encode (void);
		void set_kmer(string);
};

string CL_kmer::reverse_complement(void) {
	string return_string;

	for (string::reverse_iterator i = strKmer.rbegin(); i != strKmer.rend(); i++) {
		CL_nucleotide c(*i);
		return_string.push_back(c.complement());
	}

	return return_string;
}

CL_kmer::CL_kmer (void) { 
	//strKmer = input;

	fwd_lookup.resize(255, 255);
	rev_lookup.resize(4);

	fwd_lookup['A'] = 0;
	fwd_lookup['a'] = 0;
	fwd_lookup['C'] = 1;
	fwd_lookup['c'] = 1;
	fwd_lookup['G'] = 2;
	fwd_lookup['g'] = 2;
	fwd_lookup['T'] = 3;
	fwd_lookup['t'] = 3;
	fwd_lookup['N'] = 9;
	fwd_lookup['n'] = 9;

	rev_lookup[0] = 'A';
	rev_lookup[1] = 'C';
	rev_lookup[2] = 'G';
	rev_lookup[3] = 'T';
}

unsigned long CL_kmer::encode (void) { 
	unsigned long power = 1;
	unsigned long ulKmer = 0;
	string revcompKmer;
	//unsigned long ulRevKmer = 0;
	for (string::iterator i = strKmer.begin(); i != strKmer.end(); i++) { 
		ulKmer += fwd_lookup[*i] * power;
		power *= 4;
	}
/*
	for (string::reverse_iterator i = strKmer.rbegin(); i != strKmer.rend(); i++) {
		CL_nucleotide c(*i);
		revcompKmer.push_back(c.complement());
	}

	for (string::iterator i = revcompKmer.begin(); i != revcompKmer.end(); i++) { 
		ulRevKmer += fwd_lookup[*i] * power;
		power *= 4;
	}
//	cout << strKmer << endl << revcompKmer << endl;
*/
/*
	if (ulKmer < ulRevKmer) { return ulKmer; }
	cerr << "Badly formed jellyfish result " << strKmer << endl;
	return ulRevKmer;
	*/
	return ulKmer;
}

void CL_kmer::set_kmer (string kmer_to_set) { 
	strKmer = kmer_to_set;
}

class CL_buffer {
		vector <string> accession;
		vector <string> fasta;
		vector <string> qual;
		unsigned long max_size;
		int good_flag;
		unsigned long last_index;
	public:
		CL_buffer(void);
		void push_acc (string s) { accession.push_back(s); };
		void push_fasta (string s) { fasta.push_back(s); last_index++; };
		void push_qual (string s) { qual.push_back(s); }
		int is_good (void) { return good_flag; } ;
		void set_good (int i) { good_flag = i; };
		int is_full (void);
		void clear(void);
		unsigned long get_last_index(void) { return last_index; };
		string get_fasta(unsigned long l) { return fasta[l]; };
		string get_qual(unsigned long l) { return qual[l]; };
		string get_acc(unsigned long l) { return accession[l]; };
		void set_fasta(unsigned long l, string s) { fasta[l] = s; };
		void set_qual(unsigned long l, string s) { qual[l] = s; };
		void init (void) { good_flag = 1; };
};

void CL_buffer::clear(void) { 
	fasta.clear();
	accession.clear();
	qual.clear();
	last_index = 0;
}

int CL_buffer::is_full(void) {
	if (accession.size() >= max_size) { return 1; }
	return 0;
}

CL_buffer::CL_buffer(void) { 
	last_index = 0;
	good_flag = 1; 
	max_size = MAX_BUFFER_SIZE;
} 

int main(int argc, char* argv[]) {

	//if (argc != 6) { printUsage(); return 0; }
	char c;
	//char *FastaFile = NULL;
	char *KmerFile = NULL;

	unsigned int iLowerBoundry = 0;
	unsigned int iUpperBoundry = 0;
	unsigned int iKmerSize = 0;
	int QualityOffset = 0;

	unsigned int min_length = 1;
	bool progress = 0;
	vector <string> FastaFiles;

	while ((c = getopt (argc, argv, "ghdk:l:u:s:n:")) != -1) {
		switch(c) {
			//case 'f':
			//	FastaFile = optarg;
			//	break;
			case 'k':
				KmerFile = optarg;
				break;
			case 'l':
				unsigned int iLowerBoundryTmp;
				iLowerBoundryTmp = atoi(optarg);
				if (iLowerBoundryTmp == 0) { iLowerBoundryTmp = 1; }
				iLowerBoundry = iLowerBoundryTmp;
				break;
			case 'u':
				iUpperBoundry = atoi(optarg);
				break;
			case 's':
				iKmerSize = atoi(optarg);
				break;
			case 'n':
				min_length = atoi(optarg);
				break;
			case 'g':
				progress = 1;
				break;
			case 'd':
				DEBUG = 1;
				break;
			case 'h':
				printUsage();
				return 0;
			case 'q':
				QualityOffset = atoi(optarg);
				break;
			default:
				printUsage();
				return 0;
		}
	}


	time_t timer;
	if (progress) { 
		cerr << endl<< print_time(&timer)<< "Starting Program .." << endl;
	}

//#if ID_TYPE==1
	if (!(KmerFile && iLowerBoundry)) { 
//#else
//	if (!(KmerFile && iLowerBoundry && iUpperBoundry && iKmerSize)) { 
//#endif
		printUsage();
		return 0;
	}

	if (optind > argc) { 
		cerr << "Provide at least one file to mask!" << endl;
		printUsage();
		return 0;
	}

//cout << optind << " " << argc <<  " " << argv[optind] << endl;
	while (optind < argc) {
		FastaFiles.push_back(argv[optind]);
		optind++;
	}

	for (vector <string>::iterator it = FastaFiles.begin(); it != FastaFiles.end(); it++) { 
		ifstream fhSeq;
		fhSeq.open((*it).c_str());
		if (!fhSeq.is_open()) { 
			cout << "Couldn't open sequence file " << 
					 *it << endl; exit(-1); 
		}
	}

	ifstream fhKmer;
	fhKmer.open(KmerFile);

//	unsigned int iLowerBoundry = atoi(argv[1]);
//	unsigned int iUpperBoundry = atoi(argv[2]);
//	unsigned int iKmerSize = atoi(argv[5]);
	//unsigned int threads = 16;
	
	if (iKmerSize > 32) {
		cerr << "Choose a kmer size below 32" << endl;
		return 1;
	}

	if (!fhKmer.is_open()) { 
		cout << "Couldn't open kmer file " << KmerFile << endl; 
		return 1;
	}
	if (progress) { 	
		cerr << print_time(&timer) << "Using lower threshold: " << iLowerBoundry << endl;
		cerr << print_time(&timer) << "Using upper threshold: " << iUpperBoundry << endl;
		//cerr << print_time(&timer) << "Using kmer size: " << iKmerSize << endl;
		cerr << print_time(&timer) << "Initializing kmer hash.." << endl;
	}

//	std::ios::sync_with_stdio(false);
	//boost::unordered_set<unsigned long> setKmers;
//	boost::unordered_set<unsigned long>::iterator kmer_it;

	if (iKmerSize == 0) { 
		string KmerTmp;
		getline (fhKmer, KmerTmp);
		getline (fhKmer, KmerTmp);
		iKmerSize = KmerTmp.size();
		if (progress) { 
			cerr << endl << print_time(&timer) << "Guessing kmer size of " << iKmerSize << endl;
		}
		fhKmer.seekg(0);
	}


	unsigned long numLines = 0;
	unsigned int lastSize = 1;
	set <long> MaskKmers;
	vector <bool> vbMaskKmers;
	unsigned long power = pow(4,iKmerSize);
	try {
		vbMaskKmers.resize(power,0);
  	}
	catch (exception& e)
 	{
   	 	cout << "Could not allocate memory for k-mer index. Try using a lower k-mer value." << endl;
		return 0;
  	}
//	bool *abMaskKmers;
//	abMaskKmers = new bool[17179869184];

	if (progress) { 
		cerr << print_time(&timer) << "Reading kmer file .." << endl;
		cerr << print_time(&timer)<< "Progress: 0";
	}
	CL_buffer kmer_buffer;

	vector<long> EncodedKmersToInsert;
	EncodedKmersToInsert.resize(MAX_BUFFER_SIZE,-1);

	//vector<bool> MaskKmers;


	while (kmer_buffer.is_good()) {
		string strLine, strCount;
		string strCompactKmer;
		
		getline (fhKmer, strCount);
		if (!fhKmer.good()) { kmer_buffer.set_good(0); }
		else { 
			kmer_buffer.push_acc(strCount);
			getline (fhKmer, strLine);
			//try {
			if (strLine.size() != iKmerSize) {
				string failure = "bad kmer entry found. Check all k-mer sizes and rerun.";
				cout <<  "bad kmer entry found. Check k-mers in kmer file: " << KmerFile << endl;
				return 9;
			}
			//}
			//catch (string s) {
			//	cout << s << endl;
			//	return 9;
			//}
			kmer_buffer.push_fasta(strLine);
			//cout << strCount << endl << strLine << endl;
		}

		if (kmer_buffer.is_full() || !kmer_buffer.is_good()) { 
			#pragma omp parallel
			{

			stringstream tmp;
			unsigned int uiCount;
			CL_kmer kmer;

			
			#pragma omp for schedule (guided,100) 
			for (unsigned long l = 0; l < kmer_buffer.get_last_index(); l++) { 
				string counts = kmer_buffer.get_acc(l);
				string seq = kmer_buffer.get_fasta(l);
				//cout << seq << endl << counts << endl;

				counts.erase(0,1);
				uiCount = atoi(counts.c_str());
				//tmp << counts;
			//	tmp >> uiCount;
			//	tmp.clear();

				//#if ID_TYPE == 1 // reverse kmasker
				//	if (uiCount < iLowerBoundry) { 
				//#else
					if ((uiCount <= iLowerBoundry) || ((iUpperBoundry > 0) && (uiCount >= iUpperBoundry))) {
				//#endif

				// Uncomment this line to "reverse k-mask"
				//if (uiCount < iUpperBoundry) 
						kmer.set_kmer(seq);
						unsigned long ulKmer = kmer.encode();
						vbMaskKmers[ulKmer] = 1;
						//EncodedKmersToInsert[l] = ulKmer;
	

						//vbMaskKmers[ulKmer] = 1;
		//			#pragma omp critical
		//			{
		//			setKmers.insert(ulKmer);
		//			}
					}
				}
			// End parallel pragma
			}
		//	for (unsigned long l = 0; l < kmer_buffer.get_last_index(); l++) { 
		//		long val = EncodedKmersToInsert[l];
		//		if (val != -1) { 
		//			MaskKmers.insert(val);
		//		}
		//	}
			kmer_buffer.clear();
		}


		
		if (progress) { 
			numLines++;
			if ((numLines % MAX_BUFFER_SIZE) == 0) {
				for (unsigned int k = 0; k < lastSize; k++) {
					cerr << "\b";
				}
				stringstream Output;
				string strOutput;
				Output << numLines;
				Output >> strOutput;
				lastSize = strOutput.size();
				cerr << numLines;
			}
		}
	}
	//buffer.init();
//for (kmer_it = setKmers.begin(); kmer_it != setKmers.end(); kmer_it++) { 
//	cout << *kmer_it << endl;
//}


//	cerr << endl;

	//cerr << print_time (&timer) << "Masking sequences.." << endl;




	//int *aiConvertedString = new int[MAX_SEQUENCE_LENGTH];	
	//int *aiReverseConvertedString = new int[MAX_SEQUENCE_LENGTH];	

//	bool *baMask;
//	baMask = new bool[MAX_SEQUENCE_LENGTH];
//	int *aiConvertKmer = new int[MAX_SEQUENCE_LENGTH];	
//	for (int i = 0; i < MAX_SEQUENCE_LENGTH; i++) {
//		aiConvertKmer[i] = 0;
//	}


	for (vector <string>::iterator it = FastaFiles.begin(); it != FastaFiles.end(); it++) { 
		ifstream fhSeq;
		bool is_fastq = 0;

		ofstream fhOut;
		string output_string;

		fhSeq.open((*it).c_str());
		if (!fhSeq.is_open()) { 
			cerr << "Couldn't open sequence file " << *it << endl;
			continue;
		}

		stringstream tmp;
		//string string_tmp;
		output_string.assign(*it);
		output_string.append(".k");

		if (iUpperBoundry == 0) { 
			tmp << iLowerBoundry << "+";
		}
		else { 
			tmp << iLowerBoundry << "-" << iUpperBoundry;
		}	
		output_string.append(tmp.str());
		//output_string.append("-");

		//tmp << iUpperBoundry;
		//output_string.append(tmp.str());

		fhOut.open(output_string.c_str());

		if (!fhOut.is_open()) { 
			cerr << "Couldn't open output file " << output_string << endl;
			continue;
		}

		CL_buffer buffer;
		buffer.clear();
		numLines = 0;
		lastSize = 1;
//		CL_buffer buffer;
		string check_fastq;
		getline (fhSeq,check_fastq);
		char MaskQualChar;
		if (check_fastq[0] == '@') { 
			is_fastq = 1;
			if (!QualityOffset) { 
				if (progress) { 
					cerr << endl << print_time(&timer) << "Guessing quality offset";
				}
				QualityOffset = 64;
				getline (fhSeq,check_fastq);
				getline (fhSeq,check_fastq);
				getline (fhSeq,check_fastq);
				for (unsigned int str_it = 0; str_it < check_fastq.size(); str_it++) {
					int phred = check_fastq[str_it] - 64;
			//		cout << check_fastq[str_it] << " " << phred << endl;
					if (phred < 0) { 
						QualityOffset = 33; 
						MaskQualChar= '!';
						str_it = check_fastq.size();

					}
				}
				if (progress) { 
					cerr << endl << print_time(&timer) << "Using quality offset of " << QualityOffset;
				}
			}
			if (QualityOffset == 64) { 
				MaskQualChar = 'B';
			}
			else if (QualityOffset == 33) { 
				MaskQualChar = '!';
			}
			else { 
				cerr << endl << print_time(&timer) << "Unknown quality offset used!" << endl;
				throw "Unknown quality offset";
				//return -3;
			}

			if (progress) { 
				cerr << endl << print_time (&timer) << "Writing new fastq sequences from file " << *it << ". Progress: 0";
			}
		}
		else if (check_fastq[0] == '>') { 
			if (progress) { 
				cerr << endl << print_time (&timer) << "Writing new fasta sequences from file " << *it << ". Progress: 0";
			}
		}
		else {
			cerr << endl << "Unknown file format in " << *it << endl;
			throw "Unknown file format";
		}
		// Return to beginning of file.
		fhSeq.seekg(0);

		while (buffer.is_good()) { 
			string strLine, strAcc;
			string strKmerFwd, strKmerRev;
			
	
			getline (fhSeq, strLine);
			if (!fhSeq.good()) { buffer.set_good(0); }
			else { 
				buffer.push_acc(strLine);
				getline (fhSeq, strLine);
				buffer.push_fasta(strLine);
				if (is_fastq) { 
					getline (fhSeq, strLine);
					getline (fhSeq, strLine);
					buffer.push_qual(strLine);
				}
			}
	
			if (buffer.is_full() || !buffer.is_good()) { 
				if (DEBUG) { 
					omp_set_num_threads(1);
				}
				#pragma omp parallel 
				{
				vector <int> baMask;
				vector <int> aiConvertKmer;
				baMask.resize(MAX_SEQUENCE_LENGTH);
				aiConvertKmer.resize(MAX_SEQUENCE_LENGTH);
				#pragma omp for schedule(dynamic,100)
				for (unsigned long l = 0; l < buffer.get_last_index(); l++) { 
					if (DEBUG) { 
						cout << "---------- " << l << " --------------" << endl;
					}
					vector <int> viConvertedString;
					vector <int> viReverseConvertedString;
					vector <int> tmpviReverseConvertedString;
					vector <int> viConvertKmer;
	
	
					string sequence = buffer.get_fasta(l);
					string quality;
					if (is_fastq) { 
						quality = buffer.get_qual(l);
					}
	
					unsigned int iStringSize = sequence.size();
					if ((iStringSize < min_length) || (iStringSize < iKmerSize) || (sequence.size() != quality.size())) {
						sequence.erase();
						sequence.resize(5,'N');
						buffer.set_fasta(l,sequence);
						if (is_fastq) { 
							buffer.set_qual(l,"BBBBB");
						}
						continue;
					}
	
			
					unsigned int iMaskVectorSize = 1+  iStringSize - iKmerSize;
					unsigned int i, h;
	
	
	
					viConvertedString.resize(iStringSize);
					viReverseConvertedString.resize(iStringSize);
					viConvertedString.clear();
					viReverseConvertedString.clear();
			
		
					for (i = 0, h = iStringSize - 1; i < iStringSize; i++, h--) {
						CL_nucleotide nucleotide(sequence[i]);
						viConvertedString[i] = nucleotide.nuc_to_int(sequence[i]);

						if (viConvertedString[i] > 9) { 
							cerr << "Bad nucleotide found in sequence!" << endl;
							cerr << "Entry " << l << " in buffer" << endl;
							cerr << sequence << endl;
							if (is_fastq) { 
								cerr << quality << endl;
							}
							throw "Error";
						}
	
						viReverseConvertedString[h] = nucleotide.nuc_to_int(nucleotide.complement());
						// Don't convert kmers that contain Ns.
						//cout << sequence[i];
						if (nucleotide.get_value() == 'N') { 
							int j = i - iKmerSize + 1;
							if (j < 0) { j = 0; }
							for (;j < i + 1; j++) {
								aiConvertKmer[j] = 0;
							}
						}
						else { aiConvertKmer[i] = 1; }
					}
		
					baMask[iMaskVectorSize] = 1;
					baMask[iMaskVectorSize+1] = 1;
					for (i = 0; i <  iMaskVectorSize; i++) {
						baMask[i] = 0;
						if (aiConvertKmer[i]) { 
							//CL_nucleotide nucleotide('N');
							unsigned long ulFwd, ulRev, ulKmer;
							ulFwd = 0;
			
							ulRev = 0;
							unsigned long power = 1;
							//string strKmer= sequence.substr(i,iKmerSize);
			
							unsigned int j, k;
							for (j = i, k = iMaskVectorSize - i - 1; j < i + iKmerSize; j++, k++) {
								ulFwd += viConvertedString[j] * power;
								ulRev += viReverseConvertedString[k] * power;
								power *= 4;
							}
	
			
							if (ulFwd < ulRev) { ulKmer = ulFwd; }
							else { ulKmer = ulRev;}
							
							if (vbMaskKmers[ulKmer]) { 
								baMask[i] = 1;
							}
			
				
						}
					}	
	

					if (DEBUG) { 
						cout << "1Seq: " << sequence << endl << "Mask: "; 
						for (unsigned int k = 0; k < iMaskVectorSize; k++) {
							bool bKmer = baMask[k];
							cout << bKmer; 
						}
						cout << endl;
					}

					vector <unsigned int> start_coord;
					vector <unsigned int> end_coord;
					bool last_maskKmer = baMask[0];
					// Check the beginning
					if (baMask[0]) { 
						start_coord.push_back(0);
					}
					for (unsigned int k = 1; k < iMaskVectorSize; k++) {
						bool maskKmer = baMask[k];
						// Start of a sequence
						if (maskKmer && !last_maskKmer) {
							if (k <= iKmerSize) {
								for (unsigned int m = k; m < iKmerSize; m++) {
									baMask[m] = 1;
								}
								start_coord.assign(1,0);
							}
							else if ((end_coord.size() > 0) && ((k - end_coord.back()) < iKmerSize)) {
								for (unsigned int m = end_coord.back(); m < k; m++) { 
									baMask[m] = 1;
								}
								end_coord.pop_back();
							}
							else { 
								unsigned int index = k + (iKmerSize - 1);
								if (index < iMaskVectorSize) { 
									while (!baMask[index] && (index > k)) {
										index--; 
									}
									if (index > k) { 
										for (unsigned int m = k; m < index; m++) {
											baMask[m] = 1;
										}
										start_coord.push_back(k + iKmerSize - 1);
										baMask[k] = 1;
										maskKmer = 1;
									}
									else { 
										maskKmer = 0;
										baMask[k] = 0;
									}
								}
								else { 
									start_coord.push_back(k);
								}
							}
						}
						// End of a sequence
						else if (!maskKmer && last_maskKmer && (end_coord.size() < start_coord.size())) {
							unsigned int index = k + (iKmerSize - 1);
							if (index < iMaskVectorSize) {
								while (!baMask[index] && (index >= k)) {
									index--;
								}
								if (index > k) { 
									for (unsigned int m = k; m < index; m++) { 
										baMask[m] = 1;
									}
									maskKmer = 1;
								}
								end_coord.push_back(index);
							}
							else {
								end_coord.push_back(sequence.size());
							//	baMask[k] = last_maskKmer;
							//	maskKmer = last_maskKmer;
							}
						}
						last_maskKmer = maskKmer;
					}
					if (last_maskKmer) { end_coord.push_back(sequence.size()); }

					if (DEBUG) { 
						cout << "PMsk: "; 
						for (unsigned int k = 0; k < iMaskVectorSize; k++) {
							bool bKmer = baMask[k];
							cout << bKmer; 
						}
						cout << endl;
					}

					if ((end_coord.size() > 0) && ((sequence.size() - end_coord[0]) <= iKmerSize)) { 
						for (unsigned int i = sequence.size() - 1; i > end_coord[0]; i--) {
							baMask[i] = 1;
						}
						end_coord[0] = sequence.size() -1;
					}

					for (unsigned int i = 0; i < start_coord.size(); i++) {
						unsigned int start_mask = start_coord[i];
						unsigned int end_mask = end_coord[i];
						//cout << start << " " << end << endl;

						//int length = (end - start) - iKmerSize;
						//if ((start == 0) || (end >= (iMaskVectorSize - 1))) {
							//unsigned int start_mask = start;
							//unsigned int end_mask;
							//if (start == 0) { start_mask = 0; } 
							//else {  
							//	start_mask = iKmerSize + start; 
						//	}
						//	if (end >= (iMaskVectorSize - 1)) { end_mask = sequence.size() - 1; }
						//	else { 
						//		end_mask = end; 
						//	} 
							if ((!baMask[end_mask]) && (baMask[end_mask -1])) { 
								end_mask--;
							}
							for (unsigned int j = start_mask; j <= end_mask; j++) {
								sequence[j] = 'N';
								if (is_fastq) { 
									quality[j] = MaskQualChar;
								}
							}
						//}
					}
	
					if (DEBUG) { cout << "2Seq: " << sequence << endl; }
					for (unsigned int i = 0; i < sequence.size(); i++) {
						if (sequence[i] != 'N') { 
							sequence.erase(0,i);
							if (is_fastq) { 
								if (quality.size() <= i) {
									//cerr << "Bad quality length for " << l << " with quality " << quality << endl;
									sequence.erase();
									quality.erase();
								}
								else { 
									quality.erase(0,i);
								}
							}
							i = sequence.size();
						}
					}
					if (sequence.size() > 0) { 
						for (int i = sequence.size() - 1; i >= 0; i--) { 
							if (sequence[i] != 'N') { 
								sequence.resize(i);
								if (is_fastq) { 
									quality.resize(i, MaskQualChar);
								//	if (quality.size() <= i+1) {
										//cerr << "Bad quality length for " << l << " with quality " << quality << endl;
								//		quality.resize(i+1,MaskQualChar);
								//	}
								//	quality.erase(i+1,count);
								}
								i = -1;
							}
							else if (i <= 0) { 
								sequence.erase();
								if (is_fastq) { 
									quality.erase();
								}
							}
						}
					}
	
					if (sequence.size() < min_length) { 
						sequence.erase();
						sequence.resize(5,'N');
						if (is_fastq) { 
							quality.erase();
							for (unsigned int tmp = 1; tmp <= 5; tmp++){ 
								quality.push_back(MaskQualChar);
							}
						}
					}
					if (DEBUG) { cout << "3Seq: " << sequence << endl; }
	
					buffer.set_fasta(l,sequence);
					buffer.set_qual(l,quality);
				}
				// end of pragma brackets
				}
				for (unsigned long l = 0; l < buffer.get_last_index(); l++) { 
					fhOut << buffer.get_acc(l) << endl << buffer.get_fasta(l) << endl;
					if (is_fastq) { 
						fhOut << "+" << endl << buffer.get_qual(l) << endl;
					}
				}
				buffer.clear();
			}
			if (progress) { 
				numLines++;
				if ((numLines % MAX_BUFFER_SIZE) == 0) {
					for (unsigned int k = 0; k < lastSize; k++) {
						cerr << "\b";
					}
					stringstream Output;
					string strOutput;
					Output << numLines;
					Output >> strOutput;
					lastSize = strOutput.size();
					cerr << numLines;
				}
			}
		}
		
		fhOut.close();
		if (progress) { cerr << endl << print_time(&timer) <<"Finished masking sequences from file " << *it<< endl; }
	}


	return 0;
}


void printUsage() {
#if ID_TYPE==1
	cout << "Usage:\t" <<  PROGRAM_NAME  << " -l <lower bound> -k <kmer file> -s <kmer size> [-n <minimum_sequence_length>] <seq file 1> <seq file 2> ...\n\n" 
#else
	cout << "Usage:\t" <<  PROGRAM_NAME  << " -l <lower bound> -u <upper bound> -k <kmer file> -s <kmer size> [-n <minimum_sequence_length>] <seq file 1> <seq file 2> ...\n\n" 
#endif
		 << "\tMasks all kmers which occur outside of frequency bounds (bounds must be > 0).\n\n"
	     << "\tInput:\n"
	     << "\t\t<kmer file>: file of kmers to be sequenced\n"
	     << "\t\t<seq file> : file with sequences to be masked.\n\n"
	     << "\tOutput:\n"
	     << "\t\tMulti-fasta containing masked sequences.\n"
	     << endl << endl
		 << "\tVersion 1.0b" << endl
		 << "\tCopyright Bryan Downie (bdownie@fli-leibniz.de)" << endl << endl;
}

string print_time(time_t *timer) { 
	string readableTime;
	size_t newline;

	time(timer);
	readableTime = ctime(timer);
	newline = readableTime.find('\n');
	readableTime.erase(newline);
	
	string return_time = "[ ";
	return_time.append(readableTime);
	return_time.append(" ] ");
	return return_time;
}
