//============================================================================
// Name        : Scaffolder.cpp
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
#include <omp.h>

//#include <google/dense_hash_map>
//#include <boost/unordered_map.hpp>

#include <boost/algorithm/string.hpp>

#define MAX_BUFFER_SIZE 100000
#define TRUE 1
using namespace std;

// Function declarations
//string print_time();
void printUsage();
string get_working_directory (unsigned int LastDir);
//string print_time();

unsigned int parse_graph(char *, vector<string> *);
unsigned int parse_scaffold_file(char *, unsigned int *, vector<unsigned int> *, vector<set <unsigned int> > *, vector<string> *);
unsigned int parse_seq(string , vector<string> *);
vector<string> parse_sam (string , string , vector<unsigned int> *, unsigned int , string , bool );
unsigned int parse_contig(char *, unsigned int *, vector<unsigned int> *, vector<set <unsigned int> > *, vector<string> *);
//typedef boost::unordered_map<string, string> hashmap;
string rev_comp(string strNucleotides);
int parse_graph_for_scaffolds (vector <string> *, map <string, unsigned int> *, map <string,int> *, map <string, unsigned int> *, map <string, unsigned int> *);
string build_scaffold (map <string, unsigned int> *, map <string,int> *, map <string, unsigned int> *, map <string, unsigned int> *, string , vector <string> *);
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
	char *ContigFile = NULL;
	//char *SeqFile = NULL;
	//char *SamFile = NULL;
	char *GraphFileLoc = NULL;
	char *LibConfig = NULL;

	// All purpose index variables
	//unsigned int i, j; 
	unsigned int scaffold_number;
	int threads = 1;
	//unsigned long numLines = 0;
	//int lastSize = 0;

	// Get the different options necessary to run
	while ((c = getopt (argc, argv, "m:t:i:c:r:g:fu")) != -1) {
		switch(c) {
			// Scaffold file prepared by scaffolder_int
			case 'i':
				ScaffoldFile = optarg;
				break;
			// Number of threads to use (OpenMP)
			case 't':
				threads = atoi(optarg);
				break;
			// Maximum number of threads to use when processing multiple libraries simultaneously (not recommended)
			case 'm':
				sub_max_threads = atoi(optarg);
				break;
			// Library file (see wrapper.pl for documentation)
			case 'r':
				LibConfig = optarg;
				break;
			// Start assembly file - requires integer accessions (i.e. no Contig_1)
			case 'c':
				ContigFile = optarg;
				break;
			// Scaffold graph file (obsolete)
			case 'g':
				GraphFileLoc = optarg;
				break;
			// Obsolete
			case 'f':
				fill_gap_only = 1;
				break;
			// Obsolete
			case 'u':
				output_unmatched = 1;
				break;
			default:
				printUsage();
				return 0;
		}
	}

	// If a req'd parameter isn't given
	if (NULL == ContigFile) { 
		printUsage();
		return 0;
	}
	// If no library, perform scaffolding only
	if (NULL == LibConfig) { no_lib = 1; }

	// Start the timer
	cerr << endl<< "Starting Program .." << endl;

	// This section gives the number of the last contig accession. Should be replaced to not require system call, 
	// but it low priority because it happens only once.
	stringstream GetReadNum;
	string GetReadNumCmd, Line;
	ifstream File;
	GetReadNumCmd =  "tail -n 2 ";
	GetReadNumCmd.append(ContigFile);
	GetReadNumCmd.append(" | head -n 1 | sed 's/>//' > .contig.num\n");
	system(GetReadNumCmd.c_str());
	File.open(".contig.num");
	getline(File, Line);
	File.close();

	// # of contigs in initial assembly
	unsigned long MaxContigs = atol(Line.c_str());
	MaxContigs++;

	// Storage constructs
	// Which contigs belong to which scaffolds and vice versa
	vector<unsigned int> ContigsToScaffold;
	vector<set <unsigned int> > ScaffoldsToContigs;

	// Sequences of reads. Index of reads belongs to the scaffold (contig) to which it aligns in sam file.
	// ReadSeqs will store the unmasked sequence (found in libraries file).
	vector<string> ReadSeqs;

	// Sequence of contigs/scaffolds, indexed by accession number
	vector<string> ContigSeqs;

	// Obsolete (I think)
	vector<vector<string> > ScaffoldsToSAM;

	// For storing strings of scaffold and graph file. Graph file is only used for scaffolding-only mode
	// (i.e. not gap filling)
	vector<string> GraphStrings;
	vector<string> ScaffStrings;

	// Because the scaffold/contig accessions start at 1, placeholder is put in to make all push_backs
	// create a direct relationship between index and accession.
	set<unsigned int> PlaceHolder;
	ScaffoldsToContigs.push_back(PlaceHolder);

	// Output that should be silencable
	cerr << "Initializing...\n";
	cerr << "Contigs detected: " << MaxContigs << endl;
	cerr << "Allocating memory. If reads or contigs don't match detected, "
		 << "check last lines of contig and sam file." << endl;

	// Allocate memory for vectors
	ContigsToScaffold.resize(MaxContigs,0);
	ScaffoldsToSAM.resize(MaxContigs);
	GraphStrings.resize(MaxContigs);
	ScaffStrings.resize(MaxContigs);
	ContigSeqs.resize(MaxContigs);

	unsigned int return_value = 0;
	// Get listing of contig distance and relative orientation for use during scaffolding-only mode.
	if (!(NULL == GraphFileLoc)) { 
		return_value = parse_graph(GraphFileLoc, &GraphStrings);
	}
	map <string, unsigned int> Counts;
	map <string, int> Distances;
	map <string, unsigned int> FirstDirection;
	map <string, unsigned int> SecondDirection;

	// Come back to this later to comment subroutine purpose.
	parse_graph_for_scaffolds(&GraphStrings, &Counts, &Distances, &FirstDirection, &SecondDirection);

	// If in gap-filling mode...
	if (NULL==ScaffoldFile)  { 
		// Store scaffold sequences in ContigSeqs, do some other stuff 
		return_value = parse_contig(ContigFile, &scaffold_number, &ContigsToScaffold, &ScaffoldsToContigs, &ContigSeqs);
	}
	// Otherwise in scaffold only mode
	else {
		return_value = parse_scaffold_file(ScaffoldFile, &scaffold_number, &ContigsToScaffold, &ScaffoldsToContigs, &ScaffStrings);
		// Store contigs in ContigSeqs, do some other stuff
		return_value = parse_contig(ContigFile, &scaffold_number, &ContigsToScaffold, &ScaffoldsToContigs, &ContigSeqs);
		// Not sure about this - need to make sure there's no redundancy. Could be leftover from previous version/workflow
		if (!(NULL == GraphFileLoc) && no_lib) { 
			// For each scaffold, orient the contig in proper order and orientation based on 
			// the graph file and the scaffold file.
			//
			// Also perform overlap detection for gaps which are below a certain size (maybe make
			// it configurable at the command line in future versions if there is interest).
			for (unsigned int scaffold = 1; scaffold < ScaffoldsToContigs.size(); scaffold++) { 
				// Contig indexes which belong to a particular scaffold.
				set <unsigned int> contigs = ScaffoldsToContigs[scaffold];
				// If no contigs belong to a scaffold (data format error, most likely)
				if (contigs.size() == 0) { 
					cerr << "Empty scaffold found! Scaffold: " << scaffold << endl;
					return 0;
				}
				// A scaffold that contains only one contig doesn't need special processing.
				// Output it with "Contig_X" instead of "Scaffold_X" to indicate no scaffolding occurred.
				// Can be changed in future versions if there is interest.
				else if (contigs.size() == 1) { 
					set <unsigned int>::iterator it = contigs.begin();
					unsigned int contig_num = *it;
					cout << ">Contig_" << contig_num << endl << ContigSeqs[contig_num] << endl;
					
				 }
				 else {
					 // Need to properly implement exception handling here.
					try { 
						// Build the new scaffold string
						string output = build_scaffold(&Counts, &Distances, &FirstDirection, &SecondDirection, ScaffStrings[scaffold], &ContigSeqs);
						cout << output;
					}
					 catch (const char *s) {
						 cerr << "Failure to build scaffold: " << scaffold << endl;
						 cerr << "Error message: " << s << endl;
						 return 3;
					 }
					 catch (out_of_range& oor) {
						 cerr << "Out of Range error: " << oor.what() << endl;
					 }

					//cout << output;
				 }
			}
						 
			return 0;
			

			//return_value = parse_graph(GraphFileLoc, &GraphStrings);
		}
	}
	
	// Should be quietable
	cerr << "Number of scaffolds: " << scaffold_number << endl;
	cerr << endl;

	// Libraries for later use
	vector< vector<string> > libraries;

	// Read Library config file (put into a function and put later in main {}
	if (!no_lib) { 
		ifstream fhLibConfig(LibConfig);
		if (!fhLibConfig.is_open()) { 
			cerr << "Couldn't open library config file " << LibConfig << endl; return 0;
		}
	
		// Get relevant info (unmasked read file, masked read file, sam file, insert size, stddev)
		// for each lib
		while (fhLibConfig.good()) {
			string Library;
			string SAM_file;
			string insert_size;
			string stddev;
			string LibraryMasked;
	
			string Line;
			stringstream tmp;
	
			// Each line of the library is whitespace seperated.
			getline (fhLibConfig, Line);
			if (!fhLibConfig.good()) { break; }
			tmp << Line;
			tmp >> Library;
			tmp >> LibraryMasked;
			tmp >> SAM_file;
			tmp >> insert_size;
			tmp >> stddev;
	
			vector <string> output;
			output.push_back(Library);
			output.push_back(LibraryMasked);
			output.push_back(SAM_file);
			output.push_back(insert_size);
			output.push_back(stddev);
			libraries.push_back(output);
			output.clear();
		}

		cerr << "Processing reads for each library" << endl;
		// If processing libraries simulataneously.
		// Haven't found multiple libraries simultaneously makes sense, should probably just remove
		// this functionality.
		if (sub_max_threads == 0) { 
			sub_max_threads = (int)(threads/libraries.size());
		}
	}
	
	// Create hierarchical modulo directory structure based on scaffold number
	cerr << "Creating directory structure.."<< endl;
	create_working_directories(scaffold_number);
	cerr << "Finished creating directory structure.."<< endl;

	// Can run multipl multi-threaded processes at the same time, should be removed in later version.
	omp_set_nested(TRUE);

	// Need to look up again what this does.
	omp_set_dynamic(TRUE);
	// Number of libraries to run simultaneously (usually 1)
	#pragma omp parallel for num_threads(threads) 
	for (unsigned int i = 0; i <= libraries.size(); i++) {
		if (i == libraries.size()) { 
			// The last thread outputs the contig sequence to appropriate directory
			cerr << endl << "Outputting contig data" << endl;
			for (unsigned int Scaffold = 1; Scaffold < scaffold_number; Scaffold++) {
				string dir = get_working_directory(Scaffold);
				string contigs,  match, graph, scaff;
		
				contigs = dir;
				match = dir;
				scaff = dir;
		
				contigs.append("/contigs.fasta");
				fstream ContigFasta;

	
				// Obsolete (probably)
				fstream Graph;
				if (GraphFileLoc != NULL ) { 
					graph = dir;
					graph.append("/reads.graph");
					Graph.open(graph.c_str(), fstream::out | fstream::app);
				}
		
				// Output to contigs.fasta	
				ContigFasta.open(contigs.c_str(), fstream::out | fstream::app);
				for (set <unsigned int>::iterator sit = ScaffoldsToContigs[Scaffold].begin(); sit != ScaffoldsToContigs[Scaffold].end(); sit++) {
					ContigFasta << ">" << (*sit) << endl << ContigSeqs[*sit]<< endl;
					if (GraphFileLoc != NULL ) { 
						Graph << GraphStrings[*sit];
					}
				}
				// Obsolete (probably)
				if (GraphFileLoc != NULL ) { 
					Graph.flush();
					Graph.close();
				}
				ContigFasta.close();
		
				// Obsolete (probably)
				scaff.append("/reads.scaff");
				fstream Scaff;
				Scaff.open(scaff.c_str(), fstream::out | fstream::app);
				Scaff << ScaffStrings[Scaffold] << endl;
				Scaff.close();


				// Obsolete-ish (contigs.fasta and contigs.scaff.fasta are redundant, but 
				// both are required for LocalAssembly scripts)
				string scaffold_file = dir;
				scaffold_file.append("/contigs.scaff.fasta");
				fstream ScaffoldFileOutput;
				ScaffoldFileOutput.open(scaffold_file.c_str(), fstream::out);

				// Obsolete (maybe)
				set <unsigned int> ContigSet = ScaffoldsToContigs[Scaffold];
				if (ContigSet.size() == 1) { 
					set <unsigned int>::iterator it = ContigSet.begin();
					unsigned int contig_num = *it;
					ScaffoldFileOutput << ">Contig_" << contig_num << endl << ContigSeqs[contig_num] << endl;
					
				 }
				 else {
					try { 
						string output = build_scaffold(&Counts, &Distances, &FirstDirection, &SecondDirection, ScaffStrings[Scaffold], &ContigSeqs);
						ScaffoldFileOutput << output;
					}
					catch (const char *s) {
						 cerr << "Failure to build scaffold: " << scaffold_file << endl;
						 cerr << "Error message: " << s << endl;
					}
				}
				ScaffoldFileOutput.close();
				// End obsolete (maybe)
			}
		}
		else { 
			vector<string> library = libraries[i];
	
			// Declarations can be combined with assignments
			string Library;
			string LibraryMasked;
			string SAM_file;
			string insert_size;
			string stddev;
			Library = library[0];
			LibraryMasked = library[1];
			SAM_file = library[2];
			insert_size = library[3];
			stddev = library[4];
			
			// Reads and masked reads (MR are obsolete) which anchor to each scaffold
			vector<string> ScaffoldsToReads;
			vector<string> ScaffoldsToMaskedReads;
	
			// If a particular library is unpaired fragments (insert size of 0)
			if (atoi(insert_size.c_str()) == 0) { 
				try {
					// Assign individual reads to individual scaffolds in non-paired fashion (the 1 flag in last parameter indicates "fragment" mode)
					ScaffoldsToReads = parse_sam (SAM_file, Library, &ContigsToScaffold, scaffold_number, insert_size,1);
					ScaffoldsToMaskedReads = parse_sam (SAM_file, LibraryMasked, &ContigsToScaffold, scaffold_number, insert_size,1);
				}
				catch (char *c) { 
					cerr << "Failed to parse sam";
				}
			}
			else { 
				try {
					// Assign individual reads to individual scaffolds paired fashion (the 0 flag in last parameter indicates "default" mode)
					ScaffoldsToReads = parse_sam (SAM_file, Library, &ContigsToScaffold, scaffold_number, insert_size,0);
					ScaffoldsToMaskedReads = parse_sam (SAM_file, LibraryMasked, &ContigsToScaffold, scaffold_number, insert_size,0);
				}
				catch (char *c) { 
					cerr << "Failed to parse sam";
				}
			}
			// Silencable and not really important. Will get removed when removing multi-library parallelization
			#pragma omp critical 
			{
				cerr << "Finished reading " << insert_size << "bp library" << endl;
			}

			for (unsigned int Scaffold = 1; Scaffold < scaffold_number; Scaffold++) {
				string dir;
				dir = get_working_directory(Scaffold);
				if (ScaffoldsToReads[Scaffold].empty()) { continue; }
				string reads = dir;

				// "reads" string contains the directory path to a particular anchored, unmasked
				// read file. File has the format "reads.<insert size>.<insert size SD>.fast<a/q>"
				reads.append("/reads.");
				reads.append(insert_size.c_str());
				reads.append(".");
				reads.append(stddev.c_str());
				if (ScaffoldsToReads[Scaffold][0] == '@') { 
					reads.append(".fastq");
				}
				else if (ScaffoldsToReads[Scaffold][0] == '>') { 
					reads.append(".fasta");
				}
				else { 
					cerr << "Unknown read type found!" << endl;
				}
				fstream ReadFasta;
				ReadFasta.open(reads.c_str(), fstream::out | fstream::app);
				ReadFasta << ScaffoldsToReads[Scaffold];
				ReadFasta.flush();
				ReadFasta.close();

				// This section is for outputting the masked reads to the directory and is obsolete
				/*
				string maskedreads = dir;
				maskedreads.append("/reads.");
				maskedreads.append(insert_size.c_str());
				maskedreads.append(".");
				maskedreads.append(stddev.c_str());
				maskedreads.append(".masked.fasta");
				fstream ReadMaskedFasta;
				ReadMaskedFasta.open(maskedreads.c_str(), fstream::out | fstream::app);
				ReadMaskedFasta << ScaffoldsToMaskedReads[Scaffold];
				ReadMaskedFasta.flush();
				ReadFasta.close();
				*/

			}
		}
	}
	
	// All done!
	cerr << "Finished!" << endl;

	return 0;
}

vector<string> parse_sam (string SAM, string Reads, vector<unsigned int> *Contigs2Scaffold, unsigned int scaffold_num, string insert_size, bool single_read) { 

	// Use single_read = 1 to parse single reads only.
	ifstream fhSAM(SAM.c_str());
	if (!fhSAM.is_open()) { 
		cerr << "Couldn't open SAM file " << SAM << endl; exit(-1); 
	}
	
	cerr << Reads << endl;
	// Buffers for multi-threaded processing
	vector <string> Buffer;
	vector <string> SeqBuffer;

	// Declaration for vector that will be passed back. Probably more efficient to declare in calling function 
	// and pass as pointer.
	vector <string> Scaffolds2Reads;
	Scaffolds2Reads.resize(scaffold_num);
	unsigned int Lines = 0;
	bool good_SAM = 1;
	unsigned long numLines = 0;
	unsigned long lastSize = 1;

	ifstream fhSeq(Reads.c_str());
	if (!fhSeq.is_open()) {  
		throw "Couldn't open sequence file ";
	}
	bool is_fastq = 0;
	string test_fastq;
	// Check to see if the first entry starts with "@" character. If so, it's a FASTQ file.
	// Otherwise it's FASTA
	getline(fhSeq, test_fastq);
	if (test_fastq[0] == '@') { 
		is_fastq = 1;
	}
	fhSeq.seekg(0);


	// This output shows up twice, so the second copy should be eliminated.
	cerr << "Processing " << SAM << "  ";
	cerr <<endl;
	cerr << "Using " << sub_max_threads << " threads per subthread..." << endl;
	// Obsolete or unimplemented function (output unmatched)
	fstream Unmatched;
	if (output_unmatched) { 
		string UnmatchedReadsFile;
		UnmatchedReadsFile = Reads;
		UnmatchedReadsFile.append(".unmatched");
		Unmatched.open(UnmatchedReadsFile.c_str(), fstream::out );
	}


	vector<string> Scaffolds2ReadsSeq;
	vector<unsigned int> Scaffolds2ReadsScaffold;

	Scaffolds2ReadsSeq.resize(MAX_BUFFER_SIZE);
	Scaffolds2ReadsScaffold.resize(MAX_BUFFER_SIZE);

	vector<string> UnmatchedBuffer;
	vector<bool> UnmatchedBufferBool;
	UnmatchedBuffer.resize(MAX_BUFFER_SIZE);
	UnmatchedBufferBool.resize(MAX_BUFFER_SIZE, 0);

	// While the SAM file is not EOF
	while (good_SAM) { 
		string Line;

		getline (fhSAM, Line);
		if (Line[0] == '@') { continue; }
		if (!fhSAM.good()) { good_SAM = 0; }

		Buffer.push_back(Line);
		Lines++;
		numLines++;

		// Some kind of progress bar, but seems not implemented. Should probably be removed.
		if ((numLines % 100000) == 0) {
			for (unsigned long k = 0; k < lastSize; k++) {
			}
			stringstream Output;
			string strOutput;
			Output << numLines;
			Output >> strOutput;
			lastSize = strOutput.size();
		}

		// If the buffer is full or the sam file is EOF.
		if ((Lines >= MAX_BUFFER_SIZE) || (good_SAM == 0)) { 
			if (good_SAM == 0) { Buffer.pop_back(); }

			unsigned int SeqLines = 0;
			// Populate sequence buffer
			// Here we assume that the SAM file (excluding header) has the same # of lines as the #
			// of accessions in the Sequence file. Whole system will break otherwise, but error checking every
			// line would probably cause a large performance problem. Maybe check every 100k or so.
			while (fhSeq.good() && (SeqLines < MAX_BUFFER_SIZE )) {
				string Line1, Line2, Line3, Line4, WorkingLine;
				getline (fhSeq, Line1);
				if (!fhSeq.good()) { break; }
				// Make sure that the number of lines for each FASTA/Q entry is correct.
				if ((Line1[0] != '>') && (Line1[0] != '@')) { cerr << Line1 << endl << Line2 << endl; throw "Error"; }
				getline (fhSeq, Line2);
				// Append the insert size to  each sequence accession number
				Line1.append("_");
				Line1.append(insert_size);
				Line1.append("\n");
				Line2.append("\n");

				Line1.append(Line2);
				if (is_fastq) { 
					getline (fhSeq, Line2);
					Line2.append("\n");
					Line1.append(Line2);
					getline (fhSeq, Line2);
					Line2.append("\n");
					Line1.append(Line2);
				}

				SeqBuffer.push_back(Line1);	
				SeqLines++;
			}

			// Fragments only
			if (single_read) { 
				// Performance tweaks can be made in the schedule clause, maybe
				#pragma omp parallel for schedule(guided, 1000) num_threads(sub_max_threads)
				for (unsigned int i = 0; i < Buffer.size(); i++) { 
					string Line1, Cont1, trash;
					unsigned int Contig1, Scaffold;
					unsigned long ReadName1;
					stringstream tmp;
	
					Line1 = Buffer[i];
		
					tmp << Line1;
					tmp >> ReadName1;
					tmp >> trash;
					tmp >> Cont1;
	
					if (Cont1[0] != '*')  { 
						Contig1 = atoi(Cont1.c_str());
						Scaffold = (*Contigs2Scaffold)[Contig1];
						#pragma omp critical
						{
							Scaffolds2Reads[Scaffold].append(SeqBuffer[i]);
						}
					}
				}
			}
			// Paired reads. Note the i+=2 in for loop.
			else { 
				// Performance tweaks can be made in the schedule clause, maybe
				#pragma omp parallel for schedule(guided, 1000) num_threads(sub_max_threads)
				for (unsigned int i = 0; i < Buffer.size(); i+=2) {
					string Line1, Line2;
					unsigned int Contig1, Contig2;
					unsigned long ReadName1, ReadName2;
					unsigned int Scaffold, Scaffold2;
					stringstream ssScaffoldName;
					string ScaffoldName;
					string MatchFlag;
					set <string> ContigSet;
	
					stringstream tmp, tmp2;
					string Read1, Read2, Cont1, Cont2, trash;
	
					// Buffer of SAM file. First mapping in Line1, second in Line2.
					Line1 = Buffer[i];
		
					tmp << Line1;
	
					tmp >> ReadName1;
					tmp >> trash;
					tmp >> Cont1;
	
					Line2 = Buffer[i+1];
	
					tmp2 << Line2;
					tmp2 >> ReadName2;
					tmp2 >> trash;
					tmp2 >> Cont2;
	
					tmp.clear();
					tmp2.clear();
			
					Contig1 = 0;
					Contig2 = 0;
					char *end;
	

					if (Cont1[0] == '*') { 
						if (Cont2[0] == '*') { 
							if (output_unmatched) {
								UnmatchedBuffer[i]  = SeqBuffer[i];
								UnmatchedBuffer[i+1]  = SeqBuffer[i+1];
								UnmatchedBufferBool[i] = 1;
								UnmatchedBufferBool[i+1] = 1;
							}
							continue;
						}
						// If only Line2 has a valid mapping, then assign both reads to Contig2
						else { 
							Contig2 = strtol(Cont2.c_str(), &end, 10);
							Contig1 = Contig2; 
						}
					}
					// If only Line1 has a valid mapping, then assign both reads to Contig1
					else if (Cont2[0] == '*') { 
						Contig1 = strtol(Cont1.c_str(), &end, 10);
						Contig2 = Contig1; 
					}
					// Both Lines have valid mappings.
					else {
						char *end2;
						if (fill_gap_only && (Cont1 == Cont2)) { 
							continue; 
						}
						Contig1 = strtol(Cont1.c_str(), &end, 10);
						Contig2 = strtol(Cont2.c_str(), &end2, 10);
						// Make sure there are no non-digit characters in Cont2
						if (*end2!=0) { cerr << "Bad contig at line " << numLines << " with contig2 end character " << *end2 << endl; }
					}
					// Make sure there are no non-digit characters in Cont. This should probably be in the previous
					// logical block, not sure why it's here. Either way, it shouldn't matter.
					if (*end!=0) { cerr << "Bad contig at line " << numLines << " with contig end character " << *end << endl; }
	
					Scaffold = 0;
	
					// If contigs did not convert properly.
					if (Contig1 == 0) { 
						cout << "Bad contig name (" << Cont1 << ") or (" << Cont2 << ") (must be integers)" << endl;
	
						cout << "Line " << Line1 << endl;
						cout << "Line " << Line2 << endl;
						cout.flush();
						continue; 
					}
					// If the contigs don't map to any scaffolds
					if ((*Contigs2Scaffold)[Contig1] == 0) { 
						continue;
					}
	
	
					// If the reads align to the same contig
					if (Contig1 == Contig2) { 
						Scaffold = (*Contigs2Scaffold)[Contig1];
					}
					else {
						set <string> ContigSet, ContigSet2;
	
						Scaffold = (*Contigs2Scaffold)[Contig1];
						Scaffold2 = (*Contigs2Scaffold)[Contig2];
	
						// If the reads align to different scaffolds, we will discard the reads later.
						if (Scaffold != Scaffold2) {
							Scaffold = 0;
						}
					}

					// If both reads align to a single scaffold
					if (Scaffold) { 
						string tmp_string;
						tmp_string.assign(SeqBuffer[i]);
						tmp_string.append(SeqBuffer[i+1]);

						// Assign the reads to that scaffold.
						Scaffolds2ReadsSeq[i] = tmp_string;
						Scaffolds2ReadsScaffold[i] = Scaffold;
						tmp_string.clear();
					}
				}
			}
			// This consolidates the reads->scaffolds link from the previous loop block to Scaffolds2Reads
			for (unsigned int i = 0; i < Scaffolds2ReadsSeq.size(); i++) { 
				Scaffolds2Reads[Scaffolds2ReadsScaffold[i]].append(Scaffolds2ReadsSeq[i]);
			}
			// Obsolete/unsupported
			if (output_unmatched) {
				for (unsigned int i = 0; i < UnmatchedBufferBool.size(); i++) { 
					if (UnmatchedBufferBool[i]) { 
						Unmatched << UnmatchedBuffer[i];
					}
				}
				UnmatchedBuffer.clear();
				UnmatchedBufferBool.clear();
				UnmatchedBuffer.resize(MAX_BUFFER_SIZE);
				UnmatchedBufferBool.resize(MAX_BUFFER_SIZE, 0);
			}

			// Clear the working reads->scaffolds data structures
			Scaffolds2ReadsScaffold.clear();
			Scaffolds2ReadsSeq.clear();
			Scaffolds2ReadsScaffold.resize(MAX_BUFFER_SIZE);
			Scaffolds2ReadsSeq.resize(MAX_BUFFER_SIZE);
			Lines = 0;
			Buffer.clear();
			SeqBuffer.clear();
		}
	}
	cerr << endl;

	return Scaffolds2Reads;
}


void printUsage() {
	// Help file needs to be updated.
	cerr << "FinalFinishScaffolder -i <scaffold file> -c <contig file> -r <library file> -g <graph file> -t <number of files to process simultaneously> -m <threads per file>\n"  << endl;
}

unsigned int parse_contig(char *contig_file, unsigned int *scaffold_num, vector<unsigned int> *Contigs2Scaffold, 
						  vector<set <unsigned int> > *Scaffolds2Contigs, vector<string> *ContigSequences) { 
	ifstream fhContig(contig_file);
	if (!fhContig.is_open()) { 
		cerr << "Couldn't open contig file " << contig_file<< endl; exit(-1); 
	}
	cerr << "Processing contig file: " << contig_file << endl;

	unsigned long numLines = 0;
	unsigned long lastSize = 0;
	bool contigs_only_flag = 0;
	unsigned int tmp_scaff_num = *scaffold_num + 1;
	if (*scaffold_num == 0) { 
		contigs_only_flag = 1;
		cerr << "Outputting 1 contig per scaffold" << endl;
	}
	cerr << "Progress: ";

	while (fhContig.good()) {
		string Line1, Line2, Line3, Line4, WorkingLine;
		string dir;
		unsigned int ContigNum;
		unsigned int Scaffold;
		set<unsigned int> setContigs;
		getline (fhContig, Line1);
		if (!fhContig.good()) { break; }
		getline (fhContig, Line2);
		//cout << numLines << "\t";

		WorkingLine = Line1;
		WorkingLine.erase(0,1);

		ContigNum = atoi(WorkingLine.c_str());

		setContigs.clear();
	
		Scaffold = (*Contigs2Scaffold)[ContigNum];

		(*ContigSequences)[ContigNum] = Line2;




		if (Scaffold == 0) { 
			if (no_lib || contigs_only_flag) { 
				stringstream ssScaffoldName;
				string ScaffoldName;

				(*Contigs2Scaffold)[ContigNum] = tmp_scaff_num;
				setContigs.insert(ContigNum);
				(*Scaffolds2Contigs).push_back(setContigs);
	
				//(*scaffold_num) += 1;
			}
			//if (no_lib) { 
			else { 
				fstream OrphanContigs;
				OrphanContigs.open("contigs.orphan.fasta", fstream::out | fstream::app);
				OrphanContigs << Line1 << endl << Line2 << endl;
				OrphanContigs.close();
			}
			tmp_scaff_num++;
		}

		numLines++;

		if ((numLines % 100) == 0) {
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
	for (unsigned int k = 0; k < lastSize; k++) {
		cerr << "\b";
	}
	cerr << numLines;
	cerr << endl;
	*scaffold_num += tmp_scaff_num;

	return 1;
}

unsigned int parse_seq(string sequence_file, vector<string> *ReadSequences) {
	ifstream fhSeq(sequence_file.c_str());
	if (!fhSeq.is_open()) { 
		cerr << "Couldn't open sequence file " << sequence_file<< endl; exit(-1); 
	}
	//cerr << "Processing sequence file: " << sequence_file << " (will take some time)" << endl;
//	cerr << "Progress: ";
//	unsigned long numLines = 0;
//	unsigned long lastSize = 0;
	while (fhSeq.good()) {
		string Line1, Line2, Line3, Line4, WorkingLine;
		//unsigned int Scaffold;
		unsigned long SeqNum;
		getline (fhSeq, Line1);
		if (!fhSeq.good()) { break; }
		getline (fhSeq, Line2);

		WorkingLine = Line1;
		WorkingLine.erase(0,1);
		
		SeqNum = atol(WorkingLine.c_str());
		(*ReadSequences)[SeqNum].assign(Line2);
/*
		numLines++;

		if ((numLines % 1000000) == 0) {
			for (unsigned long k = 0; k < lastSize; k++) {
				cerr << "\b";
			}
			stringstream Output;
			string strOutput;
			Output << numLines;
			Output >> strOutput;
			lastSize = strOutput.size();
			cerr << numLines;
		}
		*/
		//cout << WorkingLine << endl;
	}
	
	return 1;
}

unsigned int parse_scaffold_file(char *scaffold_file, unsigned int *scaffold_num, vector<unsigned int> *Contigs2Scaffold,
								 vector<set <unsigned int> > *Scaffolds2Contigs, vector<string> *ScaffoldStrings) { 
	ifstream fhScaffold(scaffold_file);
	if (!fhScaffold.is_open()) { 
		cerr << "Couldn't open scaffold file " << scaffold_file<< endl; exit(-1); 
	}
	*scaffold_num = 1;
	std::ios::sync_with_stdio(false);


	cerr << "Reading scaffold file .." << endl;
	while (fhScaffold.good()) {
		string Line, ScaffoldName;
		//stringstream ssScaffoldName;
		set <unsigned int> Contigs;
		getline (fhScaffold, Line);
		if (!fhScaffold.good()) { break; }
		vector<string> SplitContigs;

		Contigs.clear();
		//vector<unsigned long> ulSplitContigs;

		//ssScaffoldName << scaffold_num++;
		//ssScaffoldName >> ScaffoldName;

		boost::split(SplitContigs, Line, boost::is_any_of(" \t"));
		for (unsigned int i = 0; i < SplitContigs.size(); i++) {
			unsigned int uiContig = atoi(SplitContigs[i].c_str());
			(*Contigs2Scaffold)[uiContig] = *scaffold_num;
			Contigs.insert(uiContig);

		}

		(*Scaffolds2Contigs).push_back(Contigs);
		

		//for (j = 0; j < SplitContigs.size(); j++) { 
		//	Contigs.insert(atol(SplitContigs[j].c_str()));
		//}
		(*ScaffoldStrings)[*scaffold_num].append(Line);
		*scaffold_num += 1;
	}
	return 1;
}



unsigned int parse_graph(char *GraphFileLocation, vector<string> *Graph) { 
	cerr << "Reading graph file.." << endl;

	ifstream fhGraph(GraphFileLocation);
	if (!fhGraph.is_open()) { 
		cerr << "Couldn't open graph file " << GraphFileLocation<< endl; exit(-1); 
	}

	while (fhGraph.good()) {
		string Line;
		string trash;
		stringstream ssLine;
		unsigned int contig;

		getline (fhGraph, Line);
		if (!fhGraph.good()) { break; }

		if (Line[0] == 'e') { 
			ssLine << Line;
			ssLine >> trash;
			ssLine >> contig;
			(*Graph)[contig].append(Line);
			(*Graph)[contig].append("\n");
		}
		else {
			(*Graph)[contig].append(Line);
			(*Graph)[contig].append("\n");
		}
	}
	return 1;
}

string create_working_directories (unsigned int LastDir) {
	string dir, dir2, dir3, dir4, dir5, dir6;
	stringstream ddir, ddir2, ddir3, ddir4, ddir5, ddir6;
	unsigned int ScaffoldDir1, ScaffoldDir2, ScaffoldDir3, ScaffoldDir4, ScaffoldDir5;
	struct stat st;

	if (stat("working",&st) != 0) {
		mkdir("working",0777);
	}
	for (unsigned int i = 1; i <= LastDir; i++) { 
		ScaffoldDir1 = i % 10;
		ScaffoldDir2 = i % 100;
		ScaffoldDir3 = i % 1000;
		ScaffoldDir4 = i % 10000;
		ScaffoldDir5 = i % 100000;

		ddir.clear();
		ddir2.clear();
		ddir3.clear();
		ddir4.clear();
		ddir5.clear();
		ddir6.clear();


		ddir6 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2 << "/" << ScaffoldDir3 
		   	  << "/" << ScaffoldDir4 << "/" << ScaffoldDir5 << "/" << i;
		ddir6 >> dir6;
		if (stat(dir6.c_str(),&st) == 0) { continue; }
	//	if (stat(dir6.c_str(),&st) == 0) { return dir6; }


		ddir << "working/" << ScaffoldDir1;
		ddir >> dir;
	
		ddir2 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2;
		ddir2 >> dir2;

		ddir3 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2 << "/" << ScaffoldDir3;
		ddir3 >> dir3;

		ddir4 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2 << "/" << ScaffoldDir3 
			  << "/" << ScaffoldDir4;
		ddir4 >> dir4;
	
		ddir5 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2 << "/" << ScaffoldDir3 
			  << "/" << ScaffoldDir4 << "/" << ScaffoldDir5;
		ddir5 >> dir5;


		if (stat(dir.c_str(),&st) != 0) {
			mkdir(dir.c_str(),0777);
		}
		if (stat(dir2.c_str(),&st) != 0) {
			mkdir(dir2.c_str(),0777);
		}
		if (stat(dir3.c_str(),&st) != 0) {
			mkdir(dir3.c_str(),0777);
		}
		if (stat(dir4.c_str(),&st) != 0) {
			mkdir(dir4.c_str(),0777);
		}
		if (stat(dir5.c_str(),&st) != 0) {
			mkdir(dir5.c_str(),0777);
		}
		if (stat(dir6.c_str(),&st) != 0) {
			mkdir(dir6.c_str(),0777);
		}
		if (stat(dir6.c_str(),&st) != 0) {
			cerr << "Couldn't build directory: "<< dir6 << endl;
			dir6.assign("error");
		}
	}
	return dir6;
}

string get_working_directory (unsigned int LastDir) {
	string dir, dir2, dir3, dir4, dir5, dir6;
	stringstream ddir, ddir2, ddir3, ddir4, ddir5, ddir6;
	unsigned int ScaffoldDir1, ScaffoldDir2, ScaffoldDir3, ScaffoldDir4, ScaffoldDir5;
	struct stat st;


	ScaffoldDir1 = LastDir % 10;
	ScaffoldDir2 = LastDir % 100;
	ScaffoldDir3 = LastDir % 1000;
	ScaffoldDir4 = LastDir % 10000;
	ScaffoldDir5 = LastDir % 100000;


	ddir6 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2 << "/" << ScaffoldDir3 
	   	  << "/" << ScaffoldDir4 << "/" << ScaffoldDir5 << "/" << LastDir;
	ddir6 >> dir6;
	if (stat(dir6.c_str(),&st) == 0) { return dir6; }


	ddir << "working/" << ScaffoldDir1;
	ddir >> dir;

	ddir2 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2;
	ddir2 >> dir2;

	ddir3 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2 << "/" << ScaffoldDir3;
	ddir3 >> dir3;

	ddir4 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2 << "/" << ScaffoldDir3 
		  << "/" << ScaffoldDir4;
	ddir4 >> dir4;
	
	ddir5 << "working/" << ScaffoldDir1 << "/" << ScaffoldDir2 << "/" << ScaffoldDir3 
		  << "/" << ScaffoldDir4 << "/" << ScaffoldDir5;
	ddir5 >> dir5;


	if (stat(dir.c_str(),&st) != 0) {
		mkdir(dir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);
	}
	if (stat(dir2.c_str(),&st) != 0) {
		mkdir(dir2.c_str(),0777);
	}
	if (stat(dir3.c_str(),&st) != 0) {
		mkdir(dir3.c_str(),0777);
	}
	if (stat(dir4.c_str(),&st) != 0) {
		mkdir(dir4.c_str(),0777);
	}
	if (stat(dir5.c_str(),&st) != 0) {
		mkdir(dir5.c_str(),0777);
	}
	if (stat(dir6.c_str(),&st) != 0) {
		mkdir(dir6.c_str(),0777);
	}
	if (stat(dir6.c_str(),&st) != 0) {
		cerr << "Couldn't build directory: "<< dir6 << endl;
		dir6.assign("error");
	}


	if (stat(dir6.c_str(),&st) != 0) {
		cout << "Couldn't find directory: "<< dir6 << endl;
	}

	return dir6;
}

unsigned long get_last_read_num(string sam_file) {
	unsigned long last_read_num = 0;

	ifstream file(sam_file.c_str());
	if (!file.is_open()) { 
		cerr << "Couldn't open file: " << sam_file << endl;
		exit(-1);
	}

	file.seekg(0, ios::end);
	file.unget();
	file.unget();
	char c = '\0';
	while (c != '\n') { 
		file.unget();
		file.unget();
		c = file.get();
	}
//	file.unget();
	
	string line;
	getline(file,line);

	//cout << "Line " << sam_file << " " << line << endl;
	stringstream tmp;
	tmp << line;
	tmp >> last_read_num;


	return last_read_num;
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

int parse_graph_for_scaffolds (
    vector <string> *graph_strings, 
    map <string, unsigned int> *counts, 
    map <string,int> *distances, 
    map <string, unsigned int> *first_direction, 
    map <string, unsigned int> *second_direction) { 

	int dist, direct1, direct2, count;
	vector <string>::iterator gs_it;
	string edge;

	for (gs_it = graph_strings->begin(); gs_it != graph_strings->end(); gs_it++) { 
		istringstream lines(*gs_it);
		string line;
		while (getline(lines, line)) { 
			string tag;
			stringstream linestream, linestream2;
			linestream << line;
			linestream >> tag;
			string contig1, contig2;
			
			 
			switch (tag[0]) {
				case 'e':
					linestream >> contig1;
					linestream >> contig2;
					linestream >> count;
	
					linestream2 << contig1 << "-" << contig2;
					linestream2 >> edge;
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
						counts->insert(pair <string, int> (edge, count));	
						distances->insert(pair <string, int> (edge, dist));	
						first_direction->insert(pair <string, unsigned int> (edge, direct1));	
						second_direction->insert(pair <string, unsigned int> (edge, direct2));	
					}
					edge.clear();
					break;
				default:
					cerr << "Flag " << tag << " not recognized." << endl;
					break;
			}
		}

	}
	return 0;
}

string build_scaffold (map <string, unsigned int> *counts, 
					map <string,int> *distances, 
					map <string, unsigned int> *first_direction, 
					map <string, unsigned int> *second_direction,
					string scaff_string,
					vector <string> *contig_seqs) {

	string FinalSequence;

	bool FirstContig = 1;
	string Sequence, edge;

	unsigned int last_contig;

	stringstream scaffolds(scaff_string);
	scaffolds >> last_contig;
	stringstream Scaffold_Acc;
	Scaffold_Acc << ">Scaffold-" << last_contig;

	while (scaffolds.good()) {
		unsigned int current_contig;
		int first_direction_val, second_direction_val;
		scaffolds >> current_contig;
		Scaffold_Acc << "_" << current_contig;

	
		Sequence = (*contig_seqs)[current_contig];
	

		stringstream tmp;
		map <string, unsigned int>::iterator first_direction_it;
		map <string, unsigned int>::iterator second_direction_it;
		map <string, int>::iterator distances_it;

		unsigned int contig1 = last_contig;
		unsigned int contig2 = current_contig;


		if (contig1 < contig2) { 
			tmp << contig1 << "-" << contig2;
			tmp >> edge;
			first_direction_it = first_direction->find(edge);
			second_direction_it = second_direction->find(edge);
			if ((first_direction_it == first_direction->end()) || second_direction_it == second_direction->end()) {
				cerr << "edge " << edge << endl;
				throw "Bad edge ";
			}
		}
		else { 
			tmp << contig2 << "-" << contig1;
			tmp >> edge;
			// Switching first and second here is necessary to maintain the ordering
			// information of the contigs
			first_direction_it = second_direction->find(edge);
			second_direction_it = first_direction->find(edge);
			if ((first_direction_it == first_direction->end()) || second_direction_it == second_direction->end()) {
				cerr << "edge " << edge << endl;
				throw "Bad edge";
			}
		}

		first_direction_val = first_direction_it->second;
		second_direction_val = second_direction_it->second;

		if (FirstContig) { 
			if (first_direction_val == 0) {
				FinalSequence = rev_comp((*contig_seqs)[last_contig]);
			}
			else { 
				FinalSequence = (*contig_seqs)[last_contig];
			}
			FirstContig = 0;
		}

		if (second_direction_val == 0) {
			Sequence = (*contig_seqs)[current_contig];
		}
		else { 
			Sequence = rev_comp((*contig_seqs)[current_contig]);
		}

		distances_it = distances->find(edge);

		int gap = distances_it->second;
		if (gap >= 0) { 
			for (int j = 0; j <= gap; j++) { 
				FinalSequence.push_back('N');
			}
		}
		else if (gap < -10) {
			
			unsigned int positive_gap = 0 - gap;

			unsigned int gap_seed_size = positive_gap / 4;
			unsigned int gap_seed_offset = positive_gap * 0.75;
			if (gap_seed_size > 30) { 
				gap_seed_size = 30; 
			}
			if (positive_gap < 20) { 
				gap_seed_size = positive_gap/2;
				gap_seed_offset = gap_seed_size;
			}
			string seed_seq;
			for (int erase_offset = gap_seed_offset - gap_seed_size; erase_offset > 0; erase_offset--, gap_seed_offset--) {
				if (positive_gap < FinalSequence.size()) {
					seed_seq.assign(FinalSequence, FinalSequence.size() - gap_seed_offset, gap_seed_size);

					for (unsigned int k = 0; k < positive_gap * 2; k++) { 
						string seed_seq2;
						seed_seq2.assign(Sequence, k, gap_seed_size);
						if (seed_seq2.size() < gap_seed_size) { 
							positive_gap = 0;
							gap = 0;
						}
						else if (seed_seq == seed_seq2) { 
							Sequence.erase(0,k+gap_seed_size);
							FinalSequence.erase(FinalSequence.size() - erase_offset, erase_offset);
							gap = 0;
							positive_gap = 0;
							erase_offset = 0;
						}
					}
				}
			}
			if (gap < 0) { 
				unsigned int clip = int(1 + positive_gap/2);
				Sequence.erase(0,clip);
				FinalSequence.erase(FinalSequence.size() - clip, clip);
				FinalSequence.append("NN");
			}
		}
		FinalSequence.append(Sequence);

		last_contig = current_contig;
	}
	string tmp;
	Scaffold_Acc >> tmp; 
	tmp.push_back('\n');
	
	FinalSequence.insert(0,tmp);
	FinalSequence.push_back('\n');
	return FinalSequence;
}



