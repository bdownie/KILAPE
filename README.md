#LICENSE
The software was developed at the Leibniz Institute on Aging - Fritz Lipmann Institute (FLI; http://www.leibniz-fli.de/) under a mixed licensing model. Researchers at academic and non-profit organizations can use KILAPE using the included license, while for-profit organizations are required to purchase a license (http://genome.leibniz-fli.de/software/index.php). By downloading the package you agree with conditions of the FLI Software License Agreement for Academic Non-commercial Research (LICENSE.pdf).

# KILAPE
K-masking and Iterative Local Assembly of Paired Ends

KILAPE (K-masking and Iterative Local Assembly of Paired Ends) is an automated scaffolding and gap filling software pipeline which predicts repetitive elements in Next Generation Sequencing read libraries without resorting to a reference sequence.  The package of KILAPE consists of pre-compiled C++ program modules as well as a set of Perl scripts for data preparation and for processing within the pipeline itself.

==================================================================================================

INSTALLATION:
Extract the downloaded archive to a local directory. Currently, KILAPE binaries have been compiled for standard 64-bit Linux distributions. 

1) tar xzf kilape_v1.0b.tar.gz

2) Copy files to installation directory. Be sure to maintain directory structure (i.e. "bin" directory
   containing binaries in the directory containing KILAPE.pl).

A sample configuration file along with explanations can be found in the "demo" directory.

To test your installation:

1) cd demo

2) Edit kilape.conf to account for binary locations (e.g. path to bowtie2)

 - scaffolds.demo.fasta should be similar/identical to scaffolds.velvet.fasta (depending on velvet version)
 
3) ../KILAPE.pl

==================================================================================================

DEPENDENCIES

REQUIRED:

Boost Libraries for C++ (www.boost.org)

Jellyfish (http://www.cbcb.umd.edu/software/jellyfish/)

BWA (http://bio-bwa.sourceforge.net/), Bowtie (http://bowtie-bio.sourceforge.net/index.shtml), or Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

OPTIONAL (recommended):

Parallel::ForkManager perl module (http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm)
- Without this, perl scripts will need to be altered and local assembly will not be run in parallel.
- Also, mapping of read libraries will not be performed in parallel.

velvet (http://www.ebi.ac.uk/~zerbino/velvet/)

wgs-assembler (http://wgs-assembler.sourceforge.net)
- Required for local assembly functionality
- Local assembly function is currently in beta and not supported.
