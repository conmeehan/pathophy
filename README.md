# pathophy
Scripts for aiding in pathogen phylogenetics analysis

### <b>Phylogenetic Inclusion python script:</b><br/>
This script is used to take a transmission cluster (e.g. defined by SNP cut-offs) and extend it to include all other related strains that share the common ancestor in a phylogenetic tree.<br/>
####Inputs:<br/>

A cluster file: one cluster per line, separated by tabs<br/>
A tree file in newick format<br/>
NOTE: requires dendropy to be installed<br/>
usage<br/>
python phylogeneticInclusion.py -t treeFile -c clusterFile<br/>

See comments in script for full instructions
<br/>

### <b>Mixed population detection by comparing VCF files</b><br/>
This script is used to try estimate mixed populations from whole genome sequencing VCF files.

#### Input<br/>
2 vcf files output from a SNP calling pipeline such as snippy but at different minumum proportion for variant evidence (fractions)<br/>
E.g. run snippy with --minfrac 0.3 and compare to default (0.9)

See comments in script for full instructions
<br/>

### <b>Snippy output to Phylogenetic input</b><br/>
This pipeline converts outputs from snippy into a whole genome alignment, a SNP alignment and an invariant count of positions for use in ascertainment bias correction in phylogenetics

#### <b>Input:</b><br/>
A folder in which there are snippy output folders. In particular the following files must be in each<br/>
	snps.vcf	This is used to get all snps and indels<br/>
	snps.aligned.fa	This is used for the Ns and large areas of either deletions or unaligned regions (e.g. insertion sites; multiple mapping sites). These two latter cannot be distinguished within this script.<br/>
	These files can be named other things (e.g. not snps but another prefix, as done by the --prefix in snippy) but must be PREFIX.snps and PREFIX.aligned.fa with no other aligned.fa file<br/>
The fasta reference sequence used to create the snippy outputs (assumed same for all). E.g. this is the ref.fa file in the reference folder of a snippy output.<br/>
An interval masking file (optional; see below)<br/>
A set of characters not to count as variable (optional; see below)<br/>
A flag for whether to include the refrence sequence (optional)<br/>

See comments in script for full instructions

### <b>Cutting BEAST output to desired step</b><br/>
Script takes a MCMC log or trees file from BEAST and will cut them at a specific step
E.g. the file may have more than 40,000 steps and you wish to cut at that step (step requested will be included)

#### Input<br/>
File to cut
File type (log or trees)
Last step to include

#### Output<br/>
New log or trees file. Will be named same as input with _cut<number>
e.g. if input it Dataset.log and the last step is 40000, the output is Dataset_cut40000.log

#### Usage<br/>
python cut_BEAST_output.py --file FileToCut --type <log/trees> --step <finalStepToInclude>

### <b>Creating core genome alignments from isolate lists of aligned genes</b><br/>
This consists of two scripts: isolateFilesToGeneFiles.py and genesToConcatAlign.py

#### isolateFilesToGeneFiles.py
This script is used to combine files of individual gene sequences, per isolate, into gene-specific alignment files

##### Input:
A folder of isolate files where the isolate name is everything before the first _.
NOTE: assumes all end in .fasta
Inside each file is a set of fasta formatted sequences where the gene name is xxx_xxx (i.e. based around 1st _ only)
e.g. >STENO_1_6, the gene name is STENO_1

##### Output:
In a folder named geneAlignments:
Per gene files where the file name is the gene name (e.g. STENO_1.fasta)
Inside each file the sequences are named as isolate_geneName (e.g. 677_STENO_1)

NOTE: uses biopython
NOTE: assumes only 1 entry per gene in each isolate. Will warn if double and only keep last one

##### Usage:
python isolateFilesToGeneFiles.py --folder folderName

#### genesToConcatAlign.py

This script is used to combine files of gene alignments into a concatenated alignment (assumes all species/isolates in all gene alignments)

##### Input:
A folder of fasta gene alignmentds where the sequence names is the species/isolate name (identical between files)
NOTE: assumes all end in .fasta

##### Output:
A concatenated alignment of all genes, using the species name as the sequence header (pulled from the first input file)
A file that outlines the gene name and the start/stop positions of the gene in the concatenated alignment
e.g. STENO_1, 1-1345

NOTE: uses biopython
NOTE: Assumes all genes have all species/isolates. If not, will mess up the alignment

##### Usage:
python genesToConcatAlign.py --folder folderName


### <b>Creating SNP cut off-based clusters from an alignment</b><br/>
This consists of two scripts: pairwiseDist.py and pairwiseCluster.py
NOTE: inputs to these scripts are position specific so must be given in a specific order

#### pairwiseDist.py
Script to take an alignment and get pairwise distance between all sequences as a matrix, with allowed exclusions
##### Input:
A fasta alignment (fastaAlignment)
Characters to not count as divergent if encountered (excludedCharacters)
e.g. if the excluded characters are N- then if a comparison has a N or - this will not be counted as different
NOTE: is not case sensitive so N and n are the same thing

##### Output:
A mnatrix of distances between the samples.

##### Usage:
python pairwiseDist.py fastaAlignment excludedCharacters

#### pairwiseCluster.py
This script is used to cluster sequences based on a pairwise distance file. This takes the output of pairwiseDist.py as input

##### Input:
The pairwise distance matrix output from pairwiseDist.py (pairwiseDistFile)
A SNP cut-off number (cut-off)
##### Output:
A file listing all sequences that are within the cut-off from another sequence
A file where each line is a list of sequences within cut-off distance of each other (tight cluster)
A file where each line is a list of sequences where at each sequence is within the cut-off from at least 1 other (e.g. if cutoff is 5 then seq1 and seq2 could be 5 apart and seq3 could be 10 from seq2 and 3 from seq1 but all on same line) (loose cluster)

NOTE: cut-off is used as <= so if cutoff is 5 then those a distance of 5 will be counted as clustered
NOTE2: NA can be an entry (e.g. for a tip distance that is not supported by bootstrap or unequal lengths of patterns) so if it is, count as not within cutoff

##### Usage:
python pairwiseCluster.py pairwiseDistFile cut-off


### <b>Detecting gene absences from MTBseq outputs</b><br/>
Script takes a folder of Position Tables created by MTBseq and searches for genes that have no coverage (suggesting they are absent).

#### Input<br/>
Folder containing Position Tables from MTBseq (default is Position_Tables)
The genes file from the reference used for MTBseq (can be found in the var/ref/ folder of the MTBseq installation)
The minimum proportion of positions absent to be counted as an absent gene (default is 0.95)
The minimum number of reads to consider a position to be covered (default is 8, same as MTBseq)

#### Output<br/>
A presence/absence map of the genes in a matrix. Names are the position table names with everythign after the first _ removed (similar to MTBseq)
The samples are listed in the first column; the genes are listed in the first row. There is a 1 for a present gene and a 0 for an absent gene

#### Usage<br/>
python MTBseq_to_genePresAbs.py --folder Position_TablesFolder --genes genesFile --cutoff proportionAbsent --minreads coverageMinumum

