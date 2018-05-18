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



