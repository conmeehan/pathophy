# pathophy
Scripts for aiding in pathogen phylogenetics analysis

<b>Phylogenetic Inclusion python script:</b><br/>
This script is used to take a transmission cluster (e.g. defined by SNP cut-offs) and extend it to include all other related strains that share the common ancestor in a phylogenetic tree.<br/>
Inputs:<br/>

A cluster file: one cluster per line, separated by tabs<br/>
A tree file in newick format<br/>
NOTE: requires dendropy to be installed<br/>
usage<br/>
python phylogeneticInclusion.py -t treeFile -c clusterFile<br/>

<b>Mixed population detection by comparing VCF files</b><br/>
This script is used to try estimate mixed populations from whole genome sequencing VCF files.

Input is 2 vcf files output from a SNP calling pipeline such as snippy but at different minumum proportion for variant evidence (fractions)<br/>
E.g. run snippy with --minfrac 0.3 and compare to default (0.9)

Output are the SNPs (genome position, ref base, alt base) which are present in the lower fraction cut-off vcf but not the upper

User may supply a specific set of positions to output into a separate filtered list. These may be drug resistant sites, lineage defining SNPs etc. original output file will still contain these, these positions willm just appear in both files.<br/>
Each line should be a different should be start of interval<tab character>end of interval per line (inclusive). Alternatively, if only one position on a line, that single site will be filtered.


NOTE: in theory, this should work with all vcf files but only has been tested with those output from snippy<br/>
NOTE: only SNPs are output. Indels, mnp and complex are skipped (for now)<br/>
Usage:<br/>
python compareVcfMixedPopulations.py --lower vcfLowerFraction --upper vcfUpperFraction --filter positionFile (optional)

