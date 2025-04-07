
README for Sex Determination of the two-spotted goby Pomatoschistus flavescens
--------------------------------------------------------------------------------
%¤\¤%\%¤\¤%\%¤\¤%\o/%¤/%¤/%¤/%¤/%¤/%¤/o\%¤\¤%\%¤\¤%\%¤\¤%\o/%¤/%¤/%¤/%¤/%¤/%¤/o\

List of data files and scripts/code, organised by section of the project.

1. Analysis of field data of Adult Sex Ratio (ASR)
--------------------------------------------------

TO DO

2. Identification of sex-specific loci using RADseq data aligned to the reference genome.
-----------------------------------------------------------------------------------------

Data:
> ddRAD reads will be accessible at: http://www.ncbi.nlm.nih.gov/bioproject/1224570
> reference genome will be accessible at: XXXX
> popmap_hap1.txt      # Population map, indicates the sex and population of origin of each individual in the ddRAD dataset. For alignment to haplotype 1 of the reference genome.
> popmap_hap2.txt      # Population map, indicates the sex and population of origin of each individual in the ddRAD dataset. For alignment to haplotype 2 of the reference genome.
> popmap_radsex.tsv    # Population map used to identify the sex of individuals in the Radsex pipeline

Code:
> index_flav.sh      # the reference genome is indexed before alignment of ddRAD reads
> align_flav.sh      # Uses the ddRAD reads and the indexed reference genome and produces sorted alignment files, for haplotypes 1 and 2 of the reference genome.
> gstacks_flav.sh    # Used the sorted alignment files and the population map and feeds them to the Stacks pipeline, first step gstacks to produce a catalog.
> population_flav.sh # Uses the catalog and the population map to produce a VCF file.
> assoc_flav.sh      # Filters the vcf, then uses plink to test for assosciation of allele frequency and SNP missingness with sex
> radsex_flav.sh     # Runs the Radsex pipeline on the ddRAD reads, finds sex specific markers, maps them to the reference genome
> markers_list.R     # Extracts the list of markers associated with sex from the output of Stacks and Radsex


3. Identification of the sex of eggs and larvae from genotype
--------------------------------------------------------------

TO DO