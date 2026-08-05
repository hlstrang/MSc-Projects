# BIOL60201 Programming Skills
## University of Manchester MSc Bioinformatics and Systems Biology
### Project description:

microRNAs (miRNAs) are short non-coding RNA that regulate gene expression. They achieve this by binding to loci within a target messenger RNA (mRNA) to degrade it, to destabilise it, or to repress its subsequent translation. It has been found that a single miRNA can target multiple mRNAs/genes, and this pleiotropic nature is of interest to pharmaceutical companies as drug targets for diseases with multifactorial origin.

To begin identifying a candidate miRNA for therapy, experiemnets can be done to identify its target genes and subsequently identify diseases associated with that gene. Various databases or repositories store information on such associations. These databases can be set up using various methodologies, such as manual curation with input from literature or computational curation with predicted targets of association between miRNA and genes as an example.

The dataset provided contains several columns of information, but the information of interest would be gene, disease, miRNA, and locus information. The goal is to output basic statistics before and after quality control for confidence scores alongside number of PMIDs and number of SNPs. Also include the counts of the top 3 organisms following quality control. There should also be figures generated post QC.

### Task:

The script should be able to handle the following:

    - Argument to read in csv file
    - Argument to output the tsv file
    - Argument to accept two parameters for quality control
    - Arguments to generate figures
    - When reading in the file, basic data descriptors should be printed out
    - After quality control, basic data descriptors should be printed out
    - Argparse should have written help to describe simply how to execute the script and its basic functionality

### Grade achieved: 100
