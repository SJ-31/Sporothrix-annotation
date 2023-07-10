# Overview
This repository contains the code used in the paper ... It consists of Nextflow workflows for genome and transcriptome assembly, then subsequent genome annotation with maker followed by variant calling with GATK. The maker and variant calling pipelines were adapted from code by Card et al. (2019) and Khalfan (2020) respectively. This file provides a brief overview of each workflow steps, their relevant shell commands, the tools involved and the external data used (links included in the references). More on the rationale behind each steps (without the shell commands), such as why a certain tool was used, can be found [here](./info/README_Extended.md)

## Genome/transcriptome assembly
- **Input:**
  - Raw paired-end fastq files
  - Reference genome fasta and gff for *Sporothrix schenckii*, GCF_000961545.1
  - mtDNA reference
- 1. FastQC, MultiQC: Evaluate raw paired-end fastq files
    * FastQC provides information about read quality, overrepresented sequences, adapter content and other statistics on a per-fastq-file basis; MultiQC collates the information from multiple FastQC report files for convenience.
```bash
fastqc <reads1> <reads2>
multiqc <directory_containing<reads>
```
- 2. FastP: Trim adapters and poly g (identified as a contaminant in prior step).
    * FastP is an all-purpose quality control tool that can perform adapter removal, read trimming and base correction
```bash
fastp -i <reads1> -I <reads2> \
    --detect_adapter_for_pe \ # Detect adapters in paired-end reads
    -c \ # Enable overlap analysis for paired-end reads, which corrects mismatched bases if the overlap meets a certain length (default 30) and the mistmach exceeds a certain number of bases (default 5)
    -Q \ # Disable phred-based quality filtering - the reads for all given samples were already high quality (> 25).
    --trim_poly_g \ # Trim poly-g tails
    -o <reads1_out> -O <reads2_out> # Specify names for the cleaned output reads
```
- 3. BBduk: Filter mtDNA
  - A kmer-based quality control tool, BBduk was used here for the function of filtering out reads that match strongly to a set of query sequences (mtDNA in this case).
  - The kmer size setting for BBduk was set to the maximum of `k=31` to make filtering as stringent as possible
```bash
bbduk.sh in1=<reads1> in2=<reads2> \
    out1=<reads1_out> out2=<reads2_out> \
    ref=<reference> \ # The reference sequences to filter out from, in fasta format
    outm=<flagged> \ # Name for file containing flagged reads that got matched to ref
    k=31 \
    -Xmx1G -Xms16M \ # Allocate memory to prevent "out of memory" issues
    ```
- 4. Spades, Megahit: Assemble reads
  - Both ran in default settings
```bash
spades.py \
    -o <sample_name> \ # The directory to store resulting files.
    # The final assembly is denoted "scaffolds.fasta"
    --pe-1 <reads1> --pe-2 <reads2> \
megahit \
    -o <sample_name> # The directory to store resulting files
    # The final assembly is denoted "final.contigs.fa"
    -1 <reads1> -2 <reads2>
```
- 5. BUSCO, Quast: Assess assembly quality
  - BUSCO was set to use the `sordariomycetes_odb10` lineage in `genome` mode, while Quast was provided with the reference genome
      - BUSCO needs to be set in `offline` mode or else it will downloads its lineage database automatically in the directory of each run
      - Quast can assess several assemblies at once for comparison purposes
```bash
# First manually download the sordariomycetes_odb10 lineage dataset
busco --download sordariomycetes_odb10

busco -i <assembly_fasta> \
    -l sordariomycetes_odb10 \
    -o <busco_output> \ # The busco output directory
    -m genome \
    --download_path <path> # The path to the "busco_downloads" folder obtained above

quast.py <assemblies> \ # The paths to all assemblies
    --fungus \ # Specifies to analyze as a fungal genome
    -r <reference> \ # path to reference genome fasta file
    -g # path to reference genome gff file
```
  - Transcriptome assembly was performed in much the same way as genome assembly, except mtDNA was not filtered, RNAspades was used as the assembler and finally BUSCO was set in `transcriptome` mode
- 6. Ragout:  Assemble contigs into scaffolds
  - Ragout was provided with the GCF_000961545.1 reference sequence to scaffold from, and run on default settings
      * The paths to files are given to Ragout with an `rcp` file, (an example can be found [here](./info/recipe.rcp)). A python script was used to automate writing this file from the command line
```bash
ragout recipe.rcp \
-o <output_directory> # Name of the output directory; the final scaffolds will be found here under the name of "target" specified in the rcp file
```
- 7. Minimap2: Align scaffolds to reference sequence
```bash
minimap2 -a <reference> <scaffolds> > aligned.sam # Align the scaffolds
```
- 8. Samtools: Extract scaffolds from sam file into separate fasta files for each chromosome
```bash
samtools sort aligned.sam -o aligned.bam # Sort, then index the aligned bam file
samtools index aligned.bam # Necessary for using samtools view with a range specifier
samtools view -h aligned.bam <chromosome> | samtools fasta > <chromosome>.fasta
# "chromosome" specifies which chromosome is to be extracted, which requires adding the -h (header) flag
# The pipe to "samtools fasta" extracts sequences aligned to that chromosome from the aligned.bam file in fasta format
```
- **Output:** genome assemblies for each sample, split into fasta files for each chromosome present in the reference sequence

## Genome annotation
- **Input:**
  - Genemarks hmm files for each sample, trained in-house with GenemarksES
  - Augustus species model trained with the Augustus web server on the GCF_000961545.1 reference
  - A database of repeat elements, generated by combining the TREP database with models produced in-house by RepeatModeler on several *Sporothrix* species
  - est evidence: includes the previously assembled transcriptome in addition to known *Sporothrix* genes collected from the NCBI nucleotide database with query `txid29907[Organism:exp] AND biomol_mrna[PROP]`
  - Repeat proteins provided with the Maker download
  - gff evidence from the GCF_000961545.1 reference
  - Protein evidence collected from UniProt Release 2023_02 with query `Sporothrix`
```bash
# Train Genemarks
gmes_petap.pl -ES -sequence <scaffolds> -fungus

# Predict repeat elements from reference fasta file
BuildDatabase -name <species>_db <reference>
RepeatModeler -database <species>_db -LTRStruct
# The -LTRStruct flag runs the LTR structural discovery pipeline in addition to the standard discovery setting
# Concatenate the RepeatModeler output (the file called consensi.fa.classified) from each species into a single file to be used by maker
```
- 1. Maker: First round of Maker annotation, using est, gff, protein, repeat evidence and Genemarks model
    * Example of options [here](./info/maker_opts.ctl_Round1), a comment with "##" marks where an option has been changed for this round
```bash
maker -CTL # Generate the maker control files in the current directory
# Manually edit maker_opts.ctl following the example above, or use the script "maker_cli.py" found in the bin directory
maker -fix_nucleotides # Run maker in the directory containing the control files. Need the "-fix_nucleotides" flag to be compatible with some of the evidence
gff3_merge -d <maker_output_directory>/<round1>_master_datastore.log # Combine all the maker annotations from each scaffold into a single gff file
```
- 2. SNAP: Train SNAP on the output from the first round of annotation.
```bash
maker2zff <round1_gff> # This maker script converts the gff output into a zff file (<round1>.ann) and extracts the sequences in fasta format (<round1>.dna)
fathom genome.ann genome.dna -gene-stats > ${name}_SNAP-gene-stats.log 2>&1 # Obtain gene log file for statistics
fathom genome.ann genome.dna -validate > ${name}_SNAP-validate.log 2>&1 # Obtain log file for gene validations.
fathom -categorize 100 <round1>.ann <round1>.dna # Categorize genes for training: genes with errors overlapping genes etc.
#   The number after specifies how much intergenic sequence to place on either side of the gene
fathom -export 100 -plus uni.* # Only the unique genes will get exported for training
forge export.ann export.dna
hmm-assembler.pl <round1> . > <round1>.snap_hmm # Build the snap model from the files in the directory
```
- 3. Maker: Second round of annotation, disabling all evidence-based annotation and using Snap and Augustus for gene prediction only. This was set up so Maker would  populate its annotations with the gff file from the first round
    * Example of options [here](./info/maker_opts.ctl_Round2)
- 4. Snap: Second round of SNAP training, using the Maker's second round output
    * The commands are the same as used in training SNAP in the first round
- 5. Maker: Final round of annotation with the most recent SNAP model
    * The maker options are the same as in the second round, just with the updated SNAP model
- 6. Maker, Seqkit: Extract Maker's final gff and fasta transcripts using the tools packaged with Maker, combine the each scaffolds' annotations into a single file and remove duplicates with Seqkit
```bash
fasta_merge -d <maker_output_directory>/<round3>_master_datastore.log
seqkit rmdup round3.fasta > round3_deduplicated.fasta
```
- 7. BUSCO: Assess BUSCO completeness of final transcripts with BUSCO in `transcriptome` mode
- **Output:** gff files containing annotations for each sample, with predicted transcripts and proteins in fasta format.

## Variant calling and identification
- **Input:**
  - Cleaned reads (after trimming and filtering as in genome assembly)
  - Reference sequence GCF_000961545.1
- The variant calling steps are identical to that of Khalfan (2020)'s [implementation](https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/Khalfan) (an explanation can be found in the link), with covariate analysis disabled (due to an R bug) and updated from the DSL 1 (which is deprecated in the most recent version of Nextflow) to DSL 2
  - The short reads were aligned to the reference
- **Output:** vcf files for each sample
- Genes of interest were first identified by manually cross-referencing BUSCO output with the GCF_000961545.1 gff file
  - A script that automates the process using bcftools, when given a plain-text file of gene locations is included in the `bin` diretory.
    - It also counts the number of variants present within each range.
```bash
bcftools view <contig_name>:<start>-<stop> <indexed_vcf>
```
- The `bin` directory contains python scripts for extracting single-copy BUSCO gene sequences from BUSCO output and combining them across samples.
- **Note:** If you intend to predict variant effects with SnpEff and have used the GCF_000961545 reference for variant calling, downloading SnpEff's database of the reference GCF_000961545 won't work. Clarification and further instruction in the [extended](./info/README_Extended.md)

### BUSCO gene extraction
- **Input**
    - GCF_000961545 reference gff and fasta file
    - BUSCO output tables for each assembly
    - Cleaned reads
- **Output**
    * Multiple sequence alignments for all common single-copy BUSCO genes
    * KALLISTO tables describing the abundance of each single-copy gene
- 1. A script (`busco_to_gff.sh`) was written to map every single-copy BUSCO gene to a gene described in the gff file, obtaining a more precise range in the process.
    - The gene table was filtered to leave only the genes that were shared by every sample (3186 genes). A list of these shared genes was created by processing the output tables with a python script
* 2. Liftoff: Lift over genome annotations from GCF_000961545 reference gff onto scaffolds, generating a lifted gff file
```bash
liftoff <scaffolds> <reference_fasta> -g <reference_gff> -o <scaffolds>_lifted.gff
```
* 3. awk: Use the BUSCO-gff mapping from step 1. to filter the single-copy BUSCO genes from the lifted gffs, generating a tsv file with their new locations on the lifted gff*
    + *The original BUSCO output table describes the BUSCO gene locations specific for the GCF_000961545 reference. Their locations may be different with each assembly
* 4. gffread: Extract the sequences of the BUSCO genes in fasta format using the coordinates specified from the previous file. Every gene is stored in a separate fasta file
* 5. MAFFT: Combine the fasta files of the same genes for each sample into a single file and perform a multiple sequence alignment using MAFFT
```bash
mafft --preservecase --auto <combined_gene_fastas> > aligned.fasta
```
* 6. kallisto: Combine the per-sample fasta files, generating a set of per-sample gene sequences. Pass the cleaned reads for the sample as input to kallisto for transcript quantification
```bash
kallisto index -i index  <sample_combined_fasta> # Create a kallisto index (the "-i" flag specifies the index filename)
kallisto quant -i index -o <sample_quantified> <reads1> <reads2> # The "-o" flag specifies output directory name
```
* 7. R: Using the per-gene multiple sequence alignments from step 5, construct a phylogenetic tree with the neighbor joining method in R, then calculate the Generalized Robinson Foulds distance between the gene's tree and a reference tree


# References
- Assembly [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. GCF_000961545.1, S_schenckii_v1reference; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000961545.1/
- Assembly [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. GCF_000820605.1 S_brasiliensis_5110_v1; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000820605.1/
- Assembly [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. GCA_016097105.2, SVAR_v1.1; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/assembly/GCA_016097105.2/
- Assembly [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. GCA_001630445.1, Genome assembly ASM163044v1; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001630445.1/
- Assembly [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. GCA_001630435.1, Genome assembly ASM163043v1; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001630435.1/
- Assembly [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. GCA_021396225.1, Genome assembly ASM2139622v1; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_021396225.1/
- Assembly [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. GCA_021396235.1, Genome assembly ASM2139623v1; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_021396235.1/
- Assembly [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. GCA_016097115.2, Genome assembly SPRO_v1.1; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_016097115.2/
- Assembly [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. GCA_021396245.1, Genome assembly ASM2139624v1; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_021396245.1/
- Schlagenhauf, Edith and Wicker, Thomas (2016). The TREP platform: A curated database of transposable elements. Available at https://trep-db.uzh.ch.
- Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology, 34(5), Article 5. https://doi.org/10.1038/nbt.3519
- Bushnell, B., Rood, J., & Singer, E. (2017). BBMerge – Accurate paired shotgun read merging via overlap. PLOS ONE, 12(10), e0185056. https://doi.org/10.1371/journal.pone.0185056
- Campbell, M. S., Holt, C., Moore, B., & Yandell, M. (2014). Genome Annotation and Curation Using MAKER and MAKER-P. Current Protocols in Bioinformatics / Editoral Board, Andreas D. Baxevanis ... [et Al.], 48, 4.11.1-4.11.39. https://doi.org/10.1002/0471250953.bi0411s48
- Card, D. C., Adams, R. H., Schield, D. R., Perry, B. W., Corbin, A. B., Pasquesi, G. I. M., Row, K., Van Kleeck, M. J., Daza, J. M., Booth, W., Montgomery, C. E., Boback, S. M., & Castoe, T. A. (2019). Genomic Basis of Convergent Island Phenotypes in Boa Constrictors. Genome Biology and Evolution, 11(11), 3123–3143. https://doi.org/10.1093/gbe/evz226
- Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics (Oxford, England), 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560
- Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008
- Eilbeck, K., Moore, B., Holt, C., & Yandell, M. (2009). Quantitative measures for the management and comparison of annotated genomes. BMC Bioinformatics, 10(1), 67. https://doi.org/10.1186/1471-2105-10-67
- Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
- Flynn, J. M., Hubley, R., Goubert, C., Rosen, J., Clark, A. G., Feschotte, C., & Smit, A. F. (2020). RepeatModeler2 for automated genomic discovery of transposable element families. Proceedings of the National Academy of Sciences of the United States of America, 117(17), 9451–9457. https://doi.org/10.1073/pnas.1921046117
- GenBank [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. NC_015923.1; Sporothrix schenckii mitochondrion, complete genome; [cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/nucleotide/NC_015923.1
- Goubert, C., Craig, R. J., Bilat, A. F., Peona, V., Vogan, A. A., & Protasio, A. V. (2022). A beginner’s guide to manual curation of transposable elements. Mobile DNA, 13(1), 7. https://doi.org/10.1186/s13100-021-00259-7
- Holt, C., & Yandell, M. (2011). MAKER2: An annotation pipeline and genome-database management tool for second-generation genome projects. BMC Bioinformatics, 12(1), 491. https://doi.org/10.1186/1471-2105-12-491
- Katoh, K., Kuma, K., Toh, H., & Miyata, T. (2005). MAFFT version 5: Improvement in accuracy of multiple sequence alignment. Nucleic Acids Research, 33(2), 511–518. https://doi.org/10.1093/nar/gki198
- Khalfan, M. (2020, March 25). Variant Calling Pipeline using GATK4 – Genomics Core at NYU CGSB. https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
- Kolmogorov, M., Raney, B., Paten, B., & Pham, S. (2014). Ragout—A reference-assisted assembly tool for bacterial genomes. Bioinformatics, 30(12), i302–i309. https://doi.org/10.1093/bioinformatics/btu280
- Korf, I. (2004). Gene finding in novel genomes. BMC Bioinformatics, 5(1), 59. https://doi.org/10.1186/1471-2105-5-59
- Li, D., Liu, C.-M., Luo, R., Sadakane, K., & Lam, T.-W. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics (Oxford, England), 31(10), 1674–1676. https://doi.org/10.1093/bioinformatics/btv033
- Li, H. (2018). Minimap2: Pairwise alignment for nucleotide sequences. Bioinformatics (Oxford, England), 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191
- Lomsadze, A., Ter-Hovhannisyan, V., Chernoff, Y. O., & Borodovsky, M. (2005). Gene identification in novel eukaryotic genomes by self-training algorithm. Nucleic Acids Research, 33(20), 6494–6506. https://doi.org/10.1093/nar/gki937
- McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., & DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110
- Prjibelski, A., Antipov, D., Meleshko, D., Lapidus, A., & Korobeynikov, A. (2020). Using SPAdes De Novo Assembler. Current Protocols in Bioinformatics, 70(1), e102. https://doi.org/10.1002/cpbi.102
- Seppey, M., Manni, M., & Zdobnov, E. M. (2019). BUSCO: Assessing Genome Assembly and Annotation Completeness. Methods in Molecular Biology (Clifton, N.J.), 1962, 227–245. https://doi.org/10.1007/978-1-4939-9173-0_14
- Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), e0163962. https://doi.org/10.1371/journal.pone.0163962
- Smith, M.R. (2020a). Information theoretic Generalized Robinson-Foulds metrics for comparing phylogenetic trees. Bioinformatics 36: 5007–5013. doi: 10.1093/bioinformatics/btaa614
- SRA [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. SRX6367452; GSM3905740: Sample 1: hyphal form of S. schenckii; Sporothrix schenckii; RNA-Seq ;[cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/sra/SRX6367452
- Stanke, M., & Morgenstern, B. (2005). AUGUSTUS: A web server for gene prediction in eukaryotes that allows user-defined constraints. Nucleic Acids Research, 33(suppl_2), W465–W467. https://doi.org/10.1093/nar/gki458
- Ter-Hovhannisyan, V., Lomsadze, A., Chernoff, Y. O., & Borodovsky, M. (2008). Gene prediction in novel fungal genomes using an ab initio algorithm with unsupervised training. Genome Research, 18(12), 1979–1990. https://doi.org/10.1101/gr.081612.108
- The UniProt Consortium. (2023). UniProt: The Universal Protein Knowledgebase in 2023. Nucleic Acids Research, 51(D1), D523–D531. https://doi.org/10.1093/nar/gkac1052
- Wingett, S. W., & Andrews, S. (2018). FastQ Screen: A tool for multi-genome mapping and quality control. F1000Research, 7, 1338. https://doi.org/10.12688/f1000research.15931.2
- Yandell, M., & Ence, D. (2012). A beginner’s guide to eukaryotic genome annotation. Nature Reviews Genetics, 13(5), 329–342. https://doi.org/10.1038/nrg3174
