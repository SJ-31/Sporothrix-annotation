## Genome/transcriptome assembly
- **Input:**
  - Raw paired-end fastq files
  - Reference genome fasta and gff for *Sporothrix schenckii*, GCF_000961545.1
  - mtDNA reference
- 1. FastQC, MultiQC: Evaluate raw paired-end fastq files
    * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides information about read quality, overrepresented sequences, adapter content and other statistics on a per-fastq-file basis; MultiQC collates the information from multiple FastQC report files for convenience.
      - The quality evaluation here determines the corrections you should make to the raw files before assembly (the next step). Some common corrections include read trimming, as quality* drops heavily as the length increases, and adapter trimming.
        - *For [Illumina](https://www.illumina.com/Documents/products/technotes/technote_Q-Scores.pdf) devices, quality is measured with the Phred scoring system, and the higher the score for a given base, the lower the probability that the base was identified/called incorrectly. A score of 30 (Q30) is a commonly used threshold,  corresponding to a 99.9% base call accuracy.
    * FastQC produces an HTML report file providing either a green (passed), yellow (warning) or red (fail) flag at each assessment. While some red flags can be addressed (e.g. fails for read quality, adapter content), those that cannot (e.g. Per sequence GC content, Sequence Duplication levels) aren't necessarily indicative of sample preparation errors. See the [Galaxy Project's tutorial](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html#assess-quality-with-fastqc---short--long-reads) for more details.
     
- 2. FastP: Trim adapters and poly g (identified as a contaminant in prior step).
    * FastP is an all-purpose quality control tool that can perform adapter removal, read trimming and base correction
- 3. BBduk: Filter mtDNA
  - A kmer-based quality control tool, BBduk was used here for the function of filtering out reads that match strongly to a set of query sequences (mtDNA in this case). The [tool's documentation](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) has more usage examples
    - A kmer is simply a subsequence of DNA sequence that is *k* nucleotides long, and is a structure used in many bioinformatics tools. With BBduk, longer kmers correspond to higher stringency, as two sequences must share more bases to match. 
  - The kmer size setting for BBduk was set to the maximum of `k=31` to make filtering as stringent as possible. 
  - mtDNA filtering was added when in a previous iteration of the workflow, some of the output contigs were found to align to *Sporothrix* mtDNA. Initially, filtering was carried out post-assembly by using Minimap2 alignments to remove sequences that had mapped to the mtDNA. This broke up the contigs though, as they had mistakenly incorporated the mtDNA sequences into the assembly of valid chromosomes.
    * BBduk's filtering efficiency was validated by aligning the mtDNA with the filtered assembly with Minimap2 as before, and obtaining no alignments, indicating that all mtDNA sequences had been removed
- 4. Spades, Megahit: Assemble reads
  - Both ran in default settings.
- 5. BUSCO, Quast: Assess assembly quality
  - BUSCO ([Benchmarking Universal Single-Copy Orthologue](https://busco.ezlab.org/)) is a tool that assesses an assembly's completeness by checking for the presence of universal conserved genes in the assembly. The BUSCO genes can either be found complete and single-copy in an assembly, complete and duplicated, fragmented or absent. The significance of these results are context-dependent, see the [BUSCO user guide](https://busco.ezlab.org/busco_userguide.html#interpreting-the-results) for help with interpretation
    - Various sets of BUSCO genes (derived from the [OrthoDB](http://www.orthodb.org/) database) are available, but the most specific lineage to the study organism should be chosen to so that the presence of lineage-specific genes is tested.
  - [Quast](https://github.com/ablab/quast) meanwhile, assesses assembly quality by comparing the assembly against a reference. It provides measures of contiguity (e.g. N50 score) and similarity of assembly to the reference.
    - One way of interpreting contiguity is that for the same genome size, a more contiguous assembly is made up of less, but longer contigs. 
  - BUSCO was set to use the `sordariomycetes_odb10` lineage in `genome` mode, while Quast was provided with the GCF_000961545.1 reference genome
      - BUSCO needs to be set in `offline` mode or else it will download the given lineage database automatically in the directory of each run
      - Quast can assess several assemblies at once for comparison purposes
  - The Megahit assemblies produced more contiguous sets of contigs* and had higher BUSCO completeness, so these were selected for downstream analysis.
    - *A contig is a sequence that has been assembled by joining together overlapping short reads.
    - Transcriptome assembly was performed in much the same way as genome assembly, except mtDNA was not filtered, RNAspades was used as the assembler and finally BUSCO was set in `transcriptome` mode
- 6. Ragout:  Assemble contigs into scaffolds*
  - Sometimes the final contigs of an assembler are highly fragmented, made up of hundreds or even thousands of contigs; scaffolding is the process of linking these contigs together into longer units (scaffolds). With the assemblies in this study, scaffolding improved BUSCO completeness (up to ~5%, increasing the number of complete BUSCO genes) but also greatly sped up the annotation process.
  - Ragout was provided with the GCF_000961545.1 reference sequence to scaffold from, and run on default settings
      * The paths to files are given to Ragout with an `rcp` file, (an example can be found [here](./info/recipe.rcp)). A python script was used to automate writing this file from the command line
  - Several scaffolding tools - ntJoin, sspace, CAP3, RagTag and Ragout - were tested on the highly fragmented sample 2 Megahit contigs (7,639 sequences) before implementing for this step. Regardless of the flags used, Sspace and ntJoin showed no reduction in contig number; CAP3 reduced it by ~ 600 contigs. The reference-based RagTag showed promise, reducing down to  529 sequences. However, Ragout had the best performance (at the cost of speed), generating 13 scaffolds and was chosen.
      * Surprisingly, the addition of more reference genomes to Ragtag did not benefit the scaffolding and often increased the number of unplaced contigs
      * Different Ragout flags (`--refine`, `--solid-scaffolds` and `--repeats`), adding phylogenetic information (in the form of a Newick tree for the samples) and changing the synteny block program made no difference either.
- 7. Minimap2: Align scaffolds to reference sequence
- 8. Samtools: Extract scaffolds from sam file into separate fasta files for each chromosome
    - *Scaffolding and separating the results by chromosome was necessary to take advantage of nextflow's innate parallelization and reduce the time taken for genome annotation.
- **Output:** genome assemblies for each sample, split into fasta files for each chromosome present in the reference sequence
     
## Genome annotation
Yandell & Ence (2012)  provide a friendly introduction to the concepts and procedures involved in genome annotation and most of the notes here are based on their [paper](https://eaton-lab.org/slides/genomics/readings/Yandell-etal-2012.pdf).
- **Input:**
  - Genemarks hmm files for each sample, trained in-house with GenemarksES
  - Augustus species model trained with the Augustus web server* on the GCF_000961545.1 reference
    - Augustus and Genemarks are ab-initio gene predictors, using sequence characteristics (such as the presence of conserved gene motifs e.g. exon/intron boundaries) for gene prediction
    * *Training Augustus dynamically with the scaffolds of each sample would simply take too long, especially during the optimization step.
  - A database of repeat elements, generated by combining the TREP database with models produced in-house by RepeatModeler on 8 *Sporothrix* species (*S. brasiliensis*, *S. globosa* CBS strain, *S. globosa* strain SS01, *S. humicola*,
*S inflata*, *S. protearum*, *S. schenckii* and *S. variecibatus*)
    - With this information, genome sequences homologous to any repeat elements will be masked prior to annotation, which "hides" them from gene predictors. This will prevent them from being recognized as genes and producing false positive predictions
    * I did not plan on using 8 species for this step, but decided upon it after having received only 18 repeat element predictions from the *S. schenckii* genome. While this (and using the TREP database) may seem gratuitous, I believe lack of access to RepBase and the highly confounding effects that repeat elements have on ab initio gene predictions justify it (Yandell & Ence, 2012).
    - These repeat elements were also supplemented with repeat proteins provided with the Maker download
  - Data for evidence-driven prediction
    - Genome annotation tools also make use of existing sequence data (termed "evidence") to make alignment-based gene predictions
    - est evidence: includes the previously assembled transcriptome in addition to known *Sporothrix* genes collected from the NCBI nucleotide database with query `txid29907[Organism:exp] AND biomol_mrna[PROP]`
    - gff evidence from the GCF_000961545.1 reference
      - A generic feature format (gff) file stores information about genome features (e.g. genes, introns, exons), mapping each feature to a specific chromosome range. See [here](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) for details about the file structure
    - Protein evidence collected from UniProt Release 2023_02 with query `Sporothrix`
- 1. Maker: First round of Maker annotation, using est, gff, protein, repeat evidence and Genemarks model
- 2. SNAP: Train SNAP on the output from the first round of annotation.
  * SNAP is another ab initio gene predictor, and is used as part of the pipeline (rather than an input) because it is fast and can be trained using gff files (i.e. the maker output) as input (unlike Genemarks). Augustus is capable of this as well, but as alluded to earlier, it is prohibitively slow with the hardware that this pipeline was developed on.
- 3. Maker: Second round of annotation, disabling all evidence-based annotation and using Snap and Augustus for gene prediction only. This was set up so Maker would  populate its annotations with the gff file from the first round
    * Example of options [here](./info/maker_opts.ctl_Round2)
- 4. SNAP: Second round of SNAP training, using the Maker's second round output
    * The commands are the same as used in training SNAP in the first round
- 5. Maker: Final round of annotation with the most recent SNAP model
    * The maker options are the same as in the second round, just with the updated SNAP model
- 6. Maker, Seqkit: Extract Maker's final gff and protein + RNA transcripts (in fasta format) using the tools packaged with Maker, combine the each scaffolds' annotations into a single file and remove duplicates with Seqkit
- 7. BUSCO: Assess BUSCO completeness of final RNA transcripts with BUSCO in `transcriptome` mode
- The two rounds additional of Maker are only for retraining SNAP to improve its predictions, given that the Genemarks and Augustus models do not change. This three-round structure is the standard protocol to avoid overtraining the predictors and because additional rounds beyond this are usually redundant (Campbell et al., 2014).
- **Output:** gff files containing annotations for each sample, with predicted transcripts and proteins in fasta format.

## Variant calling and identification
- **Input:**
  - Cleaned reads (after trimming and filtering as in genome assembly)
  - Reference sequence GCF_000961545.1
- The variant calling steps are identical to that of Khalfan (2020)'s [implementation](https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/Khalfan) (an explanation can be found in the link), with covariate analysis disabled (due to an R bug) and updated from the DSL 1 (which is deprecated in the most recent version of Nextflow) to DSL 2
  - The short reads were aligned to the reference sequence
- **Output:** VCF files for each sample
- Genes of interest were first identified by manually cross-referencing BUSCO output with the GCF_000961545.1 gff file
  - Their regions were extracted from the VCF file using Bcftools
  - This was necessary because the BUSCO gene ranges specified in the output table are far larger than the actual sequence of the gene they contain.
- The `bin` directory contains python scripts for extracting single-copy BUSCO gene sequences from BUSCO output and combining them across samples

#### Building SnpEff database
- SnpEff is a tool that categorizes the variants in a VCF file and predicts their effects on gene function, taking as input a database constructed from the organism's reference genome. 
  - SnpEff's official downloadable *S. schenckii* dataset (with the GCF_000961545 reference) uses different chromosome headers (e.g. `Cont39`, `Cont41`) to that of the NCBI (e.g. `NW_015971139.1)`. 
    - Since the VCF file was generated using the latter convention, attempting to run SnpEff without formatting the headers first will not work
      - Even after doing this however, I still received ~70% `CHROMOSOME_NOT_FOUND` errors 
      - From looking at the SnpEff data directory, I suspect that not all chromosomes are present, or perhaps the same name refers to different chromosomes). 
- Fortunately, SnpEff does provide a way for users to build their own database. After rerunning SnpEff with this custom database, all errors disappeared
  - **Note:** The last time this was updated was July 7 2023, so this might not be necessary if SnpEff has updated their downloads
- **Steps**
  - Download the GenBank file for the GCF_000961545 reference
    - The easiest way to do this is with the NCBI's `datasets` CLI tool (downloading directly from the web page would always give me a corrupted zip file)
```bash
datasets download genome accession GCF_000961545.1 --include gbff
```
  - Update the SnpEff config file, located in the SnpEff installation directory with the new genome
```config
#---
# Databases are stored here
# E.g.: Information for 'hg19' is stored in data.dir/hg19/
# Custom GCF_000961545 
customGCF000961545.genome : Sporothrix schenckii # The format is <custom_genome_name>.genome : <species>
#
# You can use tilde ('~') as first character to refer to your home directory. 
# Also, a non-absolute path will be relative to config's file tar
# 
#---
data.dew = ./data/ # This is where you specify the data directory, in the SnpEff installation directory by default
```
 - Make a directory in the SnpEff data directory with the name used for the custom genome in the config, e.g. `mkdir customGCF000961545` in my example
 - Rename the genome gbff to `genes.gbk`
 - Run `java -jar <snpEff jar> build -genbank -v <custom_genome_name>`
 - You can now use the custom genome with snpEff by specifying `<custom_genome_name>` in any commands

### BUSCO gene extraction
The aim of this procedure was to locate each single-copy BUSCO gene in the genome, using them to construct a multiple sequence alignment to determine the phylogenetic distance between each sample. BUSCO's output table ("full_table.tsv") accurately specify the locations of each gene in the assembly, though the length specified refers to the length of the gene's CDS, so will be shorter than the specified range. For example, BUSCO reported that a 2581 bp range contained the Duf1741 gene, but the length was noted as 784 bp. 
- The steps below describe how this was confirmed in a separate procedure of mapping the NCBI and BUSCO genes, then lifting over the annotations. Doing it in this way is only necessary if you want the exact gene structure (which will be known once you find the gene it corresponds to on the gff file). If you just need the sequence in the range, use `seqkit -subseq <start>:<stop> <fasta_file>` 
  - A script (`busco_to_gff.sh`) was written for this process: after running BUSCO on the reference fasta file, the ranges in the BUSCO output tables are used to identify genome features in the  with the reference gff file. NCBI genes within the BUSCO range are recorded, then a liftover tool is used to find the new ranges of the NCBI genes (and thus the BUSCO genes) in the sample assemblies.
- **Input**
    - GCF_000961545 reference gff and fasta file
    - BUSCO output tables for each assembly
    - Mapping file describing which gene in the NCBI reference a given BUSCO gene corresponds to
      - Also, the gene table was filtered to leave only the genes that were shared by every sample (3186 common genes)
    - Cleaned reads
- **Output**
    * Multiple sequence alignments for all common single-copy BUSCO genes
    * KALLISTO tables describing the abundance of each single-copy gene
* 2. Liftoff: Lift over genome annotations from GCF_000961545 reference gff onto scaffolds, generating a lifted gff file
    - "Lifting over" annotations from a reference to an assembly means finding the locations of genome features specified by the reference in the assembly (if they are present), then annotating the assembly based on this correspondence. The Amidase gene (SPSK 01061) for example, was found on chromosome NW015971149 between bases 1105526-1107646 on the second *S. schenckii* sample, but is between bases 1091665-1093785 on the reference (for reference, BUSCO reported this gene to be )
* 3. awk: Use the BUSCO-gff mapping from step 1. to filter the single-copy BUSCO genes from the lifted gffs, generating a TSV file with their new locations on the lifted gff
    + *The original BUSCO output table describes the BUSCO gene locations specific for the GCF_000961545 reference. Their locations may be different with each assembly
* 4. gffread: Extract the sequences of the BUSCO genes in fasta format using the coordinates specified from the previous file. Every gene is stored in a separate fasta file
    - Initially, I tried obtaining the sequences of the common BUSCO genes for multiple sequence alignment by using the BUSCO transcripts from the Maker assessment (in the `<busco_output_folder>/<run_lineage>/busco_sequences/single_copy_busco_sequences` folder) directly, but the result was odd: ~40% of the genes had identical sequences between samples.
* 5. MAFFT: Combine the fasta files of the same genes for each sample into a single file and perform a multiple sequence alignment using MAFFT
    + MAFFT was chosen for multiple sequence alignment because of its speed and low memory use (~3200 genes were being aligned in parallel). An earlier version of this step used MUSCLE but I quickly ran into memory issues on my laptop and could not get the process to complete.
* 6. kallisto: Combine the per-sample fasta files, generating a set of per-sample gene sequences. Pass the cleaned reads for the sample as input to kallisto for transcript quantification
* 7. R: Using the per-gene multiple sequence alignments from step 5, construct a phylogenetic tree with the neighbor joining method in R, then calculate the Generalized Robinson Foulds distance between the gene's tree and a reference tree
- This step was implemented to identify which gene(s) was responsible for producing the variation that led to the observed reference tree, which was constructed after a run-through of this pipeline using all multiple sequence alignments from step 5.

## GO Annotation and enrichment analysis
The enrichment analysis here is to determine if any GO terms are associated with genes found to have significant variants
- 1. SnpSift: Filter VCF file to keep "high" effect variants and to remove variants with warnings
  + Variant effects in SnpEff are described using the [Sequence Ontology](http://www.sequenceontology.org/browser/obob.cgi) terminology, and the impact categories are defined by the SnpEff developers themselves
  + The exact documentation for SnpEff's VCF files can be found [here](http://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf)
  + "High" effects refers to variants that cause a significant degree of sequence divergence, including copy number variants, exon losses, frameshifts, rare amino acids, stop/start codon variants, splice site variants, transcript ablation.
- 2. R: Add GO terms into the annotation field of the VCF file
  + The mapping file generated by the Python script obtains GO annotations for the gene by querying UniProt for the gene (necessary because the NCBI doesn't store associated GO annotations in the Gene database). 
  + Genes with uncharacterized protein products are unlikely to have any GO annotations
    + 38% of *S. schenckii* genes in the VCF file have no annotations
- 3. R: Format list of genes and their GO annotations for enrichment analysis
  + This is simply a script to prepare the input files for Ontologizer, which filters out genes from the samples' variant stats table file
- 4. Ontologizer: Perform enrichment analysis on the gene list
  + The gene population was defined as all genes that had variants and GO annotations
  + The study set was defined was genes in the upper quartile of "high" variants
  + Ontologizer was run on the samples using the parent-child variant of Fisher's exact test, which, unlike Fisher's, is able to take into account the relational structure of the GO ontology (Bauer et al., 2008). 
    * Since a gene annotated to a given term *t* is also annotated to the parent terms of *t*, the probability that *t* is enriched is increased if one of *t*'s parents is enriched. This interaction between parent and sibling terms violates the assumption of independence in Fisher's test.
  + Given the high number of genes being tested in my sample, (6314 in the population, 1524 in the study set), I added the flag for multiple testing correction in Ontologizer (using Bonferroni correction). The unadjusted p-values will appear in the column besides the corrected ones.

# References
* Ashburner et al. Gene ontology: tool for the unification of biology. Nat Genet. 2000 May;25(1):25-9. DOI: 10.1038/75556 [abstract | full text]
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
* Sebastian Bauer and others, Ontologizer 2.0—a multifunctional tool for GO term enrichment analysis and data exploration, Bioinformatics, Volume 24, Issue 14, July 2008, Pages 1650–1651, https://doi.org/10.1093/bioinformatics/btn250
- Seppey, M., Manni, M., & Zdobnov, E. M. (2019). BUSCO: Assessing Genome Assembly and Annotation Completeness. Methods in Molecular Biology (Clifton, N.J.), 1962, 227–245. https://doi.org/10.1007/978-1-4939-9173-0_14
- Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), e0163962. https://doi.org/10.1371/journal.pone.0163962
- Smith, M.R. (2020a). Information theoretic Generalized Robinson-Foulds metrics for comparing phylogenetic trees. Bioinformatics 36: 5007–5013. doi: 10.1093/bioinformatics/btaa614
- SRA [Internet]. Bethesda (MD): National Library of Medicine (US), National Center for Biotechnology Information; [1988] – . Accession No. SRX6367452; GSM3905740: Sample 1: hyphal form of S. schenckii; Sporothrix schenckii; RNA-Seq ;[cited 2023 Jun 23]. Available from: https://www.ncbi.nlm.nih.gov/sra/SRX6367452
- Stanke, M., & Morgenstern, B. (2005). AUGUSTUS: A web server for gene prediction in eukaryotes that allows user-defined constraints. Nucleic Acids Research, 33(suppl_2), W465–W467. https://doi.org/10.1093/nar/gki458
- Ter-Hovhannisyan, V., Lomsadze, A., Chernoff, Y. O., & Borodovsky, M. (2008). Gene prediction in novel fungal genomes using an ab initio algorithm with unsupervised training. Genome Research, 18(12), 1979–1990. https://doi.org/10.1101/gr.081612.108
* The Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023. Genetics. 2023 May 4;224(1):iyad031. DOI: 10.1093/genetics/iyad031 [abstract | full text]
- The UniProt Consortium. (2023). UniProt: The Universal Protein Knowledgebase in 2023. Nucleic Acids Research, 51(D1), D523–D531. https://doi.org/10.1093/nar/gkac1052
- Wingett, S. W., & Andrews, S. (2018). FastQ Screen: A tool for multi-genome mapping and quality control. F1000Research, 7, 1338. https://doi.org/10.12688/f1000research.15931.2
- Yandell, M., & Ence, D. (2012). A beginner’s guide to eukaryotic genome annotation. Nature Reviews Genetics, 13(5), 329–342. https://doi.org/10.1038/nrg3174
