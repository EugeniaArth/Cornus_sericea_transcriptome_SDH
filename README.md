# Cornus_sericea_transcriptome_SDH
Transcriptome assembly and functional annotation of Cornus sericea, with a focus on identifying and characterizing shikimate dehydrogenase (SDH) enzymes.


**Sample preparation**
Samples of young and mature leaves, immature fruits, buds and flowers of Cornus sericea L. were collected in July 2024 in Moscow, Russia
Total RNA was extracted using RNeasy Mini Kit (QIAGEN, #74104) following manufacturer’s instructions. Frozen plant material was homogenized
with a mortar and a pestle in liquid nitrogen. After isolation total RNA was treated by turbo DNAse (TURBO DNA-free Kit, ThermoFisher Scientific, #AM1907). 
RNA concentration was quantified with Qubit RNA HS Assay Kit (ThermoFisher Scientific, #Q32852).  RNA integrity was assessed using RNA 6000 Pico Kit 
for Bioanalyzer (Agilent Technologies, #5067-1513) before and after turbo DNAse treatment.

*Short read (SR) sequencing*

A total 200.00 ng of total RNA was used for library preparation using TruSeq Stranded mRNA kit (Illumina). 
Samples were paired-end (PE) sequenced using NovaSeq 6000, 2 x 101 bp (Illumina Technologies, San Diego, USA) platform

*Long read (LR) sequencing*

For ONT libraries RNA was converted to cDNA using Mint cDNA synthesis kit (Evrogen, Moscow, Russia) with 20-25 cycles of amplification with the following modification: 
1) cDNA synthesis is primed with oligonucleotides that contain not only dT part and adapter part but also custom barcode, specific for each sample; 
2) optimal number of cycles for amplification was selected based on the results of qPCR with primers annealing at 3’ and 5’ ends. 
The number of cycles falling to the upper part of linear phase was considered to be optimal. Amplified cDNA (150 ng) was used as input for library preparation using 
standard protocol for genomic DNA with SQK-LSK114 kit. ONT read processing Basecalling was performed by Guppy 6.5.7 (Oxford Nanopore Technologies, Oxford, United Kingdom). 
After that, reads were processed by a custom pipeline NTproc (https://github.com/shelkmike/NTproc/, 8).
This work was performed using the core facilities of the LOPUKHIN FRCC PCM “Genomics, proteomics, metabolomics” 


**Data analysis**

*Quality check*

Quality check of the short-read (SR) and long-read (LR) reads was performed using fastqc tool for SR and NanoPlot for LR

```bash

fastqc -o fastqc_reports -t 6 *SR.fastq.gz

NanoPlot -t 10 --fastq LR.fastq.gz -o NanoPlot

```

Poly(A) tails were removed from ONT reads using Cutadapt v5.0 with the parameter “-a A{100}”.
Illumina reads were trimmed for adapters and poly(A) tails using fastp:

```bash
cutadapt -a A{100} -o Trimmed/LR.trimmed.fastq.gz LR.fastq.gz 

fastp -i SR_1.fastq.gz -I SR_trimmed_1.fastq.gz  -o SR_2.fastq.gz  -O SR_trimmed_2.fastq.gz  \
--trim_poly_g --detect_adapter_for_pe --trim_poly_x -h fastp.html -j fastp.json

```

*De novo assembly and assessing the assembly quality*

*De novo* assembly was performed using trans2express (). 
Before that, forward and reverse reads for Illumina  fastq files were concatenated forming two combined files 
for forward and reverse reads. For Nanopore samples, fastq files were also concatenated forming one combined file.


```bash

python trans2express.py  -1 concatenated_1.fastq.gz -2 concatenated_2.fastq.gz --long_reads concatenated_LR.fastq.gz \
-t 30 -m 120 -o Assembly

```

In the resulting Assembly/ folder it contains a subfolder, the most important final_assembly/ contain the files:
1) final.clust_cds_longest_iso.fasta
2) final.clust_proteins_longest_iso.fasta
3) final.clust_transcripts_longest_iso.fasta
4) final.clust_annotation_longest_iso.gff3


Initial assembly quality was evaluated using rnaQUAST using Cornus florida genome (NW_026775571.1) as a reference:

```bash

python rnaQUAST.py -c final.clust_transcripts_longest_iso.fasta -r reference/GCF_030987335.1_ASM3098733v1_genomic.fna -o rnaquast/

```

Two complementary approaches were employed to assess the completeness of the transcriptome assembly:

1) BUSCO (Benchmarking Universal Single-Copy Orthologs) analysis was used to evaluate completeness based on the presence of conserved orthologs, 
by comparing C. sericea transcripts with the Eukaryota (eukaryota_odb10), Eudicots (eudicots_odb10), and Viridiplantae (viridiplantae_odb10) datasets:

```bash

busco --offline -i final.clust_cds_longest_iso.fasta --lineage_dataset busco_db/eukaryota_odb10/ -c 30 -o busco_eukaryota \
--augustus --mode tran --download_path busco_db/

busco --offline -i final.clust_cds_longest_iso.fasta --lineage_dataset busco_db/eudicots_odb10/ -c 30 -o busco_eudicots \
--augustus --mode tran --download_path busco_db/

busco --offline -i final.clust_cds_longest_iso.fasta --lineage_dataset busco_db/viridiplantae_odb10/ -c 30 -o busco_viridiplantae \
--augustus --mode tran --download_path busco_db/

```

2) As a second approach, both SR and LR reads were mapped back to the assembled transcriptome using Minimap2, and the overall mapping rate 
was calculated. 

Mapping concatenated Illumina reads to assembly using Minimap2, covert .sam to .bam, sorting and indexing using samtools:

```bash

minimap2 -ax sr final.clust_transcripts_longest_iso.fasta concatenated_1.fastq.gz concatenated_2.fastq.gz -o total_SR_minimap2.sam -t 10

samtools view -bS total_SR_minimap2.sam | samtools sort -o total_SR_minimap2.bam

samtools index total_SR_minimap2.bam

```

Same for Nanopore reads:

```bash

minimap2 -ax map-ont final.clust_transcripts_longest_iso.fasta concatenated_LR.fastq.gz -o total_LR_minimap2.sam -t 10

samtools view -bS total_LR_minimap2.sam | samtools sort -o total_LR_minimap2.bam

samtools index total_LR_minimap2.bam

```

Get Mapping Statistics

```bash

samtools flagstat total_SR_minimap2.bam
samtools flagstat total_LR_minimap2.bam

```

The results were summarized in busco.csv file and visualized using code in Analysis.Rmd file  #1 Plotting BUSCO results 

At this step, final.clust_annotation_longest_iso.gff3 contains structural annotation of the assembly - with UTR, CDS, gene, 
exons, mRNA features. 
Therefore, on the next step I have performed functional annotation.


**Annotation**

First split fasta into batches:

```bash
seqkit split -s 15 final.clust_transcripts_longest_iso.fasta -O split_batches

```

Search for annotation sources was performed against different databases.

*Against UniProt:*

Download the UniProt Plant Protein Database

```bash
wget "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=taxonomy_id:33090 AND reviewed:true" \
  -O uniprot_plants_reviewed.fasta.gz

gunzip uniprot_plants_reviewed.fasta.gz

``` 

Create a Uniprot  Database

```bash
makeblastdb -in uniprot_plants_reviewed.fasta -dbtype prot -out uniprot_plants_db

```

Then use script Uniprot.blastx.best.sh - the best hits are saved as uniprot_best_hits


*Against Refseq:*

Download the Refseq Plant Protein Database

```bash
wget -r -nd -A "*.protein.faa.gz" ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/

gunzip *.protein.faa.gz      # Unzips all

cat *.protein.faa > all_refseq_proteins.faa
```

Create a Refseq  Database:

```bash
diamond makedb --in all_refseq_proteins.faa -d refseq_proteins

```

Then use script Refseq.blastx.best.sh  - the best hits are saved as refseq_best_hits.txt
 
*Against NCBI NR:*

Download the database
```bash
update_blastdb.pl --decompress nr

```

Create a database:
```bash
makeblastdb -in nr_viridiplantae.fasta -dbtype prot -out nr_viridiplantae_db
```
Confirm the database is functional:

```bash 
blastdbcmd -db nr -info

```

Move it in separate folder:

```bash 
mkdir viridiplantae_db
mv nr_viridiplantae* viridiplantae_db/
```

Use script NR_anno.sh  - the best hits are saved as NR_best_hits_combined.txt

*KEGG pathway annotation*

It was performed using https://www.genome.jp/kaas-bin/kaas_main
Protein sequencies were send with parameters:
Sp: ath, aly, gmx, gsj, fve, pop, jre, vvi, sly, nta, osa, zma, dct, rcn, pper, egr, brp, cit, tcc, qsu, oeu, bvg, dosa, ppp, peq, aof, atr, 
cre, mng, apro, olu, cme, gsl, ccp, psom, ini, peu, rcu


*Annotation against KOG/COG*

It was performed using EggNOG-Mapper
GO terms were parsed and extracted using python script GO_analysis_from_eggnogg_data.py

For GO terms, filtering for the highest GO level for each transcript was performed, to avoid exessive and unnesessary data.

- Addition of data to gff3 
 
 After the all source files were recieved, we have performed addition of data to basic gff3 file generated by trans2express pipeline using 
 python script gff3_annotation.py. 

 The priority of annotation was folowing: RefSeq > UniProt > NR

After annotation, the basic analysis and plots were created. The code used are present at Analysis.Rmd file
