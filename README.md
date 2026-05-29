# Identification of 3-dehydroquinate dehydratase/shikimate dehydrogenase variants from *Cornus sericea* L. by high-throughput transcriptome sequencing

# Authors

Eugenia Nikonorova

All-Russian Scientific Research Institute of Medicinal and Aromatic Plants, Moscow, Russia
e-mail: gatiatulinaer@gmail.com

# Aim of the study:
Transcriptome assembly and functional annotation of *Cornus sericea*, with a focus on identifying and characterizing 
3-dehydroquinate dehydratase/shikimate dehydrogenases.


# Sample preparation
Samples of young and mature leaves, immature fruits, buds and flowers of *Cornus sericea* L. were collected in July 2024 in Moscow, Russia
Total RNA was extracted using RNeasy Mini Kit (QIAGEN, #74104) following manufacturer’s instructions. Frozen plant material was homogenized
with a mortar and a pestle in liquid nitrogen. After isolation total RNA was treated by turbo DNAse (TURBO DNA-free Kit, ThermoFisher Scientific, #AM1907). 
RNA concentration was quantified with Qubit RNA HS Assay Kit (ThermoFisher Scientific, #Q32852).  RNA integrity was assessed using RNA 6000 Pico Kit 
for Bioanalyzer (Agilent Technologies, #5067-1513) before and after turbo DNAse treatment.

**Short read (SR) sequencing**

A total 200.00 ng of total RNA was used for library preparation using TruSeq Stranded mRNA kit (Illumina). 
Samples were paired-end (PE) sequenced using NovaSeq 6000, 2 x 101 bp (Illumina Technologies, San Diego, USA) platform

**Long read (LR) sequencing**

For ONT libraries RNA was converted to cDNA using Mint cDNA synthesis kit (Evrogen, Moscow, Russia) with 20-25 cycles of amplification with the following modification: 
1) cDNA synthesis is primed with oligonucleotides that contain not only dT part and adapter part but also custom barcode, specific for each sample; 
2) optimal number of cycles for amplification was selected based on the results of qPCR with primers annealing at 3’ and 5’ ends. 
The number of cycles falling to the upper part of linear phase was considered to be optimal. Amplified cDNA (150 ng) was used as input for library preparation using 
standard protocol for genomic DNA with SQK-LSK114 kit. ONT read processing Basecalling was performed by Guppy 6.5.7 (Oxford Nanopore Technologies, Oxford, United Kingdom). 
After that, reads were processed by a custom pipeline NTproc (https://github.com/shelkmike/NTproc/, 8).
This work was performed using the core facilities of the LOPUKHIN FRCC PCM “Genomics, proteomics, metabolomics” 


# Data analysis

**Quality check**

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

**De novo assembly and assessing the assembly quality**

Before assembly,forward and reverse reads for Illumina  fastq files were concatenated forming two combined files 
for forward and reverse reads. For Nanopore samples, fastq files were also concatenated forming one combined file.

```bash

python trans2express.py  -1 concatenated_1.fastq.gz -2 concatenated_2.fastq.gz --long_reads concatenated_LR.fastq.gz \
-t 30 -m 120 -o Assembly

```

*De novo* assembly was performed using trans2express tool (https://github.com/albidgy/trans2express
). 

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

 BUSCO results were summarized in busco.csv file and visualized using code in  #1 Plotting BUSCO results chunk available at Analysis.Rmd file.

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

To get mapping statistics, run:

```bash

samtools flagstat total_SR_minimap2.bam
samtools flagstat total_LR_minimap2.bam

```

At this step, final.clust_annotation_longest_iso.gff3 contains structural annotation of the assembly - with UTR, CDS, gene, 
exons, mRNA features. 
Therefore, on the next step I have performed functional annotation.


# Annotation

First split fasta into batches - for annotation against NCBI NR and Refseq I have splitted fasta to multiple files - 1 per transcript:

```bash

seqkit split -s 1 final.clust_transcripts_longest_iso.fasta -O split_batches

```

Search for annotation sources was performed against different databases.

**Against UniProt:**

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


**Against Refseq:**

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

Then use script Refseq.blastx.best.sh  -  the best hits are saved as refseq_best_hits.txt
 
**Against NCBI NR:**

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

**KEGG pathway annotation**

It was performed using https://www.genome.jp/kaas-bin/kaas_main
Protein sequences were send with parameters:
Sp: ath, aly, gmx, gsj, fve, pop, jre, vvi, sly, nta, osa, zma, dct, rcn, pper, egr, brp, cit, tcc, qsu, oeu, bvg, dosa, ppp, peq, aof, atr, 
cre, mng, apro, olu, cme, gsl, ccp, psom, ini, peu, rcu

KEGG pathways were visualized using chunk #3 Analyse KEGG pathway and create a figure available at Analysis.Rmd file


**Annotation against KOG/COG**

It was performed using EggNOG-Mapper tool. After installation, we must download the required database files and then run it on 

```bash 

python /eggnog-mapper-2.1.12/download_eggnog_data.py

python /eggnog-mapper-2.1.12/emapper.py -i final.clust_proteins_longest_iso.fasta \
               -o  EggNOG --cpu 30

```

The resulted file is called  final.clust_proteins_longest_iso.emapper.annotations - and it will contain multiple lines  with COG, KOG and GO 
annotations. Because GO terms are always redundant, on the next step we parse GO terms and extract them using python script GO_analysis_from_eggnogg_data.py
For GO terms, filtering for the highest GO level for each transcript was performed, to avoid exessive and unnesessary data.

Data from GO annotation may be visualized using #2 Analyse GO annotation and create a figure chunk available at Analysis.Rmd file

**Addition of data to gff3**
 
 After the all source files were recieved, we have performed addition of data to basic gff3 file generated by trans2express pipeline 
 (final.clust_annotation_longest_iso.gff3) using  python script gff3_annotation.py. 
 
```bash 

python gff3_annotation.py

```

The priority of annotation was folowing: RefSeq > UniProt > NR with addition available GO terms.

After annotation, the basic analysis and plots were created. The code used are present at Analysis.Rmd file:
Summary of the annotation tools - chunk #4 Sum of annotation results in general, from resulted best_hits, etc files

Because the priority of annotation was RefSeq > UniProt > NR, not all data could be used in final gff3 file.
To count, we use the code  #5 Count annotions in final gff3 file. To plot the result we used chunk #6 Plot the counts of annotions in final gff3 file


# Extracting 3-dehydroquinate dehydratase/shikimate dehydrogenases 

At the beginning, I have search for DQD/SDHs presence:

```bash

grep -i "bifunctional 3-dehydroquinate dehydratase/shikimate dehydrogenase" annotated_final.gff3 \
> | cut -f3 \
> | sort \
> | uniq -c
   5 CDS
   5 exon
   5 five_prime_UTR
   5 gene
   5 mRNA
   5 three_prime_UTR
```

A total 5 mRNA with corresponding UTRs and CDSs.
Further, more detailed output with transcript name and coordinates was created:

```bash
grep -i "bifunctional 3-dehydroquinate dehydratase/shikimate dehydrogenase" annotated_final.gff3 \
| awk '
BEGIN{
    OFS="\t";
    print "ID","Feature","Start","End"
}
$3=="mRNA" || $3=="five_prime_UTR" || $3=="CDS" || $3=="three_prime_UTR" {
    print $1,$3,$4,$5
}'

ID	Feature	Start	End
NODE_10980_length_2779_cov_121.715385_g1243_i5	mRNA	1	2779
NODE_10980_length_2779_cov_121.715385_g1243_i5	five_prime_UTR	1	130
NODE_10980_length_2779_cov_121.715385_g1243_i5	CDS	131	1633
NODE_10980_length_2779_cov_121.715385_g1243_i5	three_prime_UTR	1634	2779
NODE_18517_length_2241_cov_255.196168_g8755_i0	mRNA	1	2241
NODE_18517_length_2241_cov_255.196168_g8755_i0	five_prime_UTR	2018	2241
NODE_18517_length_2241_cov_255.196168_g8755_i0	CDS	452	2017
NODE_18517_length_2241_cov_255.196168_g8755_i0	three_prime_UTR	1	451
NODE_20612_length_2133_cov_1535.833013_g8101_i1	mRNA	1	2133
NODE_20612_length_2133_cov_1535.833013_g8101_i1	five_prime_UTR	1894	2133
NODE_20612_length_2133_cov_1535.833013_g8101_i1	CDS	295	1893
NODE_20612_length_2133_cov_1535.833013_g8101_i1	three_prime_UTR	1	294
NODE_21648_length_2083_cov_1398.299410_g8244_i1	mRNA	1	2083
NODE_21648_length_2083_cov_1398.299410_g8244_i1	five_prime_UTR	1951	2083
NODE_21648_length_2083_cov_1398.299410_g8244_i1	CDS	367	1950
NODE_21648_length_2083_cov_1398.299410_g8244_i1	three_prime_UTR	1	366
NODE_29040_length_1767_cov_22.706636_g12409_i1	mRNA	1	1767
NODE_29040_length_1767_cov_22.706636_g12409_i1	five_prime_UTR	1	112
NODE_29040_length_1767_cov_22.706636_g12409_i1	CDS	113	1672
NODE_29040_length_1767_cov_22.706636_g12409_i1	three_prime_UTR	1673	1767

```

Then, the mRNA and amino acid sequences were extracted and saved. Example:

```bash

awk -v '/NODE_20612_length_2133_cov_1535.833013_g8101_i1/ -v RS='>' '$1 == seq {print RS $0}' final.clust_transcripts_longest_iso.fasta 

awk -v '/NODE_20612_length_2133_cov_1535.833013_g8101_i1/ -v RS='>' '$1 == seq {print RS $0}' final.clust_proteins_longest_iso.fasta 

```

For visualization in IGV viewer and search for alternative splicing SR and LR were mapped to transcriptome assembly:

```bash

minimap2 -ax sr  final.clust_transcripts_longest_iso.fasta mRNA_Cornus/Illumina/RNA_S11398Nr4.1.fastq.gz mRNA_Cornus/Illumina/RNA_S11398Nr4.2.fastq.gz -t 10 > berries.sam 
minimap2 -ax sr  final.clust_transcripts_longest_iso.fasta mRNA_Cornus/Illumina/RNA_S11398Nr6.1.fastq.gz mRNA_Cornus/Illumina/RNA_S11398Nr6.2.fastq.gz -t 10 > leaves.sam 
 minimap2 -ax sr final.clust_transcripts_longest_iso.fasta mRNA_Cornus/Illumina/RNA_S11398Nr7.1.fastq.gz mRNA_Cornus/Illumina/RNA_S11398Nr7.2.fastq.gz -t 10 > flowers.sam 


minimap2 -ax splice -uf -k14 final.clust_transcripts_longest_iso.fasta mRNA_Cornus/Nanopore/BC01.fastq.gz -t 10  > berries_LR.sam
minimap2 -ax splice -uf -k14 final.clust_transcripts_longest_iso.fasta mRNA_Cornus/Nanopore/BC08.fastq.gz -t 10  > flowers_LR.sam
minimap2 -ax splice -uf -k14 final.clust_transcripts_longest_iso.fasta mRNA_Cornus/Nanopore/BC05.fastq.gz -t 10  > leaves_LR.sam

```
The resulted .sam files were converted to .bam, sorted and indexed using Samtools:

```bash

samtools view -bS -F 0x4 berries.sam > berries.bam
samtools view -bS -F 0x4 flowers.sam > flowers.bam
samtools view -bS -F 0x4 leaves.sam > leaves.bam

samtools sort berries.bam > berries.sorted.bam
samtools sort leaves.bam > leaves.sorted.bam 
samtools sort flowers.bam > flowers.sorted.bam

samtools index -M  *sorted.bam
 
```
After that, every DQD/SDH variant were manually inspected. Due to the inconsistent nomenclature of DQD/SDH genes among plant species, 
identifiers (DQD/SDH1–5)  were assigned according to their occurrence in the annotation. One alternatively spliced isoform generated 
through intron retention was identified. Examples of IGV screens are saved in IGV folder

Therefore, the sequence from 1586 to 2031 bp was removed and translated. As the CDS	of the initial variant was located in 131-1633 bp,
the overlap inside CDS was 1586–1633 = 48 bp (16 amino acids removed near the C-terminus without frameshift)
Run script for that (for only first sequence in file)

And save result

```bash

python3 bin/splice_interval_to_fasta.py \
  -i Files/DQD_SDH/DQD_SDH_possible_nt.fasta \
  -o Files/DQD_SDH/DQD_SDH_possible_nt.spliced_1586-2031.fasta \
  --start 1586 --end 2031 --only-first

```

Next TransDecoder was used to check the CDS of the new protein sequences:

```bash

TransDecoder.LongOrfs -t Files/DQD_SDH/DQD_SDH_possible_nt.spliced_1586-2031.fasta

TransDecoder.Predict -t  Files/DQD_SDH/DQD_SDH_possible_nt.spliced_1586-2031.fasta --no_refine_starts

```

The TransDecoder ORF was classified as 5′ partial because start-site refinement failed due to insufficient training data. Therefore, the 3′ CDS boundary predicted by TransDecoder 
was retained, while the 5′ CDS boundary was manually curated based on the presence of an upstream ATG, correct reading frame, absence of premature stop codons, and similarity 
to homologous DQD/SDH proteins. Also, NCBI ORFfinder was used for additional check.
The script was used to extract 5’UTR/CDS/3’UTR from the spliced transcript, translate CDS to protein, and write all outputs plus QC checks to files. 

```bash
python bin/evaluate_utrs_cds_transdecoder.py \
  -i Files/DQD_SDH/DQD_SDH_possible_nt.spliced_1586-2031.fasta \
  --cds-start 131 --cds-end 1627 \
  -o Files/DQD_SDH/TransDecoder_eval \
  --trim-cds --only-first
  
```
The resulted protein sequence was named as DQD/SDH1.1 variant, while DQD/SDH1.2 was used for the variant with intron retention

On the next step, DQD/SDH proteins properties were analyzed.

# Analysis of DQD/SDH proteins properties

Stability of C. sericea DQD/SDHs predicted from amino acid sequence using ProtParam (https://web.expasy.org/protparam/) 

Code align DQD_SDH_possible_pr.fasta and DQD_SDH_possible_nt.fasta with Clustal Omega
Compute gap-ignored p-distance matrices from the alignments (ignores columns with - in either seq)
Build trees with IQ-TREE (uses iqtree2 if available, else iqtree)

Percent identity (%) for nucleotide  and amino acid sequences of C. sericea DQD/SDHs was calculated based on a distance matrixes without with 
gaps generated from Clustal Omega alignments using align_distance_iqtree.py script


```bash
python3 bin/align_distance_iqtree.py \
  --protein-fasta Files/DQD_SDH/DQD_SDH_possible_pr.fasta \
  --nt-fasta Files/DQD_SDH/DQD_SDH_possible_nt.fasta \
  --outdir Files/DQD_SDH/align_tree \
  --threads 8
  
```

The R code for percent identity (%) calculation is present in chunk #7 in Analysis.Rmd. Data are saved as 
DQD_SDH_possible_pr.dist_nogaps.tsv and DQD_SDH_possible_nt.dist_nogaps.tsv and visualized using tree data and chunk #8 in Analysis.Rmd

The multiple DQD/SDHs of dicot species protein sequence alignment was performed with Clustal Omega and tree was build using iqtree
with. parameters -m MFP -bb 1000 -nt 10

```bash

python3 bin/align_and_tree_one.py \
  -i Files/Alignment/For_alignment.fasta \
  -o Files/Alignment/ClustalO.aln.fa \
  --seq-type protein \
  --threads 8

```

Finally, resulting  phylogenetic tree was constructed for six C. sericea DQD/SDHs assessed in this study and the 44 sequences available in public 
databases was visualized using iTOL (https://itol.embl.de)

# Prediction of chloroplast transit peptide presence 


Subcellular localization of DQD/SDH variants was predicted using multiple tools, including WoLF PSORT, TargetP, DeepLoc 2.0, PredSL, and ChloroP 
priority was given to the predictions generated by TargetP-2.0 and PredSL, as these tools are more recent and are based on machine learning approaches,
but the combination of 2 or more positive results from the tools were also considered as positive cpTP signal. 

Next, we additionally analyzed the composition of the first 50 amino acid residues in cpTP-containing and non-cpTP DQD/SDHs of dicot species using sequence alignment.
sequences - Logo plot of the relative occurrence and bits (express the total height of each AA position corresponds to the conservation in that position) of the AA in the N-terminal 
was created using corresponding chunk #9A in Analysis.Rmd. The mean percent of uncharged, acidic, basic amino acids and proline were calculated using chunk #9B in Analysis.Rmd


# Phylogenetic study and the comparison of key amino acids in the 3-DHD and SDH domains

Data from multiple DQD/SDHs of dicot species protein sequence alignment and phylogenetic tree, which were received at the previous step were used for 
the comparison of key amino acids in the 3-DHD and SDH domains depends on clades. The code in chunk #10 was used for extraction of residue positions of interest based on
AtSDH sequence, sequence reordering, annotation of data based on cpTP presence, QDH and GA-forming activity and visualization based on phylogenetic tree data.


# Data Archiving Statement
Sequence data that support the findings of this study have been deposited in of NCBI at https://www.ncbi.nlm.nih.gov/. The associated BioProject is PRJNA1466783 
and SRA accession numbers are SRR38625913, SRR38625914, SRR38625915, SRR38625916, SRR38625917, SRR38625918


