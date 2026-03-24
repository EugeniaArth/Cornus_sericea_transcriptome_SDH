# Cornus_sericea_transcriptome_SDH
Transcriptome assembly and functional annotation of Cornus sericea, with a focus on identifying and characterizing shikimate dehydrogenase (SDH) enzymes.

**ANNOTATION**

First split fasta into batches:
```seqkit split -s 15 final.clust_transcripts_longest_iso.fasta -O split_batches```

- Against UniProt:

Download the UniProt Plant Protein Database
```wget "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=taxonomy_id:33090 AND reviewed:true" -O uniprot_plants_reviewed.fasta.gz```
```gunzip uniprot_plants_reviewed.fasta.gz``` 

Create a Uniprot  Database
`makeblastdb -in uniprot_plants_reviewed.fasta -dbtype prot -out uniprot_plants_db`

Then use script Uniprot.blastx.best.sh - the best hits are saved as uniprot_best_hits


Against Refseq:

Download the Refseq Plant Protein Database
wget -r -nd -A "*.protein.faa.gz" ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/
gunzip *.protein.faa.gz      # Unzips all
cat *.protein.faa > all_refseq_proteins.faa

Create a Refseq  Database:
diamond makedb --in all_refseq_proteins.faa -d refseq_proteins

Then use script Refseq.blastx.best.sh  - the best hits are saved as refseq_best_hits.txt
 
Against NR:
Download the database
update_blastdb.pl --decompress nr

Create a database:
makeblastdb -in nr_viridiplantae.fasta -dbtype prot -out nr_viridiplantae_db

Confirm the database is functional:
blastdbcmd -db nr -info

Move it inti separate folder:
mkdir viridiplantae_db
mv nr_viridiplantae* viridiplantae_db/

Then use script NR_anno.sh  - the best hits are saved as nr_best_hits.txt

KEGG pathway annotation was performed using https://www.genome.jp/kaas-bin/kaas_main
Protein sequencies were send with parameters:
Sp: ath, aly, gmx, gsj, fve, pop, jre, vvi, sly, nta, osa, zma, dct, rcn, pper, egr, brp, cit, tcc, qsu, oeu, bvg, dosa, ppp, peq, aof, atr, cre, mng, apro, olu, cme, gsl, ccp, psom, ini, peu, rcu


Annotation against KOG/COG was performed using EggNOG-Mapper
GO terms were parsed and extrcted using python script GO_analysis_from_eggnogg_data.py

