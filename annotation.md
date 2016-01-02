# Trinotate annotation

## setup

```sh
# setup swissprot
wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/uniprot_sprot.trinotate_v2.0.pep.gz
mv uniprot_sprot.trinotate_v2.0.pep.gz uniprot_sprot.trinotate.pep.gz
gunzip uniprot_sprot.trinotate.pep.gz
makeblastdb -in uniprot_sprot.trinotate.pep -dbtype prot
# setup uniref
wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/uniprot_uniref90.trinotate_v2.0.pep.gz
mv uniprot_uniref90.trinotate_v2.0.pep.gz uniprot_uniref90.trinotate.pep.gz
gunzip uniprot_uniref90.trinotate.pep.gz
makeblastdb -in uniprot_uniref90.trinotate.pep -dbtype prot
# setup pfam
wget https://data.broadinstitute.org/Trinity/Trinotate_v2.0_RESOURCES/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
# predict ORFs
TransDecoder.LongOrfs -t min1renormdiag.fasta
```

## blast swissprot

```sh
#! /bin/bash
#SBATCH --job-name=uniprot_blast
#SBATCH --mail-type=all
#SBATCH --mail-user=
#SBATCH --mem=2048
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=8
cd /scratch/$USER/gal_anno/
# takes ~1.5 days
blastx -query min1renormdiag.fasta -db uniprot_sprot.trinotate.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.uniprot.outfmt6
```

## hmmer pfam

```sh
#! /bin/bash
#SBATCH --job-name=hmmscan
#SBATCH --mail-type=all
#SBATCH --mail-user=
#SBATCH --mem=2048
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=8
cd /scratch/$USER/gal_anno/
# takes ~15 h
hmmscan --cpu 8 --domtblout TrinotatePFAM.out Pfam-A.hmm longest_orfs.pep > pfam.log
```

## blast uniref

```sh
#! /bin/bash
#SBATCH --job-name=uniref_blast
#SBATCH --mail-type=all
#SBATCH --mail-user=
#SBATCH --mem=12288
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=12
cd /scratch/$USER/gal_anno/
blastx -query min1renormdiag.fasta -db uniprot_uniref90.trinotate.pep -num_threads 12 -max_target_seqs 1 -outfmt 6 > blastx.uniref.outfmt6
```

### Time outs

Timed out so chunked up remaining sequences from input fasta:
```sh
# chunk up remaining seqs according to last seq with a hit
sed -n $(grep $(tail blastx.uniref.outfmt6 -n1 | cut -f1) min1renormdiag.fasta -n | cut -f1 -d ":"),588978p min1renormdiag.fasta > min1renormdiag.sp2.fasta
# split into 20 chunks
pyfasta split -n 20 min1renormdiag.sp2.fasta
```

```sh
#! /bin/bash
#SBATCH --job-name=uniref_blast
#SBATCH --mail-type=all
#SBATCH --mail-user=
#SBATCH --mem=12288
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=12
cd /scratch/$USER/gal_anno/
blastx -query min1renormdiag.fasta -db uniprot_uniref90.trinotate.pep -num_threads 12 -max_target_seqs 1 -outfmt 6 > blastx.min1renormdiag.fasta.uniref.outfmt6
```

```sh
# make/submit a job for each chunk fasta
for i in *sp2.[0-9]*; do sed "s/min1renormdiag.fasta/$i/g" ur.sh > ur.$i.sh; done
for i in ur.min*sh; do sbatch $i; done
```

Which also timed out, so reran the remaining parts of each chunk:
```sh
# Chop each chunk according to last seq to retrieve a hit
for i in {01..11}; do 
  awk "/$(tail -n1 blastx.min1renormdiag.sp2.${i}.fasta.uniref.outfmt6 | cut -f1)/{y=1}y" min1renormdiag.sp2.${i}.fasta > min1renormdiag.sp3.${i}.fasta
done
# renaming in job scripts
for i in {01..11}; do 
  sed 's/sp2/sp3/g' ur.min1renormdiag.sp2.$i.fasta.sh > ur.min1renormdiag.sp3.$i.fasta.sh
done
# increase the time limit, since the all timed out
for i in ur.min1renormdiag.sp3.*.fasta.sh; do 
  sed -i 's/1-00:00:00/2-12:00:00/g' $i
done
# submit them all
for i in ur.min1renormdiag.sp3.*.fasta.sh; do sbatch $i; done
```

Note, the above only covers chunks 1-11 as there are running jobs - will need to do this for the rest of the running jobs (12-19) as they time out. Also ...sp2.00... failed almost immediately for unknown reasons. Re-submitted it with longer time.

# refseq annotation

## setup

```sh
update_blastdb.pl --decompress taxdb
update_blastdb.pl --decompress refseq_protein
```

## blastx refseq_protein

`trin.00.sh`:
```sh
#! /bin/bash
#SBATCH --job-name=refseq_protein.00
#SBATCH --mail-user=
#SBATCH --mail-type=all
#SBATCH --mem=3072
#SBATCH --cpus-per-task=12
#SBATCH --time=5-00:00:00
cd /scratch/$USER/gal_anno/
blastx -query min1renormdiag.fasta -db /scratch/$USER/db2/refseq_protein.00 -num_threads 12 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send stitle staxids sscinames evalue" -out Trinity.fasta.refseq_protein.00.outfmt6
```

```sh
for i in {01..21}; do sed "s/refseq_protein.00/refseq_protein.$i/g" trin.00.sh > trin.$i.sh; done
for i in trin.[0-2][0-9]*; do sbatch $i; done
```

Hit RAM limit so split near problem seq and re-run
```sh
sed -n 315411,588978p min1renormdiag.fasta > min1renormdiag.2.fasta
makeblastdb -in min1renormdiag.2.fasta -dbtype nucl
for i in {01..21}; do sed "s/refseq_protein.00/refseq_protein.$i/g" trin.00.2.sh > trin.$i.2.sh; done
for i in trin.[0-2][0-9]*.2.sh; do sbatch $i; done
```


# *Bombyx mori* orthologs

## setup

```sh
wget ftp://silkdb.org/pub/current/otherdata/Gene_ontology/silkworm_glean_gene.go
wget ftp://silkdb.org/pub/current/Gene/Glean_genes/silkworm_glean_pep.fa.tar.gz
tar xvf silkworm_glean_pep.fa.tar.gz
makeblastdb -in silkpep.fa -dbtype prot
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-30/fasta/bombyx_mori/ncrna/Bombyx_mori.GCA_000151625.1.30.ncrna.fa.gz
```

## blastx *B. mori* proteome

```sh
#! /bin/bash
#SBATCH --job-name=bmori_prot
#SBATCH --mail-user=
#SBATCH --mail-type=all
#SBATCH --mem=3072
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
cd /scratch/$USER/gal_anno/
blastx -query min1renormdiag.fasta -db /scratch/$USER/gal_anno/silkpep.fa -num_threads 8 -evalue 1e-3 -max_target_seqs 1 -outfmt 6 -out bmori_prot.outfmt6
```

```sh
#! /bin/bash
#SBATCH --job-name=tbm_prot
#SBATCH --mail-user=
#SBATCH --mail-type=all
#SBATCH --mem=3072
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
cd /scratch/$USER/gal_anno/
tblastn -query silkpep.fa -db /scratch/$USER/gal_anno/min1renormdiag.fasta -num_threads 8 -evalue 1e-3 -max_target_seqs 1 -outfmt 6 -out tbm_prot.outfmt6
```

```sh
crb-blast -q min1renormdiag.fasta -t silkpep.fa -e 1e-5 -h 12 -o annotation.tsv
crb-blast -q min1renormdiag.fasta -t Bombyx_mori.GCA_000151625.1.30.ncrna.fa -e 1e-5 -h 12 -o ncrna.tsv
```
