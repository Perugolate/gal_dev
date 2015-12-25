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
cd /scratch/perugolate/gal_anno/
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
cd /scratch/perugolate/gal_anno/
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
cd /scratch/perugolate/gal_anno/
blastx -query min1renormdiag.fasta -db uniprot_uniref90.trinotate.pep -num_threads 12 -max_target_seqs 1 -outfmt 6 > blastx.uniref.outfmt6
```

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
cd /scratch/perugolate/gal_anno/
blastx -query min1renormdiag.fasta -db /scratch/perugolate/db2/refseq_protein.00 -num_threads 12 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send stitle staxids sscinames evalue" -out Trinity.fasta.refseq_protein.00.outfmt6
```

```sh
for i in trin.[0-2][0-9]*; do sed "s/refseq_protein.00/refseq_protein.$i/g" trin.00.sh > trin.$i.sh; done
for i in trin.[0-2][0-9]*; do sbatch $i; done
```

# *Bombyx mori* orthologs

## setup

```sh
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-30/fasta/bombyx_mori/pep/Bombyx_mori.GCA_000151625.1.30.pep.all.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-30/fasta/bombyx_mori/ncrna/Bombyx_mori.GCA_000151625.1.30.ncrna.fa.gz
```
