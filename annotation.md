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

