#!/usr/bin/env bash
set -e

if [ $# -ne 2 ]
then
    echo "usage: $0 out_directory test_number

Runs scaffolders on test data made by the script scaff_test_make_test_cases.sh.
Should be run in the same directory where scaff_test_make_test_cases.sh
was run.
"
    exit
fi

outdir=$1
test_number=$2
echo "1"
contigs_prefix=$(readlink -e test.$test_number.contigs.fa)
echo "2"
reference=$(readlink -e test.$test_number.ref.fa)
reads_prefix=test.$test_number.reads
insert_min=350
insert_ave=500
insert_max=650
insert_sd=50
sopra_split=10000
genome_upper_bound=$(grep -v ">" $reference | wc | awk '{print $3-$1}')
reads_fwd=$(readlink -e $reads_prefix\_1.fq)
reads_rev=$(readlink -e $reads_prefix\_2.fq)
reads_shuffled=$(readlink -e $reads_prefix.fq)

mkdir $outdir
cd $outdir
contigs_abs_path=$contigs_prefix.tag.uniquely-tagged.contigs.fa
ln -s $contigs_abs_path contigs.fa
contigs=contigs.fa
samtools faidx $contigs
samtools faidx $reference

#____________________ MAPPING ___________________________
bowtie-build $contigs $contigs
bowtie -X $insert_max --sam $contigs -1 $reads_fwd -2 $reads_rev | samtools view -bS - > a.map.bowtie.raw.bam
samtools sort a.map.bowtie.raw.bam a.map.bowtie.sorted
scaffold-wrapper-mip-make-cov-file.pl a.map.bowtie.sorted.bam $contigs.fai a.map.bowtie.sorted.bam.mip.cov
samtools flagstat a.map.bowtie.sorted.bam > a.map.bowtie.sorted.bam.flagstat

bowtie2-build $contigs $contigs
bowtie2 -X $insert_max $contigs -1 $reads_fwd -2 $reads_rev | samtools view -bS - > a.map.bowtie2.raw.bam
samtools sort a.map.bowtie2.raw.bam a.map.bowtie2.sorted
samtools flagstat a.map.bowtie2.sorted.bam > a.map.bowtie2.sorted.bam.flagstat


# _____________________ RUN SCAFFOLDERS ________________
echo "-------------------- RUN ABySS -----------------------------------"
scaffold-wrapper-abyss.pl contigs.fa ABySS $reads_fwd $reads_rev

echo "-------------------- RUN Bambus2 ---------------------------------"
scaffold-wrapper-bambus2.pl $contigs Bambus2 a.map.bowtie2.raw.bam $insert_ave $insert_sd

echo "-------------------- RUN MIP -------------------------------------"
scaffold-wrapper-mip.pl -mapper_ops "-m 1 -v 0" $contigs bowtie $reads_fwd $reads_rev $genome_upper_bound $insert_min $insert_ave $insert_max a.map.bowtie.sorted.bam.mip.cov MIP

echo "----------------- RUN OPERA --------------------------------------"
scaffold-wrapper-opera.pl $contigs $reads_fwd $reads_rev OPERA bowtie 0

echo "----------------- RUN SCARPA -------------------------------------"
scaffold-wrapper-scarpa.pl -mapper_ops "-m 1 -v 0" $contigs $reads_shuffled SCARPA $insert_ave bowtie

echo "------------------- RUN SGA --------------------------------------"
# SGA uses first 20 contigs for some estimating. These test cases don't have
# that many, so need to change the option to sga-astat.py
number_of_contigs=$(grep -c ">" $contigs)
scaffold-wrapper-sga.pl -astat_ops "-b $number_of_contigs" $contigs a.map.bowtie2.raw.bam a.map.bowtie2.sorted.bam SGA

echo "--------------------- RUN SOAP2 ----------------------------------"
scaffold-wrapper-soap2.pl $contigs $reads_fwd $reads_rev SOAP2 $insert_ave

echo "--------------------- RUN SOPRA ----------------------------------"
scaffold-wrapper-sopra.pl -mapper_ops "-m 1 -v 0" $contigs bowtie $reads_fwd $reads_rev $sopra_split $insert_ave SOPRA

echo "-------------------- RUN SSPACE ----------------------------------"
scaffold-wrapper-sspace.pl $contigs $reads_fwd $reads_rev SSPACE $insert_ave 0.4


# _____________________ CHECK SCAFFOLDS ________________
for dir in ABySS Bambus2 MIP OPERA SCARPA SGA SOAP2 SOPRA SSPACE
do
    cd $dir
    scaff_test_check_using_tags.py $insert_ave $contigs_prefix.tag scaffolds.fa $reference.fai check_scaffolds
    scaff_test_check_using_tags_get_contig_layout.py check_scaffolds check_scaffolds.graph
    cd ..
    ln -s $dir/check_scaffolds.graph.pdf $dir.graph.pdf
done
