#!/usr/bin/env bash
set -e
seed=$1
fastaq make_random_contigs --seed $seed --name_by_letters 8 5000 all_contigs.fa
formatdb -p F -i all_contigs.fa
blastall -p blastn -i all_contigs.fa -d all_contigs.fa -m8 -e1e-10 -o all_contigs.fa.blastn
rm formatdb.log all_contigs.fa.???

insert=500
trim=40

pre=test.1
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 200 A B
ref A B
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre


pre=test.2
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 200 A B C
ref A B C
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre


pre=test.3
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 100 A B C
scaffold 100 A B D
ref A B C
ref D
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre


pre=test.4
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 180 A B C
scaffold 20 A B D
ref A B C
ref D
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre


pre=test.5
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 20 A B C
scaffold 180 A B D
ref A B C
ref D
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre


pre=test.6
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 200 A B C
scaffold 200 A D C
ref A B C
ref D
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre


pre=test.7
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 180 A B C
scaffold 20 A D C
ref A B C
ref D
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre


pre=test.8
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 20 A B C
scaffold 180 A D C
ref B
ref A D C
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre

pre=test.9
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 100 A B C D E F G
scaffold 100 A B C H E F G
ref A B C D E F G
ref H
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre


pre=test.10
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 180 A B C D E F G
scaffold 20 A B C H E F G
ref A B C D E F G
ref H
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre

pre=test.11
echo "contigs_file all_contigs.fa
outprefix $pre
scaffold 20 A B C D E F G
scaffold 180 A B C H E F G
ref A B C D E F G
ref H
trim $trim" > config.$pre
scaff_test_make_test_contigs.py config.$pre


