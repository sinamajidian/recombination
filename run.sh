(base) [smajidi1@login v11c]$ cat runner_v1a.run
#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 50G
#SBATCH --time 01:00:00


./runner_v1a.sh  t1


(base) [smajidi1@login v11c]$ cat runner_v1a.sh 
#!/bin/bash


project_folder=$1
mkdir ${project_folder}; cd ${project_folder}

snp_rate=0.001   # $2
depth=15 # $3

ref="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/archive/22.fa"
parameter_file="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/archive/parameter_file_SV"


num_thread=8


###### genome simulation mother father child #######
####################################################

mkdir genome_mother; cd genome_mother

SURVIVOR simSV ${ref} ${parameter_file} ${snp_rate} 0 sim

line_hap2=$(grep -n ">" sim.fasta | sed -n 2p  | tr : "\t" | cut -f 1)
sed -n 1,$((${line_hap2}-1))p sim.fasta  > sim1.fasta
sed -n ${line_hap2},'$p' sim.fasta  > sim2.fasta
grep -n ">" sim1.fasta
grep -n ">" sim2.fasta
sed -n 1,310000p sim1.fasta > sim1_swap.fa
sed -n '310001,$p' sim2.fasta >> sim1_swap.fa
sed -n 1,310000p sim2.fasta > sim2_swap.fa
sed -n '310001,$p' sim1.fasta >> sim2_swap.fa
sed 's/>22/>sim_Mh1/g' sim1_swap.fa > sim_Mh1.fa 
sed 's/>22/>sim_Mh2/g' sim2_swap.fa > sim_Mh2.fa



cd ..; mkdir genome_father; cd genome_father
SURVIVOR simSV ${ref} ${parameter_file} ${snp_rate} 0 sim

line_hap2=$(grep -n ">" sim.fasta | sed -n 2p  | tr : "\t" | cut -f 1)
sed -n 1,$((${line_hap2}-1))p sim.fasta  > sim1.fasta
sed -n ${line_hap2},'$p' sim.fasta  > sim2.fasta
grep -n ">" sim1.fasta
grep -n ">" sim2.fasta
sed -n 1,340000p sim1.fasta > sim1_swap.fa
sed -n '340001,$p' sim2.fasta >> sim1_swap.fa
sed -n 1,340000p sim2.fasta > sim2_swap.fa
sed -n '340001,$p' sim1.fasta >> sim2_swap.fa
sed 's/>22/>sim_Fh1/g' sim1_swap.fa > sim_Fh1.fa 
sed 's/>22/>sim_Fh2/g' sim2_swap.fa > sim_Fh2.fa

cd ..; mkdir genome_child; cd genome_child
cp ../genome_mother/sim_Mh1.fa .
cp ../genome_father/sim_Fh1.fa .



##########  read simulation child  mother father  #########
###########################################################
sample_fastq="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/archive/hifi_sample.fastq"

cd ..; mkdir reads_child;cd reads_child
pbsim --depth ${depth}  --sample-fastq  ${sample_fastq}  ../genome_child/sim_Mh1.fa  --prefix  sim_Mh1  &
sleep 5
pbsim --depth ${depth}  --sample-fastq  ${sample_fastq}  ../genome_child/sim_Fh1.fa  --prefix  sim_Fh1 
cat sim_Mh1_0001.fastq | sed '/^@S1_/ s/@S1_/@sim_Mh1_/' > sim_Mh1.fastq
cat sim_Fh1_0001.fastq | sed '/^@S1_/ s/@S1_/@sim_Fh1_/' > sim_Fh1.fastq
cat sim_Mh1.fastq  sim_Fh1.fastq  > reads.fastq
rm sim_Fh1_0001.fastq  sim_Fh1_0001.maf  sim_Fh1_0001.ref  sim_Fh1.fastq  sim_Mh1_0001.fastq  sim_Mh1_0001.maf  sim_Mh1_0001.ref  sim_Mh1.fastq

cd ..; mkdir reads_mother;cd reads_mother
pbsim --depth ${depth}  --sample-fastq  ${sample_fastq}  ../genome_mother/sim1.fasta  --prefix  sim_Mh1 &
sleep 5
pbsim --depth ${depth}  --sample-fastq  ${sample_fastq}  ../genome_mother/sim2.fasta  --prefix  sim_Mh2  
cat sim_Mh1_0001.fastq | sed '/^@S1_/ s/@S1_/@sim_Mh1_/' > sim_Mh1.fastq
cat sim_Mh2_0001.fastq | sed '/^@S1_/ s/@S1_/@sim_Mh2_/' > sim_Mh2.fastq
cat sim_Mh1.fastq  sim_Mh2.fastq  > reads_M.fastq
rm sim_Mh1_0001.fastq  sim_Mh1_0001.maf  sim_Mh1_0001.ref  sim_Mh1.fastq  sim_Mh2_0001.fastq  sim_Mh2_0001.maf  sim_Mh2_0001.ref  sim_Mh2.fastq

cd ..; mkdir reads_father;cd reads_father
pbsim --depth ${depth}  --sample-fastq  ${sample_fastq}  ../genome_father/sim1.fasta  --prefix  sim_Fh1  &
sleep 5
pbsim --depth ${depth}  --sample-fastq  ${sample_fastq}  ../genome_father/sim2.fasta  --prefix  sim_Fh2 
cat sim_Fh1_0001.fastq | sed '/^@S1_/ s/@S1_/@sim_Fh1_/' > sim_Fh1.fastq
cat sim_Fh2_0001.fastq | sed '/^@S1_/ s/@S1_/@sim_Fh2_/' > sim_Fh2.fastq
cat sim_Fh1.fastq  sim_Fh2.fastq  > reads_F.fastq
rm sim_Fh1_0001.fastq  sim_Fh1_0001.maf  sim_Fh1_0001.ref  sim_Fh1.fastq  sim_Fh2_0001.fastq  sim_Fh2_0001.maf  sim_Fh2_0001.ref  sim_Fh2.fastq



#########  mapping phasing of mother father ######
##################################################

cd ..; mkdir var_mother; cd var_mother
minimap2 -ax map-pb -t ${num_thread} ${ref} ../reads_mother/reads_M.fastq > reads.sam
samtools view -hb  -@ ${num_thread}  -q 1  reads.sam > reads_q1.bam
samtools sort      -@ ${num_thread} reads_q1.bam > reads_q1_sorted.bam
samtools index     -@ ${num_thread}  reads_q1_sorted.bam
rm  reads.sam reads_q1.bam

longshot  --bam  reads_q1_sorted.bam   --ref ${ref}  --out longshot.vcf
bgzip -c longshot.vcf  > longshot.vcf.gz
tabix longshot.vcf.gz
bcftools index longshot.vcf.gz
bcftools annotate -x INFO,^FORMAT/GT  longshot.vcf.gz  > longshot_ed.vcf

whatshap phase  --ignore-read-groups longshot_ed.vcf reads_q1_sorted.bam  --output phased.vcf --reference ${ref}  &



cd ..; mkdir var_father; cd var_father
minimap2 -ax map-pb -t ${num_thread} ${ref} ../reads_father/reads_F.fastq > reads.sam
samtools view -hb -@ ${num_thread} -q 1  reads.sam > reads_q1.bam
samtools sort     -@ ${num_thread} reads_q1.bam > reads_q1_sorted.bam
samtools index    -@ ${num_thread} reads_q1_sorted.bam
rm  reads.sam reads_q1.bam

longshot  --bam  reads_q1_sorted.bam   --ref ${ref}  --out longshot.vcf
bgzip -c longshot.vcf  > longshot.vcf.gz
tabix longshot.vcf.gz
bcftools index longshot.vcf.gz
bcftools annotate -x INFO,^FORMAT/GT  longshot.vcf.gz  > longshot_ed.vcf

whatshap phase  --ignore-read-groups longshot_ed.vcf reads_q1_sorted.bam  --output phased.vcf --reference ${ref}








#########  compare phasing ##########
#####################################

cd  ../genome_mother
grep "#" sim.vcf  > sim_e.vcf
grep -v "#"  sim.vcf|  sort -k1,1V -k2,2n >> sim_e.vcf

cd  ../genome_father
grep "#" sim.vcf  > sim_e.vcf
grep -v "#"  sim.vcf|  sort -k1,1V -k2,2n >> sim_e.vcf

cd ..
whatshap  compare var_mother/phased.vcf genome_mother/sim_e.vcf > phasing_compare_mother
whatshap  compare var_father/phased.vcf genome_father/sim_e.vcf > phasing_compare_father






#########  chopping ##########
##############################

mkdir chop; cd chop
cp ../var_mother/phased.vcf  phased_M.vcf
cp ../var_father/phased.vcf  phased_F.vcf

bgzip -c  phased_M.vcf > sim_M_ed.vcf.gz
tabix sim_M_ed.vcf.gz
bcftools index sim_M_ed.vcf.gz
bcftools consensus -f ${ref} --haplotype 1pIu sim_M_ed.vcf.gz > haplotype_M1_raw.fa
bcftools consensus -f ${ref} --haplotype 2pIu sim_M_ed.vcf.gz > haplotype_M2_raw.fa
# diff haplotype_M1_raw.fa haplotype_M2_raw.fa | wc -l
sed 's/>22/>haplotype_M1/g' haplotype_M1_raw.fa > haplotype_M1.fa 
sed 's/>22/>haplotype_M2/g' haplotype_M2_raw.fa > haplotype_M2.fa


bgzip -c  phased_F.vcf > sim_F_ed.vcf.gz
tabix sim_F_ed.vcf.gz
bcftools index sim_F_ed.vcf.gz
bcftools consensus -f ${ref} --haplotype 1pIu sim_F_ed.vcf.gz > haplotype_F1_raw.fa
bcftools consensus -f ${ref} --haplotype 2pIu sim_F_ed.vcf.gz > haplotype_F2_raw.fa
# diff haplotype_F1_raw.fa haplotype_F2_raw.fa | wc -l
sed 's/>22/>haplotype_F1/g' haplotype_F1_raw.fa > haplotype_F1.fa 
sed 's/>22/>haplotype_F2/g' haplotype_F2_raw.fa > haplotype_F2.fa

cat haplotype_M1.fa haplotype_M2.fa haplotype_F1.fa haplotype_F2.fa> haplotypes.fa

samtools faidx haplotypes.fa 
cat haplotypes.fa.fai


python3 /work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/archive/p4_chop_v1d.py  ../reads_child/reads.fastq reads_chopped.fastq 3000

minimap2 -ax map-pb -t ${num_thread}  haplotypes.fa reads_chopped.fastq > reads_chopped.sam
samtools view -h    -@ ${num_thread}  -q 1  reads_chopped.sam > reads_chopped_q1.sam

cut -f 1-5  reads_chopped_q1.sam > reads_chopped_q1_5cols.sam



