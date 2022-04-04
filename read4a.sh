#!/bin/bash


chrom=$1
folder_name=$2
mkdir ${folder_name}
cd ${folder_name}

pos=$3
threads=$4
region_size=$5
parent=$6

echo "read_eval_code  ${chrom} ${folder_name} ${pos} ${threads} ${region_size}  ${parent}"
#parent="HG004"


HG002_bam_chrom="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/real_phaseRe/phasing/HG002/${chrom}/HG002_sorted_${chrom}.bam"
parent_bam_chrom="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/real_phaseRe/phasing/${parent}/${chrom}/${parent}_sorted_${chrom}.bam"
parent_vcf_chrom="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/real_phaseRe/phasing/${parent}/${chrom}/phased_${chrom}.vcf"
ref_chrom="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/archive/ref/${chrom}.fa"


# child's read
samtools view  --threads ${threads} -h -q 1 ${HG002_bam_chrom} ${chrom}:$((${pos}-${region_size}))-$((${pos}+${region_size})) > HG002_${chrom}.sam
samtools fastq  --threads ${threads} HG002_${chrom}.sam > HG002_${chrom}.fastq
rm HG002_${chrom}.sam


# parent's read
samtools view -h -q 1 ${parent_bam_chrom} ${chrom}:$((${pos}-${region_size}))-$((${pos}+${region_size})) | samtools view -hb > ${parent}_${chrom}.bam
samtools index ${parent}_${chrom}.bam

bgzip -c ${parent_vcf_chrom} > phased_${chrom}.vcf.gz
tabix phased_${chrom}.vcf.gz


echo "whatshap haplotag started"
whatshap haplotag --ignore-read-groups -o ${parent}_${chrom}_haplotagged.bam --reference ${ref_chrom} phased_${chrom}.vcf.gz  ${parent}_${chrom}.bam
samtools view -h ${parent}_${chrom}_haplotagged.bam > ${parent}_${chrom}_haplotagged.sam
samtools view ${parent}_${chrom}_haplotagged.bam > ${parent}_${chrom}_haplotagged_no_header.sam
grep "HP:i:1" ${parent}_${chrom}_haplotagged_no_header.sam > ${parent}_${chrom}_HP1.sam
grep "HP:i:2" ${parent}_${chrom}_haplotagged_no_header.sam > ${parent}_${chrom}_HP2.sam

rm ${parent}_${chrom}_haplotagged_no_header.sam
rm ${parent}_${chrom}_haplotagged.sam
rm phased_${chrom}.vcf*




echo "*********** loop started ** numebr of batches is" $(cat ${parent}_${chrom}_HP1.sam|wc -l) "*" $(cat ${parent}_${chrom}_HP2.sam|wc -l)

i=0
for i1 in $(seq 1 $(cat ${parent}_${chrom}_HP1.sam|wc -l)) ; do 
for i2 in $(seq 1 $(cat ${parent}_${chrom}_HP2.sam|wc -l)) ; do 

i=$((${i}+1))
echo ${i}
echo -n ">"  > ${parent}_${chrom}_HP1_HP2_${i}.fa
sed -n ${i1}p ${parent}_${chrom}_HP1.sam | cut -f 1,10 |tr '\t' '\n'  >> ${parent}_${chrom}_HP1_HP2_${i}.fa
echo -n ">"  >> ${parent}_${chrom}_HP1_HP2_${i}.fa
sed -n ${i2}p ${parent}_${chrom}_HP2.sam | cut -f 1,10 |tr '\t' '\n'  >> ${parent}_${chrom}_HP1_HP2_${i}.fa
samtools faidx  ${parent}_${chrom}_HP1_HP2_${i}.fa


minimap2 -ax map-pb -t ${threads} ${parent}_${chrom}_HP1_HP2_${i}.fa  HG002_${chrom}.fastq  > HG002_${parent}_${chrom}_${i}.sam

cut -f 1-6 HG002_${parent}_${chrom}_${i}.sam  > HG002_${parent}_${chrom}_${i}_6part.sam

done
done


echo "finished for ${parent}_${chrom}_${i} " 

