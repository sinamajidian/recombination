#!/usr/bin/env python

from sys import argv

vcf_file_address= argv[1]

vcf_file = open(vcf_file_address,'r')
first_first= True

vcf_file_out = open(vcf_file_address[:-4]+"_merg.vcf",'w')


for line in vcf_file:
    line_strip = line.strip()
    if line_strip.startswith('#'):
        if line_strip.startswith('#CHROM'):
            line_parts=line_strip.split('\t')
            vcf_file_out.write("\t".join(line_parts[:-1])+"\n")
        else:
            vcf_file_out.write(line_strip+"\n")

    else:
        line_parts=line_strip.split('\t')
        #chrom = line_parts[0]
        #var_pos = int(line_parts[1])    # genomic position of variants
        #ref=line_parts[3]
        #alt=line_parts[4]
        sim1=line_parts[9]
        sim2=line_parts[10]

        if sim1 == "0|1"   and sim2 == "./.":
            gt_both = "0|1"
        elif sim1 == "./." and sim2 == "0|1":
            gt_both = "1|0"
        elif sim1 == "0|1" and sim2 == "0|1":
            gt_both = "1|1"
            #print(var_pos)

        line_parts_ed=line_parts[:7]+[".","GT",gt_both]
        line_strip_ed='\t'.join(line_parts_ed)
        #print(line_strip_ed)
        vcf_file_out.write(line_strip_ed+"\n")

vcf_file_out.close()



