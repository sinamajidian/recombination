#!/usr/bin/env python
# coding: utf-8

# In[1]:





# In[13]:


file_fasta_address= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/v10a/t5/p7/haplotype_F2.fa"
file_fasta = open(file_fasta_address,'r');

# assimption one record only
record_name= []
sequence= ''

for line in file_fasta:
    line_strip = line.strip()
    if line_strip.startswith(">"):
        record_name.append(line_strip)
    else:
        sequence+=line_strip
print(len(sequence))


# In[14]:



pos= 11001096              # 1-based genomic position

pos_list =[11001096,11001260 ]

for pos in pos_list:

    base=sequence[pos-1]  # python 0-based position

    print("The base in position",pos,"is",base,".")


# In[ ]:





# In[ ]:




