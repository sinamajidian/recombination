#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[65]:


#!/usr/bin/python3

import numpy as np
from sys import argv


#file_fq_input_addrss= "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/archive/sample_22_hifi.fastq"
# file_fq_output_addrss = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/recombination/r1a/p4/"+"_chopped.fastq"

file_fq_input_addrss =  argv[1]
file_fq_output_addrss = argv[2]

file_fq_input= open(file_fq_input_addrss,'r');
file_fq_ouput= open(file_fq_output_addrss,'w');


chunk_length = int(argv[3]) # 3000



record = []
for line in file_fq_input:
    line_strip = line.strip()
    # if line_strip.startswith("@"):   # new recrod started 
    #    record = [read_name_new_record]
        
    if len(record)<4:
        record.append(line_strip)
    else:
        read_name,sequence,pluse_line,quality  = record        
        read_length = len(sequence)
        if read_length != len(quality): 
            print("ERROR, length of sequence !=quality ",len(sequence), len(quality))
        
        chunk_num = int(read_length/chunk_length)
        for chunk_i in range(chunk_num):
            if chunk_i < chunk_num-1 :
                sequence_i =  sequence[chunk_length*chunk_i:chunk_length*(chunk_i+1)]
                quality_i  = quality[chunk_length*chunk_i:chunk_length*(chunk_i+1)]
            else: 
                sequence_i =  sequence[chunk_length*chunk_i:]
                quality_i  = quality[chunk_length*chunk_i:]
            file_fq_ouput.write(read_name+"_part_"+str(chunk_i)+"\n")
            file_fq_ouput.write(sequence_i+"\n"+pluse_line+"\n"+quality_i+"\n")


        read_name_new_record=line_strip
        record = [read_name_new_record]
        
        
# last record 
read_name,sequence,pluse_line,quality  = record        
read_length = len(sequence)
chunk_num = int(read_length/chunk_length)
for chunk_i in range(chunk_num):
    if chunk_num  != chunk_i:
        sequence_i =  sequence[chunk_length*chunk_i:chunk_length*(chunk_i+1)]
        quality_i  = quality[chunk_length*chunk_i:chunk_length*(chunk_i+1)]
    else: 
        sequence_i =  sequence[chunk_length*chunk_i:]
        quality_i  = quality[chunk_length*chunk_i:]
    file_fq_ouput.write(read_name+"_part_"+str(chunk_i)+"\n")
    file_fq_ouput.write(sequence_i+"\n"+pluse_line+"\n"+quality_i+"\n")

file_fq_input.close()
file_fq_ouput.close()
            
print("The chopped version is in ",file_fq_output_addrss)
    



    


# In[ ]:




