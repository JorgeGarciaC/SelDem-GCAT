#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys

input_file=sys.argv[1]
out_file=sys.argv[2]

# In[84]:


data=[]
#for line in open('./clean_sn.txt'):
for line in open(input_file):
    line_data = line.split()
    data.append(line_data)


# In[85]:


import collections

a_dict = collections.defaultdict(list)

for i in data:
    a_dict[i[1]].append(i[0])
    


# In[114]:


import csv

w = csv.writer(open(out_file, "w"))
for key, val in a_dict.items():
    w.writerow([key] + val)

