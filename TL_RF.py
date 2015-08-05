import dendropy
from dendropy.calculate import treecompare
from dendropy.utility.fileutils import find_files
import csv
import sys
import os
import pandas as pd
o_file = sys.argv[1]
#Path to first batch of files
i_file = sys.argv[2]
#Path to second batch of files
ext = sys.argv[3]
#Expected extension of files
ilist = find_files(top=o_file, filename_filter=ext)
olist = find_files(top=i_file, filename_filter=ext)
split1 = [os.path.split(file)[1] for file in ilist]
split2 = [os.path.split(file)[1] for file in olist]
RF = []
TLdiff = []
T1L = []
T2L = []
shared_files = []
for file in ilist:
    tree1 = dendropy.Tree.get_from_path(file, 'nexus')
    TL1 = tree1.length()
    T1L.append(TL1)
    if os.path.split(file)[1] in split2:
        shared_files.append(file)
        tree2 = dendropy.Tree.get_from_path(file, 'nexus', taxon_namespace=tree1.taxon_namespace)
        TL2= tree2.length()
        T2L.append(TL2)
        TLdiff.append(TL1-TL2)
        RF.append(treecompare.symmetric_difference(tree1,tree2))


df = pd.DataFrame(shared_files)
df['RF'] = RF
print(df)
df.to_csv('c.csv')

