#!/usr/bin/env python
import re
from collections import defaultdict
import sys


'''

Searches the .hyb file line by line for any microRNA then counts each one


Example hyb line:

365_33519       TAGCTTATCAGACTGATGTAATACTGTCTGGTAAAACCGT        .       hsa-miR-21-5p_MIMAT0000076_Homo_sapiens_miR-21-5p       1       19      1       19      2.6e-02 hsa-miR-429_MIMAT0001536_Homo_sapiens_miR-429   19      40      1       22      4.8e-04


'''

input_file = sys.argv[1]
d = defaultdict(int)
comp = re.compile(r'.*?\t([a-z]{3}-miR-.*?)_.*?') # Groups after tab mir match. Non greedy part is important because there could be two mirs
lines = 0
with open(input_file) as infile:
    for line in infile:
        lines+=1
        matches  = set(comp.findall(line))
        for item in matches:
            d[item]+=1
        if lines%1000 == 0:
            print(lines, "lines")
order = [(mir,counts) for mir, counts in d.items()]
order.sort(key=lambda x:x[1])
with open("%s.mirna_counts" % input_file, "w") as outfile:
    outfile.write("microRNA\tCounts\n")
    for mir, counts in order:
        outfile.write("%s\t%i\n" % (mir, counts))
        

