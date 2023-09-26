import sys
import FACSus as fs
import pandas as pd

marker = str(sys.argv[1])
indir = str(sys.argv[2])

fs.parse_plink(indir + marker + '.assoc.linear', True, indir + marker + '.baseline.gwas.txt')
