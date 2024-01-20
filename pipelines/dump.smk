configfile: "pipelines/config.json"

import json
import pandas as pd
import numpy as np
import sys
import os
sys.path.append("scripts")
import FACSus as fs

# Currently deprecated
rule combine_gwas_res:
    input:
        strip,
        manhattan,
        qq
    output:
        combined = "{marker}.combined.png"
    run: 
        fs.combine_gwas_res(wildcards.marker, indir = "/gpfs3/well/ansari/users/gjx698/BOSON_FACS/graphs/facs-host/baseline/all/", save = True)