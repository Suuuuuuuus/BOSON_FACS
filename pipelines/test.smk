configfile: "pipelines/config.json"
include: "auxiliary.smk"

from os.path import exists
import json
import pandas as pd
import numpy as np
import sys
sys.path.append("scripts")
import FACSus as fs

rule test:
    output:
        vcf = "results/tmp/{id}.{chr}.txt"
    shell: """
        echo {wildcards.id}
    """