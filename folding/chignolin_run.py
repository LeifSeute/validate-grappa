#! bin/bash

name = 'chignolin'

BOX_SIZE=4

import os

AMBER_IDS = ["23784", "53435", "98634"]
AMBER_JOB_NAMES = ["amber1", "amber2", "amber3"]

os.system(f"cd {name}_amber99sbildn")

for i, (id, jobname) in enumerate(zip(AMBER_IDS, AMBER_JOB_NAMES)):
    print(f"Submitting job {i} with ID {id} and jobname {jobname}...")
    os.system(f'sbatch -J {jobname} submit_continue.sh {id} {BOX_SIZE}')
