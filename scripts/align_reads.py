import argparse
import logging
import subprocess
import os
import math
from Bio.Seq import Seq
from datetime import datetime


def mapSmallRNAReadsToGenomeUsingBowtie1(options):
    """
    This function maps the reads to the genome.
    Please check the .alignment file to see how many reads mapped to the genome.
    Please do not issue a samtools flagstat command. The output isn't as accurate.
    """
    # Generate the bowtie index if one is not provided
    if options["bowtie_index"]==None:
        cmd=f"bowtie-build --threads {options["CPU"]} {options["genome"]} {options["output_directory"]}/bowtie1_index"

        # logging.info("Calling bowtie index with the command: %s", cmd)
        print(f"Calling bowtie index with the command: {cmd}")
        result = subprocess.run([cmd], shell=True, capture_output=True, text=True, check=True)
        # logging.info("Bowtie index stdout:\n %s", result.stdout)
        print(result.stdout)
        if result.stderr:
            # logging.warning("Bowtie index stderr:\n %s", result.stderr)
            print(result.stderr)

        bowtie1_index=options["output_directory"]+"/bowtie1_index"

    else:
        bowtie1_index=options["bowtie_index"]
        
    if not os.path.exists(bowtie1_index+".1.ebwtl"):
        large_index=0
    else:
        large_index=1
    
    map_limit = str(options["map_limit"])
    cmd="bowtie "
    if large_index==1:
        cmd+=" --large-index "
    cmd+="-f -m "
    cmd+=map_limit
    cmd+=" -v 0 -a -p "+options["CPU"]+" "
    cmd+=bowtie1_index+" "
    cmd+=options["consolidated_filename"]+" "
    cmd+=" "+options["output_directory"]+"/"+options["grouped_name"]+"_bowtie1.bwt "
    cmd+="2> "+options["output_directory"]+"/"+options["grouped_name"]+"_bowtie1.alignment"

    # logging.info("Calling bowtie alignment with the command: %s", cmd)
    print(f"calling bowtie with command: {cmd}")
    result = subprocess.run([cmd], shell=True, capture_output=True, text=True, check=True)
    print(result.stdout)
    # logging.info("Bowtie alignment stdout:\n %s", result.stdout)
    if result.stderr:
        # logging.warning("Bowtie alignment stderr:\n %s", result.stderr)
        print(result.stderr)



if __name__ == "__main__":
    options_dict = {
        "bowtie_index":None,
        "CPU":"8",
        "genome":"/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/genome/tomato_M82.fasta",
        "output_directory":"results/single_WTM82_rep3_S4",
        # "input_library": "/storage/scratch/s01/users/kal_grimes_tropic_bio/total_samples.fastq.gz",
        # "input_library": "/storage/scratch/s01/jubina_benny_tropic_bio/1115-phasedrna-discovery/phased-rna-discovery/resources/srnaseq/fastq/leaf_S10_R1.fastq.gz",
        "grouped_name": "WTM82_rep3_S4",
        "map_limit":"1",
        }
    # options_dict["input_filename"] = options_dict["input_library"].split("/")[-1].split(".")[0]
    options_dict["consolidated_filename"] = f"{options_dict["output_directory"]}/{options_dict["grouped_name"]}.consolidated.fasta"

    print(options_dict)
    
    print("Mapping reads using bowtie...")
    mapSmallRNAReadsToGenomeUsingBowtie1(options_dict)

    print("Script complete!")