import argparse
import logging
import subprocess
import gzip
import sys
import os
from datetime import datetime

"""
Script creates part of the phasi scoring pipeline.
This script trims adapters and then dedups the fastq file into a fasta file which is deduped by read sequence but keeps a tally of count that seq was seen
To run:
sbatch -o slurm.stdout -e slurm.err -J findphasi -c 4 --mem=64G --wrap "python scripts/consolidate_reads.py"
"""


def trimAdapters(options, log_filepath, input_filename, output_filename):
    """
    Trim adapters from sequence
    """
    # input_filename=options["input_library"]
    # output_filename=options["adapter_trimmed_filename"]
    cpu_request = str(options["CPU"])
    cmd=f"trimmomatic SE -threads {cpu_request} {input_filename} {output_filename} ILLUMINACLIP:/home/kal_grimes_tropic_bio/proj/phasi/findPhasiRNAs/resources/adapters.fasta:2:30:10"

    # logging.info("Calling Trimmomatic SE with the command: %s", cmd)
    print(f"Calling Trimmomatic SE with the command: {cmd}")
    result = subprocess.run([cmd], shell=True, capture_output=True, text=True, check=True)
    # logging.info("Trimmomatic SE stdout:\n %s", result.stdout)
    print(result.stdout)
    # if result.stderr:
    print(result.stderr)
        # logging.warning("Trimmomatic SE stderr:\n %s", result.stderr)

    # sample_name = options["input_filename"].split("/")[-1].split(".")[0]
    # fasta_outname = f"{options["output_directory"]}/{sample_name}.fa"
    # cmd=f"sed -n '1~4s/^@/>/p;2~4p' {output_filename} > {fasta_outname}"
    # os.system(cmd)


def consolidateReads(input_filename, data):
    """
    Select unique reads and combine their counts. Eliminates the quality values.
    """
    def _process(lines=None):
        ks = ['name', 'sequence', 'optional', 'quality']
        return {k: v for k, v in zip(ks, lines)}

    # input_filename=options["adapter_trimmed_filename"]
    # output_filename=options["consolidated_filename"]

    # data={}
    n = 4
    with gzip.open(input_filename, 'rt') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            # create block of 4 lines (a record) and process each
            if len(lines) == n:
                record = _process(lines)
                lines = []

                # for each record (block of 4 dict), create or add sequence to data dict
                if record["sequence"] not in data:
                    data[record["sequence"]]=1
                else:
                    data[record["sequence"]]+=1

    return data

def write_consolidated_file(data, output_filename):
    with open(output_filename, "w") as outfile:
        for seq_num,seq in enumerate(data):
            outfile.write(f">read_{str(seq_num+1)}_{str(data[seq])}\n{seq}\n")


if __name__ == "__main__":
    # input_library = [
    #     "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/M82C_1_S10_R1_lane1_trimmed.fastq.gz",
    #     "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/M82C_2_S22_R1_lane1_trimmed.fastq.gz",
    #     # "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/SRR8466936_lane1_trimmed.fastq.gz",
    #     # "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/SRR8466937_lane1_trimmed.fastq.gz",
    #     # "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/SRR8466952_lane1_trimmed.fastq.gz",
    #     # "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/SRR8466953_lane1_trimmed.fastq.gz",
    #     "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/WTM82_rep1_S8_lane1_trimmed.fastq.gz",
    #     "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/WTM82_rep2_S7_lane1_trimmed.fastq.gz",
    #     "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/WTM82_rep3_S4_lane1_trimmed.fastq.gz",
    # ]
    input_library = [
        # "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/M82C_1_S10_R1_lane1_trimmed.fastq.gz",
        # "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/M82C_2_S22_R1_lane1_trimmed.fastq.gz",
        # "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/WTM82_rep1_S8_lane1_trimmed.fastq.gz",
        # "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/WTM82_rep2_S7_lane1_trimmed.fastq.gz",
        "/storage/scratch/s01/vered/projects/2024/0925-phased-rna-tomato-M82/srna/results/trimmed/WTM82_rep3_S4_lane1_trimmed.fastq.gz",
    ]
    options_dict = {
        # "input_library": "/storage/scratch/s01/jubina_benny_tropic_bio/1115-phasedrna-discovery/phased-rna-discovery/resources/srnaseq/fastq/leaf_S10_R1.fastq.gz",
        # "input_library": "/storage/scratch/s01/users/kal_grimes_tropic_bio/total_samples.fastq.gz",
        "output_directory": "results/single_WTM82_rep3_S4",
        "grouped_name": "WTM82_rep3_S4",
        "CPU": "4",
    }
    # options_dict["input_filename"] = options_dict["input_library"].split("/")[-1].split(".")[0]
    # options_dict["adapter_trimmed_filename"] = f"{options_dict["output_directory"]}/{options_dict["input_filename"]}_adapter_trimmed.fastq"
    consolidated_filename = f"{options_dict["output_directory"]}/{options_dict["grouped_name"]}.consolidated.fasta"

    print(options_dict)
    print(f"Input library is: {input_library}")
    print(f"Consolidated_filename is: {consolidated_filename}")

    if os.path.exists(options_dict["output_directory"])==False:
        cmd=f"mkdir -p {options_dict["output_directory"]}"
        os.system(cmd)
    else:
        print(f"Output directory {options_dict["output_directory"]} already exists. Please delete or rename before re-running.")
        print("The script had to terminate prematurely. Please check log file for more details.")
        sys.exit()
    
    out_data = {}
    for trimmed_fastq_file in input_library:

        input_filename = trimmed_fastq_file.split("/")[-1].split(".")[0]
        # adapter_trimmed_filename = f"{options_dict["output_directory"]}/{input_filename}_adapter_trimmed.fastq"

        # print(f"Trimming adapters in {fastq_file}...")
        # trimAdapters(options_dict, "path/placeholder", fastq_file, adapter_trimmed_filename)

        print(f"Consolidating reads in {trimmed_fastq_file}...")
        data = consolidateReads(trimmed_fastq_file, out_data)

        # overwrite out_data to append more
        out_data = data

    print("writing consolidated file...")
    write_consolidated_file(out_data, consolidated_filename)

    print("Script complete!")

