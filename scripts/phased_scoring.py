import argparse
import logging
import subprocess
import sys
import os
import json
import math
from Bio.Seq import Seq
from datetime import datetime


def readMappedData(options,phase):
    """
    Reads in the mapped data into two dictionaries
    """
    whole_mapped_data={}
    size_mapped_data={}
    mapped_data_per_size_per_register={}
    score={}
    readcount={}
    # readseq={}
    alignment_filename = f"{options["output_directory"]}/{options["grouped_name"]}_bowtie1.bwt"
    fhr=open(alignment_filename,"r")

    for line in fhr:
        """
        Data example:
        read_12124_6138 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTT   IIIIIIIIIIIIIIIIIIIII   0
        read_3307062_86 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTAGT        IIIIIIIIIIIIIIIIIIIIIIII        0
        read_4610579_11 -       chr05   3810328 TTTGAGCAAGAAAGTCAGAGTT  IIIIIIIIIIIIIIIIIIIIII  0
        read_6964633_31 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTA  IIIIIIIIIIIIIIIIIIIIII  0
        read_7421569_42 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTAG IIIIIIIIIIIIIIIIIIIIIII 0
        read_9872381_2  +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTAGTAAGA    IIIIIIIIIIIIIIIIIIIIIIIIIIII    0
        read_12413908_9 -       chr05   3810328 TTTGAGCAAGAAAGTCAGAGTTA IIIIIIIIIIIIIIIIIIIIIII 0
        read_14527831_4 -       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTA  IIIIIIIIIIIIIIIIIIIIII  0
        read_26843794_1 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTAGTAAGACG  IIIIIIIIIIIIIIIIIIIIIIIIIIIIII  0
        read_33631320_1 -       chr05   3810327 ATTTGAGCAAGAAAGTCAGAGTTA        IIIIIIIIIIIIIIIIIIIIIIII        0
        read_38589680_2 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTAGTA       IIIIIIIIIIIIIIIIIIIIIIIII       0
        read_39409062_1 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTAGTAAGAC   IIIIIIIIIIIIIIIIIIIIIIIIIIIII   0
        read_42504578_1 -       chr05   3810326 CATTTGAGCAAGAAAGTCAGAGTT        IIIIIIIIIIIIIIIIIIIIIIII        0
        """
        try:
            read_id, strand, chromosome, coordinate, sequence, quality, mapped_times = line.strip().split()
        except ValueError:
            continue

        try:
            coordinate=int(coordinate)
            mapped_times=int(mapped_times)+1
            length=len(sequence)
        except ValueError:
            continue

        """
        If on - strand, account for 2bp overhang and then group together
            should it be positive/negative or sense and anti-sense?
            check against a phased that is reversed? Can't figure out if reversed as not identified?
        """
        if strand=="-":
            coordinate+=2
            # seq=str(Seq(sequence).reverse_complement())
        # else:
        #     seq=sequence
        
        # seperate count info per line from deduping reads in consolidation script
        if 'x' in read_id.split("_")[-1]:
            count=int(read_id.split("_")[-1][1:])
        else:
            count=int(read_id.split("_")[-1])

        """
        For each line in bwt (mapped read unique by sequence) regardless of length,
            keep a count of how many different sequences map to this coordinate
        """
        if chromosome not in whole_mapped_data:
            whole_mapped_data[chromosome]={}
        if coordinate not in whole_mapped_data[chromosome]: 
            whole_mapped_data[chromosome][coordinate]=0
        whole_mapped_data[chromosome][coordinate]+=1
        
        # AFTER "whole" data count collected, dont count anything except correct phase length (eg: 21)
        if phase!=length:
            continue

        """
        For each line in bwt (mapped read unique by sequence) only of relevant read length,
            keep a count of how many different sequences map to this coordinate
        """
        if chromosome not in size_mapped_data:
            size_mapped_data[chromosome]={}
        if coordinate not in size_mapped_data[chromosome]: 
            size_mapped_data[chromosome][coordinate]=0
        size_mapped_data[chromosome][coordinate]+=1


        # calculate a "register" id
        register=coordinate % length

        """
        For each line in bwt (mapped read unique by sequence) only of important length,
            keep a count of how many different sequences map to this register and coordinate
            An attempt to group seqs based on a coordinate keyvalue, and number of cycle calculations?
        """
        # create dict entry
        if chromosome not in mapped_data_per_size_per_register:
            mapped_data_per_size_per_register[chromosome]={}
        if register not in mapped_data_per_size_per_register[chromosome]:
            mapped_data_per_size_per_register[chromosome][register]={}
        if coordinate not in mapped_data_per_size_per_register[chromosome][register]:
            mapped_data_per_size_per_register[chromosome][register][coordinate]=0

        # overwrite dict entry
        mapped_data_per_size_per_register[chromosome][register][coordinate]+=1

        # logic check: all 21bp reads have been consolidated to a single line and therefore at this coordinate it should only exist once
        if mapped_data_per_size_per_register[chromosome][register][coordinate]>2:
            # logging.warning("Trouble with alignments: %s %s %s %s",length,chromosome,register,coordinate)
            print(f"Trouble with alignments: {length} {chromosome} {register} {coordinate}")
        
        if chromosome not in score:
            score[chromosome]={}
        if coordinate not in score[chromosome]:
            score[chromosome][coordinate]=0
        score[chromosome][coordinate]+=count
        
        if chromosome not in readcount:
            readcount[chromosome]={}
        if coordinate not in readcount[chromosome]:
            readcount[chromosome][coordinate]={}
        if strand not in readcount[chromosome][coordinate]:
            readcount[chromosome][coordinate][strand]=count

#         if chromosome not in readseq:
#             readseq[chromosome]={}
#         if coordinate not in readseq[chromosome]:
#             readseq[chromosome][coordinate]={}
#         if strand not in readseq[chromosome][coordinate]:
#             readseq[chromosome][coordinate][strand]=seq

    # with open('test_data_out/whole_mapped_data.json', 'w') as fp:
    #     json.dump(whole_mapped_data, fp)
    # with open('test_data_out/mapped_per_size_per_reg.json', 'w') as fp1:
    #     json.dump(mapped_data_per_size_per_register, fp1)
    # with open('test_data_out/score.json', 'w') as fp2:
    #     json.dump(score, fp2)
    # with open('test_data_out/readcount.json', 'w') as fp3:
    #     json.dump(readcount, fp3)

    return whole_mapped_data, size_mapped_data, mapped_data_per_size_per_register, score, readcount


def siftRegionsOfInterest(options,mapped_data_per_size_per_register,phase,cycle):
    """
    This function will look through the mappings and sift regions of the chromosomes which are of interest to us
    """
    for chromosome in sorted(mapped_data_per_size_per_register):
        # Make separate files for each chromosome
        output_filename = f"{options["output_directory_per_run"]}/{options["grouped_name"]}_{str(phase)}_{str(cycle)}_{chromosome}.regionsOfInterest.full"
        fhw=open(output_filename,"w")

        for register in sorted(mapped_data_per_size_per_register[chromosome]):
            start,end=0,0
            for coordinate in sorted(mapped_data_per_size_per_register[chromosome][register]):
                if start == 0:
                    start = coordinate
                elif end == 0:
                    if coordinate-start < phase*cycle:
                        end = coordinate
                    else:
                        start = coordinate
                else:
                    if coordinate-end < phase*cycle:
                        end = coordinate
                    else:
                        fhw.write(str(register)+"\t"+str(start)+"\t"+str(end+phase-1)+"\n")
                        end=0
                        start=coordinate
            if end!=0:
                fhw.write(str(register)+"\t"+str(start)+"\t"+str(end+phase-1)+"\n")

        fhw.close()


def parseGFF(gff_path):
    # read gff to dict chr1: [{start: 1, end: 123, strand: +}, {start: 234, end: 345, strand: -}]
    with open(gff_path, "rt") as infile:
        lines = infile.read().split("\n")
        gff_dict = {}

    for line in lines:
        if line != "" and not line.startswith("#"):
            chr = line.split("\t")[0]
            start = line.split("\t")[3]
            end = line.split("\t")[4]
            strand = line.split("\t")[6]
            if chr not in gff_dict.keys():
                gff_dict[chr]=[
                    {"start": start,
                    "end": end,
                    "strand": strand}
                ]
            else:
                gff_dict[chr].append(
                    {"start": start,
                    "end": end,
                    "strand": strand}
                    )
    return gff_dict


def filterByRegion(options, mapped_data_per_size_per_register, phase, cycle):
    if options["regions_file"]:
        # logging.info("Regions file GFF option provided: %s", options.regions_file)
        print(f"Regions file GFF option provided: {options["regions_file"]}")
        gff_data = parseGFF(options["regions_file"])
        fail_count = 0

        # filter files based off chromosome
        for chromosome in sorted(mapped_data_per_size_per_register):
            if chromosome in gff_data.keys():
                output_filename=options["output_directory_per_run"]+"/"+options["grouped_name"]+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest"
                input_filename=f"{output_filename}.full"

                # find if any data is within regions of interest and keep it
                data = []
                with open(input_filename, "r") as input:
                    for line in input:
                        line_data = line.strip("\n").split("\t")
                        for region in gff_data[chromosome]:
                            # gff start is before ROI end, and gff end is after ROI start
                            if region["start"] <= line_data[2] and region["end"] >= line_data[1]:
                                if line not in data:
                                    data.append(line)    

                if data:
                    # write data to new file ROI
                    with open(output_filename, "w") as file:
                        for line in data:
                            file.write(line)
                else:
                    # logging.warning("Chromosome %s has no regions which made it through the filtering.", chromosome)
                    print(f"Chromosome {chromosome} has no regions which made it through the filtering.")
                    fail_count += 1
            
            else:
                # add to fail count
                # logging.warning("Chromosome %s not in regions file and has been filtered out.", chromosome)
                print(f"Chromosome {chromosome} not in regions file and has been filtered out.")
                fail_count += 1

        if fail_count == len(mapped_data_per_size_per_register):
            # capture for if filtering removes all regions
            # logging.error("No regions survived filtering based on gff regions inputted.")
            print("No regions survived filtering based on gff regions inputted.")
            print("The program had to terminate prematurely....Please check log file for more details.")
            sys.exit()

    else:
        # rename/copy files from roi.full to roi if filtering not needed
        for chromosome in sorted(mapped_data_per_size_per_register):
            filename=f"{options["output_directory_per_run"]}/{options["grouped_name"]}_{str(phase)}_{str(cycle)}_{chromosome}.regionsOfInterest"
            os.rename(f"{filename}.full", filename)


def nCr(n,r):
    if (n-r)<0 or n<1 or r<1:
        return 1
    p = math.factorial(n)//(math.factorial(r)*math.factorial(n-r))
    # print(f"n value is {n}, r value is {r}")
    return p


def computePValues(options, whole_mapped_data, size_mapped_data, mapped_data_per_size_per_register, phase, cycle):
    """
    Computes the P-values using a Hypergeometric distribution
    M = total number of de-duped reads (coordinates) in the fastq that is phase length (21)
    n = number of de-duped reads that cover the phased window/phasiRNA region (phase length 21)
    N = number of all possible phasing positions per that phasiRNA region (region // phase lenth (21); floor division)
    k = number of phased positions that are covered by a read
    """
    min_reads_mapped_to_a_phased_register=3
    min_reads_in_a_window=10
    chromosome_hits=[]
    for chromosome in sorted(mapped_data_per_size_per_register):
        chromosome_hits.append(chromosome)
        try:
            fhr=open(f"{options["output_directory_per_run"]}/{options["grouped_name"]}_{str(phase)}_{str(cycle)}_{chromosome}.regionsOfInterest","r")
            fhw=open(f"{options["output_directory_per_run"]}/{options["grouped_name"]}_{str(phase)}_{str(cycle)}_{chromosome}.regionsOfInterest.concentrated","w")
        except FileNotFoundError:
            # filter by gff removes some files so if you can't identify them then move onto the next chromosome
            continue
        for line in fhr:
            register,start,end=line.strip().split()
            register=int(register)
            start=int(start)
            end=int(end)
            
            begin=start
            sys.stdout.flush()
            print("p-val info:")
            print(f"Register: {register}, start = {start}, end = {end}")
        
            while begin+(phase*min_reads_mapped_to_a_phased_register) <= end+1:
                finish=begin+(phase*cycle)-1
                num_all_reads=0
                num_all_size_reads=0
                k=0
                n=0
                m=cycle*2
                pvalue=0

                stats_N = (finish - begin) // phase
                print(f"N = {stats_N}")

                for i in range(begin,finish+1):
                    try:
                        # number of phase size reads mapped to register
                        k+=mapped_data_per_size_per_register[chromosome][register][i]
                    except KeyError:
                        pass

                stats_k = k
                print(f"k = {stats_k}")

                if k<min_reads_mapped_to_a_phased_register: 
                    begin+=phase
                    continue
                
                for i in range(begin,finish+1):
                    try:
                        num_all_reads+=whole_mapped_data[chromosome][i]
                    except KeyError:
                        pass

                if num_all_reads<min_reads_in_a_window:
                    begin+=phase
                    continue

                for i in range(begin,finish+1):
                    try:
                        num_all_size_reads+=size_mapped_data[chromosome][i]
                    except KeyError:
                        pass

                stats_M = num_all_size_reads
                print(f"M = {stats_M}")

                # register_i is an iterator different from register
                for register_i in sorted(mapped_data_per_size_per_register[chromosome]):
                    for i in range(begin,finish+1):
                        try:
                            n+=mapped_data_per_size_per_register[chromosome][register_i][i]
                        except KeyError:
                            pass

                stats_n = n
                print(f"n = {stats_n}")

                # if n/num_all_reads<0.3:
                #     begin+=phase
                #     continue

                for x in range(k,m+1):
                    numerator=nCr((phase-1)*m,n-x)*nCr(m,x)
                    pvalue+=numerator
                    # print(f"x is {x}, numerator = {numerator}, pval is {pvalue}")
                denominator=nCr(phase*m,n)
                pvalue=pvalue/denominator
                # print(f"denom is {denominator}, pval is {pvalue}")
                

                if pvalue>=float(options["pvalue_cutoff"]):
                    begin+=phase
                    continue

                stuffs_to_be_printed_to_file=[register,begin,finish,k,n,m,num_all_reads,n/num_all_reads,pvalue]
                fhw.write("\t".join(map(str,stuffs_to_be_printed_to_file))+"\n")
                sys.stdout.flush()
                begin+=phase


def generatePositivePHASLoci(options,whole_mapped_data,phase,cycle):
    """
    Generate a file with the set of positive loci on all the chromosomes
    """
    out_filename=f"{options["output_directory_per_run"]}/{options["grouped_name"]}_{str(phase)}_{str(cycle)}.positive_phase_loci"
    fhw=open(out_filename,"w")

    for chromosome in sorted(whole_mapped_data):

        filename=f"{options["output_directory_per_run"]}/{options["grouped_name"]}_{str(phase)}_{str(cycle)}_{chromosome}.regionsOfInterest.concentrated"
        flag_reg=1000
        window_start=0
        window_end=0

        try:
            fhr=open(filename,"r")
        except FileNotFoundError:
            continue

        for line in fhr:
            register,start,end=map(int,line.strip().split()[:3])

            if register==flag_reg:
                if window_end>start:
                    window_end=end
                else:
                    fhw.write(chromosome+"\t"+str(window_start)+"\t"+str(window_end)+"\n")
                    window_start=start
                    window_end=end

            else:
                if flag_reg!=1000:
                    fhw.write(chromosome+"\t"+str(window_start)+"\t"+str(window_end)+"\n")
                window_start=start
                window_end=end
                flag_reg=register

        fhr.close()
        fhw.write(chromosome+"\t"+str(window_start)+"\t"+str(window_end)+"\n")
    fhw.close()
    

def readFastaFile(filename):
    """
    Reads in a fasta file and returns a dictionary
    The keys in the dictionary is same as the fasta header
    for each sequence upto the first space.
    """
    info={}
    fhr=open(filename,"r")

    while(True):
        line=fhr.readline()
        if not line: break

        if(">" in line):
            try:
                info[line.strip()[1:].split()[0]]=fhr.readline().strip()
            except ValueError:
                pass
    return info
        

def generatePhasingScore(options, phase, cycle, score, readcount):
    """
    Generates phasing scores for the phased loci
    """
    # score,readcount,readseq=readDataForPhasingScoreComputation(options,phase)
    phased_loci_filename=f"{options["output_directory_per_run"]}/{options["grouped_name"]}_{str(phase)}_{str(cycle)}.positive_phase_loci"
    final_phase_loci=f"{options["output_directory_per_run"]}/{options["grouped_name"]}_{str(phase)}_{str(cycle)}.phasing_score_phase_loci"
    fhr=open(phased_loci_filename,"r")
    out4=open(final_phase_loci,"w")

    for line in fhr:
        chromosome,ss,ee=line.strip().split()
        ss=int(ss)
        ee=int(ee)
        phasing_score_filename=options["output_directory_per_run"]+"/"+str(phase)+"_"+str(chromosome)+"_"+str(ss)+"_"+str(ee)+".phasing_score"
        abundance_score_filename=options["output_directory_per_run"]+"/"+str(phase)+"_"+str(chromosome)+"_"+str(ss)+"_"+str(ee)+".abundance"
        out=open(phasing_score_filename,"w")
        out2=open(abundance_score_filename,"w")
        score_count={}

        for site in range(ss,ee+1):
            start=site-(phase*4)
            end=site+(phase*5)-1
            max_within_site,max_within_count,all_scores=0,0,0
            for cor in range(start,end+1):
                if cor not in score[chromosome]:continue
                all_scores+=score[chromosome][cor]
                for i in readcount[chromosome][cor]:
                    if max_within_count<readcount[chromosome][cor][i]:
                        max_within_site=cor
                        max_within_count=readcount[chromosome][cor][i]
            all_scores-=max_within_count
            P,k=0,0
            s=start

            while s<end:
                if s not in score[chromosome]:
                    s+=phase
                    continue
                if score[chromosome][s]!=0:
                    P+=score[chromosome][s]
                    k+=1
                    if s == max_within_site:
                        P-=max_within_count 
                s+=phase
            U=all_scores-P
            
            if k>=3:
                phas_score=math.log((1+(10*(P/(1+U))))**(k-2))
            else:
                phas_score=0

            out.write(str(site)+"\t"+str(phas_score)+"\n")
            out4.write(chromosome+"\t"+str(site)+"\t"+str(phas_score)+"\n")

            if chromosome not in score_count:
                score_count[chromosome]={}
            if site not in score_count[chromosome]:
                score_count[chromosome][site]=phas_score
            if site in readcount[chromosome] and '+' in readcount[chromosome][site] and readcount[chromosome][site]['+']!=0:
                out2.write(str(site)+"\t"+str(readcount[chromosome][site]['+'])+"\n")
            if site in readcount[chromosome] and '-' in readcount[chromosome][site] and readcount[chromosome][site]['-']!=0:
                out2.write(str(site)+"\t-"+str(readcount[chromosome][site]['-'])+"\n")

        out.close()
        out2.close()
    out4.close()


if __name__ == "__main__":
    phase = 21
    cycle = 22
    options_dict = {
        "output_directory": "results/wtm82_fastq",
        "regions_file": "/home/kal_grimes_tropic_bio/proj/phasi/tas4.gff",
        # "input_library": "/storage/scratch/s01/jubina_benny_tropic_bio/1115-phasedrna-discovery/phased-rna-discovery/resources/srnaseq/fastq/leaf_S10_R1.fastq.gz",
        "pvalue_cutoff": "0.01",
        "grouped_name": "wtm82",
    }
    # options_dict["input_filename"] = options_dict["input_library"].split("/")[-1].split(".")[0]

    if options_dict["regions_file"]:
        regions_basename = os.path.basename(options_dict["regions_file"]).split(".")[0]
        options_dict["output_directory_per_run"] = f"{options_dict["output_directory"]}/phase_{str(phase)}_cycle_{str(cycle)}_{regions_basename}"

    else:
        options_dict["output_directory_per_run"] = f"{options_dict["output_directory"]}/phase_{str(phase)}_cycle_{str(cycle)}"

    # logging.info("Reading mapped data for phase length %s...", phase)
    whole_mapped_data, size_mapped_data, mapped_data_per_size_per_register, score_dict, readcount_dict = readMappedData(options_dict,phase)

    # create dir for results
    cmd=f"mkdir {options_dict["output_directory_per_run"]}"
    os.system(cmd)

    # logging.info("Sifting regions of interest...")
    print("Sifting regions of interest...")
    siftRegionsOfInterest(options_dict,mapped_data_per_size_per_register,phase,cycle)

    # logging.info("Filtering results by region")
    print("Filtering results by region")
    filterByRegion(options_dict,mapped_data_per_size_per_register, phase, cycle)

    # logging.info("Calculating p values...")
    print("Calculating p values...")
    computePValues(options_dict, whole_mapped_data, size_mapped_data, mapped_data_per_size_per_register, phase,cycle)

    # logging.info("Generating phased loci...")
    print("Generating phased loci...")
    generatePositivePHASLoci(options_dict,whole_mapped_data,phase,cycle)

    # logging.info("Generating phasing scores...")
    print("Generating phasing scores...")
    generatePhasingScore(options_dict, phase, cycle, score_dict, readcount_dict)

    # logging.info("Creating plots...")
    print("Creating plots...")
    cmd=f"mkdir -p {options_dict["output_directory_per_run"]}/plots"
    os.system(cmd)

    cmd="Rscript --vanilla plot.R "
    cmd+=" "+options_dict["output_directory_per_run"]
    cmd+=" "+str(phase)
    cmd+=" "+str(cycle)
    # logging.info("Calling Rscript plot.R with the command: %s", cmd)
    print(f"Calling Rscript plot.R with the command: {cmd}")
    result = subprocess.run([cmd], shell=True, capture_output=True, text=True, check=True)
    # logging.info("R script stdout:\n %s", result.stdout)
    print(result.stdout)
    if result.stderr:
        # logging.warning("R script stderr:\n %s", result.stderr)
        print(result.stderr)

    print("Script complete!")
