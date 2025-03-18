#######################################################################################################
# Developed by: Sagnik Banerjee, Kalon Grimes
# Version: 1.1
# 
# This program will analyze locations of phasiRNAs. Please
# launch python findPhasiRNAs.py --help for more details of each possible functionality
# of the software. 
#######################################################################################################

import argparse
import logging
import subprocess
import sys
import os
import math
from Bio.Seq import Seq
from datetime import datetime

def create_logs():
    # make log dir
    try:
        os.mkdir("logs")
    except:
        pass

    timenow = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    log_filepath = f"logs/{timenow}_findphasi.log"

    # initiate logging
    logging.basicConfig(
        filename=log_filepath,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
    )

    logging.info("Logging initiated.")
    print(f"Logging initiated: {log_filepath}")
    return log_filepath



def parseCommandLineArguments():
    """
    Parses the arguments provided through command line.
    Launch python analyzePhasiRNAs.py --help for more details
    """
    parser = argparse.ArgumentParser(prog="findPhasiRNAs.py",description="findPhasiRNAs can be used to find genomic locations where phasing occurs. ")
    optional_arg = parser.add_argument_group("Optional Arguments")
    required_arg = parser.add_argument_group("Required Arguments")
    genome_mutex = parser.add_mutually_exclusive_group(required = True)
    
    required_arg.add_argument("--input_library","-i",help="Specify the name of the file which has the small-RNA reads. This option is mutually exclusive with --consolidated_library")
    genome_mutex.add_argument("--genome","-g",help="Specify the name of the genome fasta file of the organism. Please note that the program will not be able to handle multiple fasta files. ")
    genome_mutex.add_argument("--bowtie_index","-bindex",help="Provide the bowtie index. This argument is optional. If no index is provided then the software will generate one.")
    required_arg.add_argument("--output_directory","-out",type=lambda x: "results/" + x, help="Specify an output directory to which all the generated files will be housed. This includes the log file which can be later checked. Please make sure that there are sufficient permissions to create the output directory. The program will throw an error if creation of the output directory fails. If the directory already exists then its contents will be overwritten without warning. This directory will contain the summary file containing the details of the execution",required=True)
    optional_arg.add_argument("--regions_file", "-r", help="Optional GFF file to filter results to only corresponding regions that have some overlap.")
    optional_arg.add_argument("--small_rna_size","-srnasize",nargs="+",help="Specify the size of the small RNA that you wish to analyze. You can enter more than one possible size.",default=["21"])
    optional_arg.add_argument("--number_of_cycles","-numcycles",nargs="+",help="Specify the number of cycles you wish to analyze with. You can enter multiple number of number of cycles. The accepted values are 9, 10, 11, 12 and 13",default=["9"])
    optional_arg.add_argument("--pvalue_cutoff","-p",help="Enter the p-value cut off",default=0.05)
    optional_arg.add_argument("--clean_up","-c",help="Set this to 1 if you wish to clean up all the intermediate files. The program will keep all temporary files by default.",default=0)
    optional_arg.add_argument("--CPU","-n",help="Provide the number of CPUs to be used. Default is 1.",default="1")
    optional_arg.add_argument("--map_limit","-mapl",help="Specify the mapping limit. Only reads which are mapped at most -mapl times will be considered. The default is 1. The maximum number of alignments allowed for a single read is 10. ",default=1)
    optional_arg.add_argument("--force",help="Overwrite contents of output directory if it exists.",default=0)
       
    # Supressed arguments
    parser.add_argument("--input_filename","-ifname",help=argparse.SUPPRESS)
    parser.add_argument("--input_path","-ipath",help=argparse.SUPPRESS)
    parser.add_argument("--consolidated_filename","-cfname",help=argparse.SUPPRESS)
    parser.add_argument("--adapter_trimmed_filename","-atfname",help=argparse.SUPPRESS)
    parser.add_argument("--output_directory_per_run","-output_directory_per_run",help=argparse.SUPPRESS)
    
    return parser.parse_args()
    

def analyzeCommandLineArguments(options):
    """
    Performs checks on the validity of the arguments provided 
    through the command line.
    """
    flag=0
    if os.path.exists(options.output_directory)==False:
        cmd="mkdir -p "+options.output_directory
        os.system(cmd)
    else:
        if options.force==0:
            logging.error("Output directory '%s' already exists. Please re-run the program with -f 1 to force rewrite of the directory", options.output_directory)
            flag=1
        else:
            cmd="rm -rf "+options.output_directory
            os.system(cmd)
            cmd="mkdir -p "+options.output_directory
            os.system(cmd)
        
    cmd="touch "+options.output_directory+"/Log.out"
    os.system(cmd)

    if options.bowtie_index == None:
        logging.warning("No bowtie index provided. Proceeding to building index...")
        if not os.path.exists(options.genome):
            logging.error("The genome file '%s' does not exist.", options.genome)
            flag=1

    if options.input_library == None:
        logging.error("The input file '%s' does not exist", options.input_library)
        flag=1

    if options.regions_file:
        if os.path.exists(options.regions_file)==False:
            logging.error("The input file '%s' does not exist", options.regions_file)
            flag=1

    for ele in options.number_of_cycles:
        if ele not in ["9","10","11","12","13"]:
            logging.error("'%s' is not a valid choice of cycle. Valid choices are 9, 10, 11, 12 and 13.", ele)
            flag=1
    
    # throw error and exit if any safeties are not satisfied    
    if flag==1:
        logging.info("The script terminated due to errors.")
        print("The script had to terminate prematurely. Please check log file for more details.")
        sys.exit()

    if options.input_library.split(".")[-1]=="fq" or options.input_library.split(".")[-1]=="fastq":
        options.input_filename=options.output_directory+"/"+options.input_library.split("/")[-1].split(".")[0]+".fa"
    
    options.input_path="/".join(options.input_library.split("/")[:-1])
    options.input_filename=options.input_library.split("/")[-1].split(".")[0]
    options.consolidated_filename=options.output_directory+"/"+options.input_filename+".consolidated.fasta"
    options.adapter_trimmed_filename=options.output_directory+"/"+options.input_filename+"_adapter_trimmed.fastq"
        
    options.small_rna_size=list(map(int,options.small_rna_size))
    options.number_of_cycles=list(map(int,options.number_of_cycles))
    options.pvalue_cutoff=float(options.pvalue_cutoff)
    options.map_limit=int(options.map_limit)
    return options


def readFastqFile(filename):
    reads={}
    fhr=open(filename,"r")
    while True:
        line=fhr.readline()
        if not line:
            break
        reads[line.split()[0][1:]]=[fhr.readline().strip(),fhr.readline().strip(),fhr.readline().strip()]
    return reads


def trimAdapters(options, log_filepath):
    """
    Trim adapters from sequence
    """
    input_filename=options.input_library
    output_filename=options.adapter_trimmed_filename
    cmd="trimmomatic SE -threads "+str(options.CPU)+" "+input_filename
    cmd+=" "+output_filename
    cmd+=" ILLUMINACLIP:resources/adapters.fasta:2:30:10 "

    logging.info("Calling Trimmomatic SE with the command: %s", cmd)
    result = subprocess.run([cmd], shell=True, capture_output=True, text=True, check=True)
    logging.info("Trimmomatic SE stdout:\n %s", result.stdout)
    if result.stderr:
        logging.warning("Trimmomatic SE stderr:\n %s", result.stderr)

    cmd="sed -n '1~4s/^@/>/p;2~4p' "+options.adapter_trimmed_filename+" > "+options.output_directory+"/"+options.input_library.split("/")[-1].split(".")[0]+".fa"
    os.system(cmd)


def consolidateReads(options):
    """
    Select unique reads and combine their counts. Eliminates the quality values.
    """
    input_filename=options.adapter_trimmed_filename
    output_filename=options.consolidated_filename
    fhw=open(output_filename,"w")
    fhr=open(input_filename,"r")
    data={}

    while True:
        line=fhr.readline().strip()
        if not line:
            break

        id=line
        seq=fhr.readline().strip()
        useless=fhr.readline()
        quality=fhr.readline()

        if seq not in data:
            data[seq]=1
        else:
            data[seq]+=1

    for seq_num,seq in enumerate(data):
        fhw.write(f">read_{str(seq_num+1)}_{str(data[seq])}\n{seq}\n")

    fhw.close()
    

def mapSmallRNAReadsToGenomeUsingBowtie1(options):
    """
    This function maps the reads to the genome.
    Please check the .alignment file to see how many reads mapped to the genome.
    Please do not issue a samtools flagstat command. The output isn't as accurate.
    """
    # Generate the bowtie index if one is not provided
    if options.bowtie_index==None:
        cmd="bowtie-build"
        cmd+=" --threads "+options.CPU+" "
        cmd+=options.genome+" "
        cmd+=options.output_directory+"/bowtie1_index"

        logging.info("Calling bowtie index with the command: %s", cmd)
        result = subprocess.run([cmd], shell=True, capture_output=True, text=True, check=True)
        logging.info("Bowtie index stdout:\n %s", result.stdout)
        if result.stderr:
            logging.warning("Bowtie index stderr:\n %s", result.stderr)

        bowtie1_index=options.output_directory+"/bowtie1_index"

    else:
        bowtie1_index=options.bowtie_index
        
    if not os.path.exists(bowtie1_index+".1.ebwtl"):
        large_index=0
    else:
        large_index=1
    
    cmd="bowtie "
    if large_index==1:
        cmd+=" --large-index "
    cmd+=" -f -m "
    cmd+=str(options.map_limit)
    cmd+=" -v 0 -a -p "+options.CPU+" "
    cmd+=bowtie1_index+" "
    cmd+=options.consolidated_filename+" "
    cmd+=" "+options.output_directory+"/"+options.input_filename+"_bowtie1.bwt "
    cmd+=" 2> "+options.output_directory+"/"+options.input_filename+"_bowtie1.alignment "

    logging.info("Calling bowtie alignment with the command: %s", cmd)
    result = subprocess.run([cmd], shell=True, capture_output=True, text=True, check=True)
    logging.info("Bowtie alignment stdout:\n %s", result.stdout)
    if result.stderr:
        logging.warning("Bowtie alignment stderr:\n %s", result.stderr)


def readMappedData(options,phase):
    """
    Reads in the mapped data into two dictionaries
    """
    whole_mapped_data={}
    mapped_data_per_size_per_register={}
    score={}
    readcount={}
    # readseq={}
    alignment_filename = f"{options.output_directory}/{options.input_filename}_bowtie1.bwt"
    fhr=open(alignment_filename,"r")

    for line in fhr:
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

        if strand=="-":
            coordinate+=2
            # seq=str(Seq(sequence).reverse_complement())
        # else:
        #     seq=sequence
            
        if chromosome not in whole_mapped_data:
            whole_mapped_data[chromosome]={}
        if coordinate not in whole_mapped_data[chromosome]: 
            whole_mapped_data[chromosome][coordinate]=0
        whole_mapped_data[chromosome][coordinate]+=1
        
        if phase!=length:
            continue
        if chromosome not in mapped_data_per_size_per_register:
            mapped_data_per_size_per_register[chromosome]={}
        register=coordinate % length
        if register not in mapped_data_per_size_per_register[chromosome]:
            mapped_data_per_size_per_register[chromosome][register]={}
        if coordinate not in mapped_data_per_size_per_register[chromosome][register]:
            mapped_data_per_size_per_register[chromosome][register][coordinate]=0
        mapped_data_per_size_per_register[chromosome][register][coordinate]+=1
        if mapped_data_per_size_per_register[chromosome][register][coordinate]>2:
            logging.warning("Trouble with alignments: %s %s %s %s",length,chromosome,register,coordinate)
        
        if 'x' in read_id.split("_")[-1]:
            count=int(read_id.split("_")[-1][1:])
        else:
            count=int(read_id.split("_")[-1])
        
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

    return whole_mapped_data, mapped_data_per_size_per_register, score, readcount


def siftRegionsOfInterest(options,mapped_data_per_size_per_register,phase,cycle):
    """
    This function will look through the mappings and sift regions of the chromosomes which are of interest to us
    """
    for chromosome in sorted(mapped_data_per_size_per_register):
        # Make separate files for each chromosome
        output_filename=options.output_directory_per_run+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest.full"
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
    if options.regions_file:
        logging.info("Regions file GFF option provided: %s", options.regions_file)
        gff_data = parseGFF(options.regions_file)
        fail_count = 0

        # filter files based off chromosome
        for chromosome in sorted(mapped_data_per_size_per_register):
            if chromosome in gff_data.keys():
                output_filename=options.output_directory_per_run+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest"
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
                    logging.warning("Chromosome %s has no regions which made it through the filtering.", chromosome)
                    fail_count += 1
            
            else:
                # add to fail count
                logging.warning("Chromosome %s not in regions file and has been filtered out.", chromosome)
                fail_count += 1

        if fail_count == len(mapped_data_per_size_per_register):
            # capture for if filtering removes all regions
            logging.error("No regions survived filtering based on gff regions inputted.")
            print("The program had to terminate prematurely....Please check log file for more details.")
            sys.exit()

    else:
        # rename/copy files from roi.full to roi if filtering not needed
        for chromosome in sorted(mapped_data_per_size_per_register):
            filename=options.output_directory_per_run+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest"
            os.rename(f"{filename}.full", filename)


def nCr(n,r):
    if (n-r)<0 or n<1 or r<1:
        return 1
    return math.factorial(n)/(math.factorial(r)*math.factorial(n-r))


def computePValues(options,whole_mapped_data,mapped_data_per_size_per_register,phase,cycle):
    """
    Computes the P-values using a Hypergeometric distribution
    """
    min_reads_mapped_to_a_phased_register=3
    min_reads_in_a_window=10
    chromosome_hits=[]
    for chromosome in sorted(mapped_data_per_size_per_register):
        chromosome_hits.append(chromosome)
        try:
            fhr=open(options.output_directory_per_run+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest","r")
            fhw=open(options.output_directory_per_run+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest.concentrated","w")
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
            while begin+(phase*min_reads_mapped_to_a_phased_register) <= end+1:
                finish=begin+(phase*cycle)-1
                num_all_reads=0
                k=0
                n=0
                m=cycle*2
                pvalue=0

                for i in range(begin,finish+1):
                    try:
                        k+=mapped_data_per_size_per_register[chromosome][register][i]
                    except KeyError:
                        pass

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

                # register_i is an iterator different from register
                for register_i in sorted(mapped_data_per_size_per_register[chromosome]):
                    for i in range(begin,finish+1):
                        try:
                            n+=mapped_data_per_size_per_register[chromosome][register_i][i]
                        except KeyError:
                            pass

                if n/num_all_reads<0.3:
                    begin+=phase
                    continue

                for x in range(k,m+1):
                    numerator=nCr((phase-1)*m,n-x)*nCr(m,x)
                    pvalue+=numerator
                denominator=nCr(phase*m,n)
                pvalue=pvalue/denominator

                if pvalue>=options.pvalue_cutoff:
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
    out_filename=options.output_directory_per_run+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+".positive_phase_loci"
    fhw=open(out_filename,"w")

    for chromosome in sorted(whole_mapped_data):
        filename=options.output_directory_per_run+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+"_"+chromosome+".regionsOfInterest.concentrated"
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
    phased_loci_filename=options.output_directory_per_run+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+".positive_phase_loci"
    final_phase_loci=options.output_directory_per_run+"/"+options.input_filename+"_"+str(phase)+"_"+str(cycle)+".phasing_score_phase_loci"
    fhr=open(phased_loci_filename,"r")
    out4=open(final_phase_loci,"w")

    for line in fhr:
        chromosome,ss,ee=line.strip().split()
        ss=int(ss)
        ee=int(ee)
        phasing_score_filename=options.output_directory_per_run+"/"+str(phase)+"_"+str(chromosome)+"_"+str(ss)+"_"+str(ee)+".phasing_score"
        abundance_score_filename=options.output_directory_per_run+"/"+str(phase)+"_"+str(chromosome)+"_"+str(ss)+"_"+str(ee)+".abundance"
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
            

def cleanUpTemporaryFiles(options):
    """
    Performs cleanup of the directory and keeps only the files that are required.
    """
    os.system("rm "+options.output_directory_per_run+"/*.abundance")
    os.system("rm "+options.output_directory_per_run+"/*.phasing_score")
    os.system("rm "+options.output_directory_per_run+"/*regionsOfInterest*")
    os.system("mv "+options.output_directory_per_run+"/* "+options.output_directory_per_run+"/../")
    os.system("rm -rf "+options.output_directory_per_run)


def main():
    log_filepath = create_logs()

    commandLineArg=sys.argv
    if len(commandLineArg)==1:
        print("Please use the --help option to get usage information")

    logging.info("Parsing arguments...")
    options = parseCommandLineArguments()
    logging.info("Arguments passed: %s", options)

    logging.info("Checking arguments...")
    options = analyzeCommandLineArguments(options)
    logging.info("Checked arguments: %s", options)

    logging.info("Trimming adapters...")
    trimAdapters(options, log_filepath)

    logging.info("Consolidating fastq reads...")
    consolidateReads(options)

    logging.info("Mapping reads to genome using bowtie...")
    mapSmallRNAReadsToGenomeUsingBowtie1(options)

    for phase in options.small_rna_size:

        logging.info("Reading mapped data for phase length %s...", phase)
        whole_mapped_data, mapped_data_per_size_per_register, score_dict, readcount_dict = readMappedData(options,phase)

        for cycle in options.number_of_cycles:
            options.output_directory_per_run=options.output_directory+"/"+"phase_"+str(phase)+"_cycle_"+str(cycle)
            cmd="mkdir "+ options.output_directory_per_run
            os.system(cmd)

            logging.info("Sifting regions of interest...")
            siftRegionsOfInterest(options,mapped_data_per_size_per_register,phase,cycle)

            logging.info("Filtering results by region")
            filterByRegion(options,mapped_data_per_size_per_register, phase, cycle)

            logging.info("Calculating p values...")
            computePValues(options,whole_mapped_data,mapped_data_per_size_per_register,phase,cycle)

            logging.info("Generating phased loci...")
            generatePositivePHASLoci(options,whole_mapped_data,phase,cycle)

            logging.info("Generating phasing scores...")
            generatePhasingScore(options, phase, cycle, score_dict, readcount_dict)

            logging.info("Creating plots...")
            cmd="Rscript --vanilla plot.R "
            cmd+=" "+options.output_directory_per_run
            cmd+=" "+str(phase)
            cmd+=" "+str(cycle)
            logging.info("Calling Rscript plot.R with the command: %s", cmd)
            result = subprocess.run([cmd], shell=True, capture_output=True, text=True, check=True)
            logging.info("R script stdout:\n %s", result.stdout)
            if result.stderr:
                logging.warning("R script stderr:\n %s", result.stderr)

            if options.clean_up!=0:
                logging.info("Clean up prompted, removing temporary files...")
                cleanUpTemporaryFiles(options)
    

if __name__ == "__main__":
    main()
    logging.info("Script complete!")
    
    
    
    
    
    
    
    
    
    
    
    
    
    