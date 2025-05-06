"""
Data input example:
        read_12124_6138 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTT   IIIIIIIIIIIIIIIIIIIII   0
        read_3307062_86 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTAGT        IIIIIIIIIIIIIIIIIIIIIIII        0
        read_4610579_11 -       chr05   3810328 TTTGAGCAAGAAAGTCAGAGTT  IIIIIIIIIIIIIIIIIIIIII  0
        read_6964633_31 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTA  IIIIIIIIIIIIIIIIIIIIII  0
        read_7421569_42 +       chr05   3810329 TTGAGCAAGAAAGTCAGAGTTAG IIIIIIIIIIIIIIIIIIIIIII 0


# Computes the P-values using a Hypergeometric distribution
M = total number of de-duped reads (coordinates) in the fastq that is phase length (21)
n = number of de-duped reads that cover the phased window/phasiRNA region (phase length 21)
N = number of all possible phasing positions per that phasiRNA region (region // phase lenth (21); floor division)
k = number of phased positions that are covered by a read

pval_input = {
        'Pathway_A': (20000, 500, 400, 25),  # 25 genes out of 500 genes are in pathway with 400 genes
        'Pathway_B': (20000, 500, 150, 15),  # 15 genes out of 500 genes are in pathway with 150 genes
        'Pathway_C': (20000, 500, 300, 10),  # 10 genes out of 500 genes are in pathway with 300 genes
        'Pathway_D': (20000, 500, 200, 20),  # 20 genes out of 500 genes are in pathway with 200 genes
        'Pathway_E': (20000, 500, 250, 18),  # 18 genes out of 500 genes are in pathway with 250 genes
        'Pathway_F': (20000, 500, 100, 3),   # 3 genes out of 500 genes are in pathway with 100 genes
        'Pathway_G': (20000, 500, 50, 8),    # 8 genes out of 500 genes are in pathway with 50 genes
        'Pathway_H': (20000, 500, 80, 5),    # 5 genes out of 500 genes are in pathway with 80 genes
    }

process:

take alignment
check sense/antisense
count reads of length 21 for M
filter for region
count reads of length 21 for n

for each register:
    calculate N given length of region (tas4) - bear in mind cut site (trigger +10) and a new window based on each register 1bp, 2bp, 3bp etc
    calculate k - starting at "start", how many of N are covered by a read. +21 loop until end of range, keep counter of reads.

"""
def check_in_region(chr, coord, strand, regions_dict):

    print(chr)
    print(coord)
    print(regions_dict["start"], regions_dict["end"])


    if chr != regions_dict["chr"]:
        print(f"{chr} not in line")
        return False
    
    # bin together antisense and sense strands
    # TODO as this doesnt work properly. look at antisense and sense in relation to gene region
    if strand == "-":
        coord += 2
    
    if regions_dict["start"] <= coord <= regions_dict["end"]:
        print("pass")
        return True

    return False


def pval_prepare(alignment_file, regions_dict, mirna_length):

    # pval_M
    total_coords_covered = 0
    # pval_n
    region_coords_covered = 0
    # register based pval_N
    register_phased_pos_count = {}
    # register based pval_k
    register_phased_coords_covered = {}

    pval_data = {}
    
    with open(alignment_file, "r") as infile:
        for line in infile:

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

            # only want reads that are of the right length, otherwise jump to next line in file
            if length != mirna_length:
                continue

            total_coords_covered +=1

            # only keep reads/lines in ROI
            if check_in_region(chromosome, coordinate, strand, regions_dict):
                region_coords_covered +=1

                # calculate register based off of remainder from cut site
                register = (coordinate - regions_dict["trigger_cut_site"]) % mirna_length
                print(f"register is {register}")

                # bin coords covered per register into dict #TODO allow for double because of antisense?
                if register in register_phased_coords_covered.keys():
                    register_phased_coords_covered[register] +=1
                else:
                    register_phased_coords_covered[register] = 1

                # calculate register possible phased loc #TODO allow for double because of antisense?
                if register not in register_phased_pos_count.keys():
                    register_phased_pos_count[register] = (regions_dict["end"] - (regions_dict["trigger_cut_site"] + register)) // mirna_length

    # iterate through 
    for reg_pos in range(0, mirna_length):
        print(f"trying to get data for register position: {reg_pos}")
        try:
            pval_data[reg_pos] = [
                total_coords_covered,
                region_coords_covered,
                register_phased_pos_count[reg_pos],
                register_phased_coords_covered[reg_pos]
            ]
        except:
            print(f"No data for register {reg_pos}")

    print(pval_data)
    

            


if __name__ == "__main__":

    """
    region:
    chr05	slyTAS4	tasRNA	3810129	3810580	.	+	.

    trigger:
    chr05	3810317	3810339	slyTAS4_sly_miR828	0	-
    cutsite = trigger start +10
    """

    mirna_length = 21
    regions_dict = {
        "chr": "chr05",
        "start":3810129,
        "end": 3810580,
        "strand": "+",
        "trigger_cut_site": 3810327
    }
    alignment_file = "/home/kal_grimes_tropic_bio/proj/phasi/findPhasiRNAs/results/wtm82_fastq/wtm82_bowtie1.bwt"
    # alignment_file = "/home/kal_grimes_tropic_bio/proj/phasi/test_bowtie.bwt"
    # regions_dict = {
    #     "chr": "chr07",
    #     "start":2013338,
    #     "end": 2013450,
    #     "strand": "+",
    #     "trigger_cut_site": 2013349
    # }

    pval_prepare(alignment_file, regions_dict, mirna_length)