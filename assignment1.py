import mysql.connector
import pysam
import os

__author__ = 'ARNOLD Vivienne'

##
## Concept:
## TODO
##


class Assignment1:
    
    def __init__(self, gene_of_interest, genome_reference, bam_file, output_file):

        self.gene = gene_of_interest

        self.genome_reference = genome_reference

        self.bam_file = bam_file

        # open bamfile as samfile
        self.samfile = pysam.AlignmentFile(bam_file, "rb")

        self.output_file = output_file

        if os.path.exists(output_file) and os.path.isfile(output_file):
            print("Fetch data from file <"+output_file+">")

            file = open(output_file, "r")
            self.transcript_info = file.readline().strip("\n").split("\t")

        else:
            self.transcript_info = self.download_gene_coordinates()

        self.gene_symbol = self.transcript_info[0]

        self.transcript = self.transcript_info[1]

        self.chromosome = self.transcript_info[2][3:]

        self.start_position = int(self.transcript_info[3])

        self.stop_position = int(self.transcript_info[4])

        self.strand = self.transcript_info[5]

        self.number_of_exons = int(self.transcript_info[6])

        self.exon_coordinates = self.get_exon_coordinates()

    def download_gene_coordinates(self):

        print("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password',
                                      db=self.genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
        
        ## Execute query
        cursor.execute(query)

        attributes_of_all_transcripts = []
        separator = '\t'

        for row in cursor:

            if row[0] == self.gene:

                transcript_attributes = []

                for query_field in row:
                    transcript_attributes.append(query_field)

                attributes_of_all_transcripts.append(transcript_attributes)

        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print("Done fetching data")

        transcript_attributes = attributes_of_all_transcripts[0]

        ## Write to file
        with open(self.output_file, "w") as fh:
            fh.write(separator.join(str(attribute) for attribute in transcript_attributes) + "\n")

        return transcript_attributes

    def print_gene_symbol(self):

        print("Genome reference:".ljust(20, " ")+self.genome_reference)

        print("Gene symbol:".ljust(20, " ")+self.gene_symbol)

        print("Transcript:".ljust(20, " ")+self.transcript)

    def print_coordinates_of_gene(self):

        string_chromosome = "Chromosome "+ self.chromosome

        string_start_position = format(self.start_position, "08,d")

        string_stop_position = format(self.stop_position, "08,d")

        if self.strand == "-":
            string_strand = "reverse strand"
        elif self.strand == "+":
            string_strand = "forward strand"
        else:
            string_strand = "unknown orientation"

        print(("Location:").ljust(20, " ")+
                  string_chromosome+": "+string_start_position+"-"+string_stop_position+" "+string_strand)

    def print_number_of_exons(self):

        print("Number of exons:".ljust(20, " ")+str(self.number_of_exons))

    def get_exon_coordinates(self):

        start_of_exons = self.transcript_info[7].strip("b\'").strip(",\'").split(",")

        stop_of_exons = self.transcript_info[8].strip("b\'").strip(",\'").split(",")

        list_of_exon_coordinates = []

        for exon in range(self.number_of_exons):
            list_of_exon_coordinates.append([int(start_of_exons[exon]), int(stop_of_exons[exon])])

        return list_of_exon_coordinates

    def print_region_of_gene(self):

        print("Coordinates:".ljust(20, " ")+"Exon"+" "*6+"Start".ljust(15, " ")+"End")

        for exon in range(self.number_of_exons):
            print(" "*20+(str(exon+1)).ljust(10, " ")+
                  format(self.exon_coordinates[exon][0], "08,d").ljust(15, " ")+
                  format(self.exon_coordinates[exon][1], "08,d"))

    def print_sam_header(self):

        header = self.samfile.header["HD"]

        headerline = ""

        for key in header:
            headerline += key + ": " + header[key] + "\t"

        print("Samfile header:".ljust(20, " ") + headerline)

    def get_properly_paired_reads_of_gene(self):
        print("todo")
        
    def get_gene_reads_with_indels(self):
        print("todo")
        
    def calculate_total_average_coverage(self):
        print("todo")
        
    def calculate_gene_average_coverage(self):
        print("todo")
        
    def get_number_mapped_reads(self):
        print("todo")

    def print_summary(self):

        separate_blocks = "\n"+"="*80+"\n"

        print(separate_blocks)

        print("Information to selected gene\n")

        self.print_gene_symbol()

        self.print_coordinates_of_gene()

        print()

        self.print_number_of_exons()

        self.print_region_of_gene()

        print(separate_blocks)

        print("Information to selected bam-file <"+self.bam_file+">\n")

        self.print_sam_header()

        print(separate_blocks)


    
    
def main():
    print("Assignment 1\n")

    assignment1 = Assignment1("KCNE1", "hg38", "chr21.bam", "KCNE1_transcripts.txt")

    assignment1.print_summary()

    print("Done with assignment 1")

    assignment1.get_exon_coordinates()

if __name__ == '__main__':
    main()
    
    
