import mysql.connector
import pysam
import os

__author__ = 'ARNOLD Vivienne'

##
## Concept:
## Init class Assignment1 with name of gene, genome reference, a bam file and a filename for output.
## Information to gene (transcript, location and gene region) will be shown, as well as a summary of the bam file
## regarding coverage and reads.


class Assignment1:
    
    def __init__(self, gene_of_interest, genome_reference, bam_file, output_file):

        ## Define parameters from function call
        self.gene = gene_of_interest
        self.genome_reference = genome_reference
        self.bam_file = bam_file
        self.alignment_file = pysam.AlignmentFile(bam_file, "rb")
        self.output_file = output_file

        ## Get gene information from data base or existing file
        if os.path.exists(output_file) and os.path.isfile(output_file):
            print("Fetch data from file <"+output_file+">")
            file = open(output_file, "r")
            self.transcript_info = file.readline().strip("\n").split("\t")
        else:
            self.transcript_info = self.download_gene_coordinates()

        print("\nPreparing report of <"+str(bam_file)+">")

        ## Set variables for gene information
        self.gene_symbol = self.transcript_info[0]
        self.transcript = self.transcript_info[1]
        self.chromosome = self.transcript_info[2]
        self.start_position = int(self.transcript_info[3])
        self.stop_position = int(self.transcript_info[4])
        self.strand = self.transcript_info[5]
        self.number_of_exons = int(self.transcript_info[6])
        self.exon_coordinates = self.get_exon_coordinates()

        ## Set variables retrieved from functions
        self.sam_header = self.get_sam_header()
        self.number_of_properly_paired_reads = self.get_number_of_properly_paired_reads_of_gene()
        self.number_of_gene_reads_with_indels = self.get_number_of_gene_reads_with_indels()
        self.mean_gene_coverage = self.get_average_gene_coverage()
        self.mean_total_coverage = self.get_total_average_coverage()
        self.number_of_mapped_reads = self.get_number_mapped_reads()

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

        ## Get query entries for selected gene
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
        string_chromosome = "Chromosome "+self.chromosome[3:]
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

    def get_exon_coordinates(self):
        start_of_exons = str(self.transcript_info[7]).strip("b\'").strip(",\'").split(",")
        stop_of_exons = str(self.transcript_info[8]).strip("b\'").strip(",\'").split(",")
        list_of_exon_coordinates = []

        for exon in range(self.number_of_exons):
            list_of_exon_coordinates.append([int(start_of_exons[exon]), int(stop_of_exons[exon])])

        return list_of_exon_coordinates

    def print_exon_information(self):
        print("Number of exons:".ljust(20, " ")+str(self.number_of_exons))

        print("Coordinates:".ljust(20, " ")+"Exon"+" "*6+"Start".ljust(15, " ")+"End")

        for exon in range(self.number_of_exons):
            print(" "*20+(str(exon+1)).ljust(10, " ")+
                  format(self.exon_coordinates[exon][0], "08,d").ljust(15, " ")+
                  format(self.exon_coordinates[exon][1], "08,d"))

    def get_sam_header(self):
        header = self.alignment_file.header["HD"]

        headerline = ""

        for key in header:
            headerline += key + ": " + header[key] + "\t"

        return headerline

    def get_number_of_properly_paired_reads_of_gene(self):
        counter_properly_paired_reads = 0

        for read in self.alignment_file.fetch(self.chromosome, self.start_position, self.stop_position):
            if read.is_proper_pair:
                counter_properly_paired_reads += 1

        return counter_properly_paired_reads
        
    def get_number_of_gene_reads_with_indels(self):
        counter_reads_with_indels = 0

        for pileupcolumn in self.alignment_file.pileup(self.chromosome, self.start_position, self.stop_position):
            for pileupread in pileupcolumn.pileups:
                if pileupread.indel:
                   counter_reads_with_indels +=1

        return counter_reads_with_indels

    def get_total_average_coverage(self):
        coverage_sum = 0
        counter_column = 0

        for pileupcolumn in self.alignment_file.pileup(self.chromosome):
            coverage_sum += pileupcolumn.n
            counter_column += 1

        average_total_coverage = round((coverage_sum / counter_column), 1)

        return average_total_coverage

    def get_average_gene_coverage(self):
        coverage_sum = 0
        counter_column = 0

        for pileupcolumn in self.alignment_file.pileup(self.chromosome, self.start_position, self.stop_position):
            coverage_sum += pileupcolumn.n
            counter_column += 1

        average_gene_coverage = round((coverage_sum / counter_column), 1)

        return average_gene_coverage

    def get_number_mapped_reads(self):
        counter_mapped_reads = 0

        for read in self.alignment_file.fetch(self.chromosome, self.start_position, self.stop_position):
            if not read.is_unmapped:
                counter_mapped_reads += 1

        return counter_mapped_reads

    def print_summary(self):
        separate_blocks = "\n"+"="*80+"\n"
        print(separate_blocks)

        print("Information to selected gene\n")
        self.print_gene_symbol()
        self.print_coordinates_of_gene()
        print()
        self.print_exon_information()
        print(separate_blocks)

        print("Bamfile information\n")
        print("Samfile header:".ljust(25, " ")+self.sam_header)
        print("Mean total coverage:".ljust(25, " ")+str(self.mean_total_coverage)+" %")
        print()
        print("Information to selected gene")
        print("Mapped reads:".ljust(25, " ")+str(self.number_of_mapped_reads))
        print("Properly paired reads:".ljust(25, " ")+str(self.number_of_properly_paired_reads))
        print("Reads with indels:".ljust(25, " ")+str(self.number_of_gene_reads_with_indels))
        print("Mean gene coverage:".ljust(25, " ")+str(self.mean_gene_coverage)+" %")
        print(separate_blocks)


def main():
    print("Assignment 1\n")
    assignment1 = Assignment1("KCNE1", "hg38", "chr21.bam", "KCNE1_transcripts.txt")
    assignment1.print_summary()
    print("Done with assignment 1")
    assignment1.alignment_file.close()


if __name__ == '__main__':
    main()
