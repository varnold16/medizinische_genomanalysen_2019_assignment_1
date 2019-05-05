import mysql.connector
import pysam

__author__ = 'ARNOLD Vivienne'

##
## Concept:
## TODO
##


class Assignment1:
    
    def __init__(self, gene_of_interest, genome_reference, bam_file, output_file):
        ## gene of interest
        self.gene = gene_of_interest

        self.genome_reference = genome_reference

        ## open bamfile as samfile
        self.samfile = pysam.AlignmentFile(bam_file, "rb")

        self.output_file = output_file

        self.transcript_info = self.download_gene_coordinates()

    
    def download_gene_coordinates(self):
        ## TODO concept
        
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

        ## Write to file
        ## TODO this may need some work 
        with open(self.output_file, "w") as fh:

            for row in cursor:

                if row[0] == self.gene:

                    transcript_attributes = []

                    for query_field in row:
                        transcript_attributes.append(query_field)

                    attributes_of_all_transcripts.append(transcript_attributes)


            for transcript in attributes_of_all_transcripts:
                fh.write(separator.join(str(attribute) for attribute in transcript)+"\n")

        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print("Done fetching data")

        self.transcript_attributes = attributes_of_all_transcripts[0]

        return transcript_attributes

        
    def get_coordinates_of_gene(self):

        chromosome = "Chromosome "+ self.transcript_info[2][3:]

        start_position = format(self.transcript_info[3], "08,d")

        stop_position = format(self.transcript_info[4], "08,d")

        if self.transcript_info[5] == "-":
            strand = "reverse strand"
        elif self.transcript_info[5] == "+":
            strand = "forward strand"
        else:
            strand = "unknown orientation"

        print(("Location:").ljust(20, " ")+
              chromosome+": "+start_position+"-"+stop_position+" "+strand)

        
    def get_gene_symbol(self):

        gene_symbol = self.transcript_attributes[0]

        transcript = self.transcript_attributes[1]

        print("Genome reference:".ljust(20, " ")+self.genome_reference)
        print("Gene symbol:".ljust(20, " ")+gene_symbol)
        print("Transcript:".ljust(20, " ")+transcript)

                        
    def get_sam_header(self):

        header = self.samfile.header["HD"]

        headerline = ""

        for key in header:
            headerline += key+": "+header[key]+"\t"

        print("Header:".ljust(20, " ")+headerline)

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

    def get_region_of_gene(self):
        print("todo")
        
    def get_number_of_exons(self):
        print("ads")
    
    
    def print_summary(self):
        self.get_gene_symbol()

        self.get_coordinates_of_gene()

        self.get_sam_header()



        print("Print all results here")
    
    
def main():
    print("Assignment 1")
    assignment1 = Assignment1("KCNE1", "hg38", "chr21.bam", "KCNE1_transcripts.txt")
    print(assignment1.transcript_info)
    assignment1.print_summary()
    
    
    print("Done with assignment 1")
    
        
if __name__ == '__main__':
    main()
    
    
