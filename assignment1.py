import mysql.connector
import pysam

__author__ = 'ARNOLD Vivienne'

##
## Concept:
## TODO
##


class Assignment1:
    
    def __init__(self):
        ## gene of interest
        self.gene = "KCNE1"

        ## open bamfile as samfile
        self.samfile = pysam.AlignmentFile("chr21.bam", "rb")

    
    def download_gene_coordinates(self, genome_reference, file_name):
        ## TODO concept
        
        print("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
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

        gene_attributes = []
        separator = '\t'

        ## Write to file
        ## TODO this may need some work 
        with open(file_name, "w") as fh:
            for row in cursor:
                if row[0] == self.gene:
                    current_match = []
                    for query_field in row:
                        #fh.write(str(query_field) + "\n")
                        current_match.append(str(query_field))
                    gene_attributes.append(current_match)
                    fh.write(separator.join(current_match)+"\n")

        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print("Done fetching data")
        print(gene_attributes)
        return gene_attributes
        
    def get_coordinates_of_gene(self):

        print("todo")
        
    def get_gene_symbol(self):
        print("todo")
                        
    def get_sam_header(self):
        print("todo")
        
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
        print("Print all results here")
    
    
def main():
    print("Assignment 1")
    assignment1 = Assignment1()
    gene_entry = assignment1.download_gene_coordinates("hg38", "output")[0]
    print(gene_entry)
    assignment1.print_summary()
    
    
    print("Done with assignment 1")
    
        
if __name__ == '__main__':
    main()
    
    
