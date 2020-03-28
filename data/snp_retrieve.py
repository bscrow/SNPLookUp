from Bio import Entrez
import time

Entrez.email = "boshen@live.com"


def dbsnp_retrieve(chrom):
    """Retrieves SNP infromation from NCBI's SNP database.

    Retrieves all SNPs from each chromosome with GMAF of at least 0.3, and
    then stores the list of SNPs retrieved. This program is run before
    SNPLookUp is used in order to save on the time required to send and
    process requests to NCBI.

    Args:
        chrom: chromosome to retrieve information from.

    """
    output = []
    for a in range(30, 51):
        gmaf_range = f"{a / 100 - 0.005}:{a / 100 + 0.00499}"
        search_term = f"{chrom}[Chromosome] AND" + f"(k[Allele] OR m[Allele] OR r[Allele] OR s[Allele] OR w[Allele] OR y[Allele]) " + f"AND {gmaf_range}[Global Minor Allele Frequency] " + "AND \"snv\"[SNP Class] " + "NOT \"merged rs\"[Filter]"
        handle = Entrez.esearch(db="snp", term=search_term, RetMax=200000)
        record = Entrez.read(handle)
        print(
            f"{record['Count']} matches found on Chromosome {chrom} for gmaf range {gmaf_range}"
        )
        for snp in record["IdList"]:
            output.append("rs" + snp + f"\t{a / 100:.3f}\n")
        time.sleep(3)
    with open(f"snps_by_chr_hg38//chr{chrom}_snps.txt", "w") as file:
        file.writelines(output)
    print(len(output))


if __name__ == "__main__":
    chromosomes = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y".split(
        " ")
    for chromosome in chromosomes:
        dbsnp_retrieve(chromosome)
        time.sleep(5)
