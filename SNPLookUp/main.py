from query_objects import *
from visualisation import *
from exceptions import *
from helper import *


def query(query_file, exclude_regions=True, version="hg19", save_to_txt=False, dir_name="raw_query_output"):
    print("Initialising...", end="\t")
    chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                   "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
    # processing queries
    query_regions = parse_query_file(query_file)

    # if exclude regions, read exclude regions file
    exclude_regions_by_chromosome = initialise_exclude_regions(exclude_regions)
    print("Done!")
    # arrange queryregions by chromosome, create QueryRegion objects (and SubRegion objects)
    # after filtering out excluded regions. SNPs have not been added at this point
    print("Creating QueryRegions...", end="\t")
    queryregions_by_chromosome = initialise_query_regions(query_regions, exclude_regions_by_chromosome)
    print("Done!")
    # Retrieving and adding SNPs
    print("Retrieving SNPs...")
    for chrom in queryregions_by_chromosome:
        retrieve_n_add_snps(chrom, queryregions_by_chromosome[chrom], version)
    # Writing data to file, if save_to_txt is True
    if save_to_txt:
        print("Saving query to text...", end="\t")
        save(queryregions_by_chromosome, dir_name)
        print("Saved!")
    # returns the list of QueryRegion objects with all qualifying SNPs, unfiltered
    print("Done!")
    return queryregions_by_chromosome


def snp_filter(query_region, ldmatrix_access_token, gmaf, r2, max_number_of_snps, save_to_txt=False,
               dirname="filter_output_all"):
    print("Filtering snps...")
    query_region_gmaf = snp_filter_by_gmaf(query_region, gmaf)
    query_region_gmaf_ld = snp_filter_by_ld(query_region_gmaf, ldmatrix_access_token, r2)
    query_region_final = snp_filter_by_max_number(query_region_gmaf_ld, max_number_of_snps)
    if save_to_txt:
        save(query_region_final, dirname)
    return query_region_final


def snp_filter_by_gmaf(query_region, gmaf, save_to_txt=False, dirname="filter_output_gmaf"):
    print(f"Filtering by gmaf of {gmaf}...", end="\t")
    query_region = query_clone(query_region)
    if type(query_region) == dict:
        temp = []
        for qreg in query_region.values():
            temp.extend(qreg)
        query_region = temp
    for qreg in query_region:
        qreg_snps = qreg.get_snps()
        snps_to_keep = [snp for snp in qreg_snps if snp.gmaf >= gmaf]
        qreg.keep_snps(snps_to_keep)
    if save_to_txt:
        save(query_region, dirname)
    print("Done!")
    return query_region


def snp_filter_by_ld(query_region, ldmatrix_access_token, r2, save_to_txt=False, dirname="filter_output_ld"):
    print(f"Filtering by LD of {r2}...")
    query_region = query_clone(query_region)
    if "Y" in query_region:
        print("Chromosome Y cannot be filtered by LD. Removing...")
        del query_region["Y"]
    if type(query_region) == dict:
        temp = []
        for qreg in query_region.values():
            temp.extend(qreg)
        query_region = temp
    for qreg in query_region:
        qreg_snps = qreg.get_snps()
        ldmatrix = get_ldmatrix(qreg_snps, ldmatrix_access_token, qreg.get_label())
        qreg.keep_snps(ld_filter(qreg_snps, ldmatrix, r2))
    if save_to_txt:
        save(query_region, dirname)
    print("Done!")
    return query_region


def snp_filter_by_max_number(query_region, max_number_of_snps, save_to_txt=False, dirname="filter_output_maxnum"):
    print(f"Filtering by max {max_number_of_snps} number of snps...", end="\t")
    query_region = query_clone(query_region)
    if type(query_region) == dict:
        temp = []
        for qreg in query_region.values():
            temp.extend(qreg)
        query_region = temp
    for qreg in query_region:
        qreg_snps = sorted(qreg.get_snps(), key=lambda x: x.pos)
        if len(qreg_snps) > max_number_of_snps:
            qreg.keep_snps(even_selection(qreg_snps, max_number_of_snps))
    if save_to_txt:
        save(query_region, dirname)
    print("Done!")
    return query_region


def save_for_design_studio(queryregions_by_chromosome, file_name="input", dir_name="input_for_design_studio"):
    out_path = os.path.join("..", "output", dir_name)
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    if type(queryregions_by_chromosome) == list:
        queryregions_by_chromosome = {queryregions_by_chromosome[0].chromosome[3:]: queryregions_by_chromosome}
    snps = retrieve_n_sort_snps(queryregions_by_chromosome)
    print("Saving query as Design Studio input...", end="\t")
    save_as_csv(snps, os.path.join(out_path, file_name + ".csv"))
    save_as_bed(snps, os.path.join(out_path, file_name + ".bed"))
    print(f"Done!\nQueries saved in {out_path}")
