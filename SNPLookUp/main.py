"""Finds single nucleotide polymorphism (SNP) loci within an input region of the human genome.

SNPLookUp takes as input a bed file containing regions of interest in the
GRCh37 chromosome. For each region of interest, SNPs within the region are
identified and filtered based on user-defined cutoffs for minimum GMAF,
maximum pairwise LD between SNPs within the region, and/or maximum number of
SNPs. SNPs selected can be subsequently visualised and saved in either csv or
bed format that are compatible with Illumina's DesignStudio web service.

  Typical usage example:

  test_query = query(input_file_path, save_to_txt=True)
  filtered_test_query = snp_filter(test_query, LdMatrix_access_token, 0.4, 0.2, 50, save_to_txt=True)
  save_for_design_studio(filtered_test_query, file_name="test_query")
"""
from query_objects import *
from visualisation import *
from exceptions import *
from helper import *


def query(query_file,
          exclude_regions=True,
          version="hg19",
          save_to_txt=False,
          dir_name="raw_query_output"):
    """Runs the initial query of regions in the query bed file.

    Converts each input region of interest into a QueryRegion object, which
    contains information about valid sub-regions of interest and SNPs in the
    region.

    Args:
        query_file: Either a string containing the path to a .bed query file,
            or a list in the format of:
             [ ["chr1", start_pos, end_pos, label],
               ["chr5", start_pos, end_pos, label], ... ].
        exclude_regions: Boolean. If True, regions listed in .bed files in the
            /excluded_regions folder will be excluded from the query.
            (default: True)
        version: String. States the genome version of the input chromosomal
            positions. (default: "hg19")
        save_to_txt: Boolean. If True, saves each query as a text file in
            output/dir_name/{query.label}.txt. (default: False)
        dir_name: String that sets the name of the directory to save the
            output to. (default: "raw_query_output")

    Returns:
        A dictionary in the format of chromosome_number:list of QueryRegions.
        For example:

        {"1": [<QueryRegion>, <QueryRegion>],
         "2": [<QueryRegion>, <QueryRegion>, <QueryRegion>],
         "X": [<QueryRegion>]}

    Raises:
        ParseException: Invalid input in query_file.
    """
    print("Initialising...", end="\t")
    chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                   "12", "13", "14", "15", "16", "17", "18", "19", "20", "21",
                   "22", "X", "Y")
    # processing input queries
    query_regions = parse_query_file(query_file)

    # if exclude regions, read exclude regions file
    exclude_regions_by_chromosome = initialise_exclude_regions(exclude_regions)
    print("Done!")
    # arrange queryregions by chromosome, create QueryRegion objects (and SubRegion objects)
    # after filtering out excluded regions. SNPs have not been added at this point
    print("Creating QueryRegions...", end="\t")
    queryregions_by_chromosome = initialise_query_regions(
        query_regions, exclude_regions_by_chromosome)
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


def snp_filter(query_region,
               ldmatrix_access_token,
               gmaf,
               r2,
               max_number_of_snps,
               save_to_txt=False,
               dirname="filter_output_all"):
    """Filters the SNPs in the input QueryRegions.

    This is a composite function that filters for gmaf, ld and max_number.
    Function composes of snp_filter_by_gmaf(), snp_filter_by_ld(), and
    snp_filter_by_maxnum().

    Args:
        query_region. Either a list or dictionary of QueryRegion objects.
        ldmatrix_access_token. String token required to access the ldmatrix api
            via requests
        gmaf. Float between 0.3 - 0.5. Minimum gmaf that SNPs in the
            QueryRegions must have.
        r2. Float between 0 and 1. Maximum linkage disequilibrium measured in
            terms of r^2 value between any two SNPs in each Query Region.
        max_number_of_snps. Int capping the maximum number of SNPs in each
            region.
        save_to_txt. Boolean. If True, saves each query as a text file in
            output/dir_name/{query.label}.txt. (default: False)
        dir_name. String that sets the name of the directory to save the output
            to. (default: "filter_output_all")

    Returns:
        A dictionary in the format of chromosome_number:list of QueryRegions.
        For example:

        {"1": [<QueryRegion>, <QueryRegion>],
         "2": [<QueryRegion>, <QueryRegion>, <QueryRegion>],
         "X": [<QueryRegion>]}
    """
    print("Filtering snps...")
    query_region_gmaf = snp_filter_by_gmaf(query_region, gmaf)
    query_region_gmaf_ld = snp_filter_by_ld(query_region_gmaf,
                                            ldmatrix_access_token, r2)
    query_region_final = snp_filter_by_max_number(query_region_gmaf_ld,
                                                  max_number_of_snps)
    if save_to_txt:
        save(query_region_final, dirname)
    return query_region_final


def snp_filter_by_gmaf(query_region,
                       gmaf,
                       save_to_txt=False,
                       dirname="filter_output_gmaf"):
    """Filters the SNPs in the input QueryRegions based on the gmaf.

    Filters the SNPs found within each input QueryRegion object by a global
    minor allele frequency cutoff defined by the user.

    Args:
        query_region: Either a list or dictionary of QueryRegion objects.
        gmaf: Float between 0.3 - 0.5. Minimum gmaf that SNPs in the
            QueryRegions must have.
        save_to_txt: Boolean. If True, saves each query as a text file in
            output/dir_name/{query.label}.txt. (default: False)
        dirname: String that sets the name of the directory to save the output
            to. (default: "filter_output_gmaf")

    Returns:
        A dictionary in the format of chromosome_number:list of QueryRegions
        consisting of all the input QueryRegions after filtering the list of
        SNPs they contain. For example:

        {"1": [<QueryRegion>, <QueryRegion>],
         "2": [<QueryRegion>, <QueryRegion>, <QueryRegion>],
         "X": [<QueryRegion>]}
    """
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


def snp_filter_by_ld(query_region,
                     ldmatrix_access_token,
                     r2,
                     save_to_txt=False,
                     dirname="filter_output_ld"):
    """Filters the SNPs in the input QueryRegions based on ld.

    Filters the SNPs found within each input QueryRegion object based on their
    pairwise linkage disequilibrium (LD), measured in terms of r^2 value. It is
    guaranteed that after filtering, the LD between any two SNPs in each
    QueryRegion does not exceed the LD cutoff. Pairwise LD data is retrieved
    from LDMatrix web service at https://ldlink.nci.nih.gov/?tab=ldmatrix.

    Args:
        query_region: Either a list or dictionary of QueryRegion objects.
        ldmatrix_access_token: String token required to access the ldmatrix api
            via requests
        r2: Float between 0 and 1. Maximum linkage disequilibrium measured in
            terms of r^2 value between any two SNPs in each Query Region.
        save_to_txt: Boolean. If True, saves each query as a text file in
            output/dir_name/{query.label}.txt. (default: False)
        dirname: String that sets the name of the directory to save the output
            to. (default: "filter_output_ld")

    Returns:
        A dictionary in the format of chromosome_number:list of QueryRegions
        consisting of all the input QueryRegions after filtering the list of
        SNPs they contain. For example:

        {"1": [<QueryRegion>, <QueryRegion>],
         "2": [<QueryRegion>, <QueryRegion>, <QueryRegion>],
         "X": [<QueryRegion>]}
    """
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
        ldmatrix = get_ldmatrix(qreg_snps, ldmatrix_access_token,
                                qreg.get_label())
        qreg.keep_snps(ld_filter(qreg_snps, ldmatrix, r2))
    if save_to_txt:
        save(query_region, dirname)
    print("Done!")
    return query_region


def snp_filter_by_max_number(query_region,
                             max_number_of_snps,
                             save_to_txt=False,
                             dirname="filter_output_maxnum"):
    """Filters the SNPs in the input QueryRegions by a maximum number.

    Filters the SNPs in the input QueryRegions, capping number of SNPs to
    max_number. This program tries to select SNPs evenly in each QueryRegion
    to keep, such that the SNPs remaining are spread apart as evenly as
    possible within the QueryRegion.

    Args:
        query_region: Either a list or dictionary of QueryRegion objects
        max_number_of_snps: Int capping the maximum number of SNPs in each
            region.
        save_to_txt: Boolean. If True, saves each query as a text file in
            output/dir_name/{query.label}.txt. (default: False)
        dirname: String that sets the name of the directory to save the output
            to. (default: "filter_output_maxnum")

    Returns:
        A dictionary in the format of chromosome_number:list of QueryRegions
        consisting of all the input QueryRegions after filtering the list of
        SNPs they contain. For example:

        {"1": [<QueryRegion>, <QueryRegion>],
         "2": [<QueryRegion>, <QueryRegion>, <QueryRegion>],
         "X": [<QueryRegion>]}
    """
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


def save_for_design_studio(queryregions_by_chromosome,
                           file_name="input",
                           dir_name="input_for_design_studio"):
    """Saves QueryRegions and Snps contained.

    Saves the input QueryRegions into a .bed file and a .csv file that qre
    compatible as inputs to Illumina's Design Studio.

    Args:
        queryregions_by_chromosome: Either a list or dictionary of QueryRegion
            objects.
        file_name: String that sets the name of the .bed and .csv files.
            (default: "input")
        dir_name: String that sets the name of the directory to save the
            output to. (default: "input_for_design_studio")
    """
    out_path = os.path.join("..", "output", dir_name)
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    if type(queryregions_by_chromosome) == list:
        queryregions_by_chromosome = {
            queryregions_by_chromosome[0].chromosome[3:]:
            queryregions_by_chromosome
        }
    snps = retrieve_n_sort_snps(queryregions_by_chromosome)
    print("Saving query as Design Studio input...", end="\t")
    save_as_csv(snps, os.path.join(out_path, file_name + ".csv"))
    save_as_bed(snps, os.path.join(out_path, file_name + ".bed"))
    print(f"Done!\nQueries saved in {out_path}")
