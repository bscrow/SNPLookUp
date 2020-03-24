
from exceptions import *
from query_objects import *

import csv
import glob
import gzip
import os
import requests
import time
from tqdm import tqdm
import urllib3


def even_selection(snps, num):
    if num == len(snps):
        return snps
    elif num == 1:
        return [snps[len(snps) // 2]]

    left, right = [], []
    partition = (snps[0].pos + snps[-1].pos) // 2
    leftn, rightn = num - num // 2, num // 2
    for snp in snps:
        if snp.pos <= partition:
            left.append(snp)
        else:
            right.append(snp)
    if len(left) < leftn:
        leftn = len(left)
        rightn = num - len(left)
    elif len(right) < rightn:
        rightn = len(right)
        leftn = num - len(right)

    return even_selection(left, leftn) + even_selection(right, rightn)


def filter_regions(start_pos, end_pos, exclude_regions):
    subregions = []
    for i in exclude_regions:
        exclude_start, exclude_end = int(i[0]), int(i[1])
        if exclude_start < start_pos <= exclude_end < end_pos:
            start_pos = exclude_end + 1
        elif start_pos < exclude_start <= exclude_end < end_pos:
            subregions.append((start_pos, exclude_start))
            start_pos = exclude_end + 1
        elif start_pos < exclude_start < end_pos < exclude_end:
            subregions.append((start_pos, exclude_start))
            break
    subregions.append((start_pos, end_pos))
    return subregions


def get_ldmatrix(snps, ldmatrix_access_token, label):
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    matrices = []
    headers = {'content-type': 'application/json'}
    params = (('token', ldmatrix_access_token),)
    snps = [snp.rsid for snp in snps]
    count = 0
    total = len(snps)
    print(f"Retrieving ldmatrix for a total of {total} SNPs in {label}...\n")
    while len(snps) > 900:
        snps_query = "\n".join(snps[:900])
        print(f"Retrieving LDMatrix for SNPs {count}-{min(count+900, total)}")
        data = ('{"snps": "' + snps_query + '", "pop": "ALL", "r2_d": "r2"}').encode('unicode_escape')
        response = requests.post('https://ldlink.nci.nih.gov/LDlinkRest/ldmatrix', headers=headers, params=params,
                                 data=data, verify=False, timeout=300)
        if response.ok:
            matrices.append(response.text)
        else:
            raise LdmatrixException(response.text)
        snps = snps[450:]
        count += 450
        time.sleep(5)
    return matrices


def initialise_exclude_regions(is_exclude):
    chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                   "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
    exclude_regions_by_chromosome = {}
    for chromosome in chromosomes:
        exclude_regions_by_chromosome[chromosome] = []

    if is_exclude:
        print("Excluding unwanted regions...")

        exclude_region_files = glob.glob(os.path.join("..", "data", "excluded_regions", "*.bed"))

        for exclude_region_file in exclude_region_files:
            with open(exclude_region_file, "r") as f:
                exclude_data = f.readlines()
            for line in exclude_data:
                entry = line.strip().split("\t")
                exclude_regions_by_chromosome[entry[0][3:]].append(entry[1:])

    for chrom in exclude_regions_by_chromosome:
        exclude_regions_by_chromosome[chrom].sort(key=lambda x: x[0])
    return exclude_regions_by_chromosome


def initialise_query_regions(query_regions, exclude_regions_by_chromosome):
    chromosomes = ("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                   "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
    queryregions_by_chromosome = {}
    for chromosome in chromosomes:
        queryregions_by_chromosome[chromosome] = []

    for region in query_regions:
        if region[0][3:] not in chromosomes:
            raise ParseException(f"{region[0][3:]} is not a valid chromosome\n" +
                                 f"In query {region}")
        chrom, start, end, label = region[0], int(region[1]), int(region[2]), region[3]
        queryregions_by_chromosome[region[0][3:]] \
            .append(QueryRegion(chrom, start, end, label,
                                filter_regions(start, end, exclude_regions_by_chromosome[region[0][3:]])))
    for chrom in chromosomes:
        if not queryregions_by_chromosome[chrom]:
            del queryregions_by_chromosome[chrom]
    return queryregions_by_chromosome


def ld_filter(snps, matrices, r2):
    snps_rsid_to_keep = []  # good snps, to be kept
    keep = [True] * 900
    first = True
    for matrix in matrices:
        if not first:
            keep = keep[450:] + [True] * 450
        first = False
        snp_rsid_list = []
        LDMatrix = []
        lines = matrix.split("\n")
        for line in lines:
            entry = line.split("\t")
            snp_rsid_list.append(entry[0])
            LDMatrix.append([float(a) for a in entry[1:]])
        for i in range(len(LDMatrix)):
            if keep[i]:
                for j in range(i + 1, len(LDMatrix[i])):
                    if LDMatrix[i][j] > r2:
                        keep[j] = False
        for k in range(len(snp_rsid_list)):
            if keep[k]:
                snps_rsid_to_keep.append(snp_rsid_list[k])
    for snp in snps:
        if snp.rsid not in snps_rsid_to_keep:
            snps.remove(snp)
    return snps


def parse_query_file(query_file):
    query_regions = []
    if type(query_file) == str and os.path.isfile(query_file):
        with open(query_file, "r") as datafile:
            data = datafile.readlines()
        for line in data:
            query_regions.append(line.strip().split("\t")[:4])

    elif type(query_file) == list:
        query_regions = query_file
    else:
        raise ParseException("Please ensure that the input file is either a list of queries, " +
                             "or a path string to the query file")
    return query_regions


def query_clone(query_region):
    if type(query_region) == list:
        query_region_clone = [region.clone() for region in query_region]
    elif type(query_region) == dict:
        query_region_clone = dict.fromkeys(query_region.keys())
        for chrom in query_region:
            query_region_clone[chrom] = [region.clone() for region in query_region[chrom]]
    else:
        raise ParseException("function only accepts list or dict objects as input")
    return query_region_clone


def retrieve_n_add_snps(chrom, queryregions_in_chrom, version):
    if version == "hg38":
        raise UnderDevelopmentException("Support for hg38 is still under development.")
    elif version != "hg19":
        raise ParseException(f"Genome version {version} not recognised.")
    print(f"Retrieving SNPs from Chromosome {chrom}...")
    subregions = []
    for queryregion in queryregions_in_chrom:
        subregions.extend(queryregion.get_subregions())

    snps_to_gmaf = {}
    with open(os.path.join("..", "data", "snps_by_chr_hg38", f"chr{chrom}_snps.txt"), "r") as file:
        for entry in file:
            snp, gmaf = entry.strip().split("\t")
            snps_to_gmaf[snp] = float(gmaf)

    with gzip.open(os.path.join("..", "data", "snps_by_chr_hg19", f"chr_{chrom}.txt.gz"), 'rt') as file:
        for a in range(7):  # remove blanks
            file.readline()
        for line in tqdm(file):
            line = line.split("\t")
            # checking if snp entry is accurate
            if not (line[20] == "151" or line[21] == "GRCh37.p13"):
                continue
            try:
                snp_rsid, snp_pos = "rs" + line[0], int(line[11])
            except ValueError:
                continue
            for subregion in subregions:
                if snp_pos in subregion and snp_rsid in snps_to_gmaf:
                    subregion.add_snp(Snp(chrom, snp_pos, snp_rsid, snps_to_gmaf[snp_rsid]))
    print(f"Chromosome {chrom} done!")


def retrieve_n_sort_snps(queryregions_by_chromosome):
    all_snps = []
    snps = []
    for queryregions in queryregions_by_chromosome.values():
        for queryregion in queryregions:
            all_snps.extend(queryregion.get_snps())
    for a_snp in all_snps:
        add = True
        for snp in snps:
            if a_snp.rsid == snp.rsid or (a_snp.chromosome == snp.chromosome and a_snp.pos == snp.pos):
                add = False
                break
        if add:
            snps.append(a_snp)

    snps.sort(key=lambda x: x.pos)
    snps.sort(key=lambda x: x.chromosome)
    return snps


def save(queryregions_by_chromosome, dir_name):
    if not os.path.exists(os.path.join("..", "output", dir_name)):
        os.mkdir(os.path.join("..", "output", dir_name))
    if type(queryregions_by_chromosome) == list:
        queryregions_by_chromosome = {queryregions_by_chromosome[0].chromosome[3:]: queryregions_by_chromosome}
    for chrom in queryregions_by_chromosome:
        if not os.path.exists(os.path.join("..", "output", dir_name, "chr" + chrom)):
            os.mkdir(os.path.join("..", "output", dir_name, "chr" + chrom))
        file_names = {}
        for queryregion in queryregions_by_chromosome[chrom]:
            # handle duplicate region labels
            if queryregion.label in file_names:
                label = queryregion.label + f" ({file_names[queryregion.label]})"
                file_names[queryregion.label] += 1
            else:
                file_names[queryregion.label] = 1
                label = queryregion.label

            with open(os.path.join("..", "output", dir_name, f"chr{chrom}", label + ".txt"), "w") as file:
                file.write(str(queryregion))


def save_as_bed(snps, file_path):
    with open(file_path, "w") as bedfile:
        for snp in snps:
            bedfile.write(f"chr{snp.chromosome}\t{snp.pos-1}\t{snp.pos}\t{snp.rsid}\n")


def save_as_csv(snps, file_path):
    with open(file_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Type", "Name", "Chromosome", "Start", "End"])
        for snp in snps:
            writer.writerow(["Snp", snp.rsid])
