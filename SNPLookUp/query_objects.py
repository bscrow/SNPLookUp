from exceptions import *


class QueryRegion:

    def __init__(self, chromosome, start_pos, end_pos, label, regions, version="hg19"):
        self.chromosome = chromosome
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.label = label
        self.version = version
        if type(regions[0]) == tuple:
            self.regions = []
            for region in regions:
                self.regions.append(SubRegion(chromosome, region[0], region[1]))
        elif isinstance(regions[0], SubRegion):
            self.regions = regions
        else:
            raise ParseException("Incorrect format of SubRegions list")

    def __str__(self):
        ret_str = f"\nRegion in chromosome {self.chromosome}, between positions {self.start_pos}-{self.end_pos}\n " \
                  + f"Region contains the following {len(self.regions)} sub-regions: \n\n"
        for region in self.regions:
            ret_str += str(region) + "\n"
        return ret_str

    def add_snps(self, snps):
        for region in self.regions:
            region.add_snps(snps)

    def clone(self):
        region_clones = [subreg.clone() for subreg in self.regions]
        return QueryRegion(self.chromosome, self.start_pos, self.end_pos, self.label, region_clones, self.version)

    def get_label(self):
        return self.label

    def get_subregions(self):
        return self.regions

    def get_snps(self):
        snps = []
        for region in self.regions:
            snps.extend(region.get_snps())
        return snps

    def keep_snps(self, snps):
        for region in self.regions:
            region.keep_snps(snps)

    def remove_snps(self, snps):
        for region in self.regions:
            region.remove_snps(snps)


class SubRegion:

    def __init__(self, chromosome, start_pos, end_pos):
        self.chromosome = chromosome
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.snps = []

    def __contains__(self, pos):
        return self.start_pos <= pos <= self.end_pos

    def __str__(self):
        ret_str = f"Sub-region in chromosome {self.chromosome}, between positions {self.start_pos}-{self.end_pos}\n " \
                  + f"Sub-region contains the following {len(self.snps)} SNPs: \n\n"
        for snp in self.snps:
            ret_str += str(snp)
        return ret_str

    def add_snp(self, snp):
        assert self.end_pos >= snp.pos >= self.start_pos
        self.snps.append(snp)

    def add_snps(self, snps):
        for snp in snps:
            if self.end_pos >= snp.pos >= self.start_pos:
                self.snps.append(snp)

    def clear_snps(self):
        self.snps = []

    def clone(self):
        snp_clones = [snp.clone() for snp in self.snps]
        subregion_clone = SubRegion(self.chromosome, self.start_pos, self.end_pos)
        subregion_clone.add_snps(snp_clones)
        return subregion_clone

    def get_snps(self):
        return self.snps

    def keep_snps(self, snps):
        self.clear_snps()
        self.add_snps(snps)

    def remove_snps(self, snps_to_remove):
        for snp in self.snps:
            if snp in snps_to_remove:
                self.snps.remove(snp)


class Snp:

    def __init__(self, chromosome, pos, rsid, gmaf, alleles="/"):  # alleles retrieval not yet implemented
        self.chromosome = chromosome
        self.pos = pos
        self.rsid = rsid
        self.gmaf = gmaf
        self.alleles = alleles

    def __str__(self):
        return f"SNP: {self.rsid}, at chr{self.chromosome}:{self.pos}, with GMAF {self.gmaf}\n"

    def __eq__(self, other):
        if self.chromosome == other.chromosome and self.pos == other.pos and self.rsid == other.rsid:
            return True
        return False

    def clone(self):
        return Snp(self.chromosome, self.pos, self.rsid, self.gmaf, self.alleles)


