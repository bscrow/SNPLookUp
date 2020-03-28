"""Contains custom objects to represent each genomic region of interest.

This module consists of the custom objects to represent each genomic region of
interest, namely the following 3 classes:
    QueryRegion: represents one continuous genomic region of interest from the
    user.
    SubRegion: represents one continuous section of a QueryRegion, after
    excluding excluding regions within each QueryRegions from the user.
    Each QueryRegion contains one or more SubRegion.
    Snp:  represents one Snp locus. Each SubRegion can contain any number of Snps.

"""
from exceptions import *


class QueryRegion(object):
    """Represents one continuous genomic region of interest from the user.

    Attributes:
        chromosome: A string indicating which chromosome this region falls in.
        start_pos: An integer indicating the start position of this region.
        end_pos: An integer indicating the end position of this region.
        label: A string from the user that names this region.
        regions: A list of SubRegions contained in this QueryRegion.
        version:  A string indicating the human genome version used.
    """

    def __init__(self,
                 chromosome,
                 start_pos,
                 end_pos,
                 label,
                 regions,
                 version="hg19"):
        """Inits QueryRegion instance"""
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
        """Defines the string representation of the QueryRegion class."""
        ret_str = f"\nRegion in chromosome {self.chromosome}, between positions {self.start_pos}-{self.end_pos}\n " + f"Region contains the following {len(self.regions)} sub-regions: \n\n "
        for region in self.regions:
            ret_str += str(region) + "\n"
        return ret_str

    def add_snps(self, snps):
        """Adds the input snps.

        Adds the input snps into a SubRegion that contains its locus.

        Args:
            snps: A list of Snp objects to be added.
        """
        for region in self.regions:
            region.add_snps(snps)

    def clone(self):
        """Creates a copy of this QueryRegion class instance.
        Returns:
            A deep copy of this QueryRegion object with the same attributes.
        """
        region_clones = [subreg.clone() for subreg in self.regions]
        return QueryRegion(self.chromosome, self.start_pos, self.end_pos,
                           self.label, region_clones, self.version)

    def get_label(self):
        return self.label

    def get_subregions(self):
        return self.regions

    def get_snps(self):
        """Gets the list of Snps present in this QueryRegion.

        Gets and concatenates the lists of Snps from each of the SubRegions
        contained in this QueryRegion.

        Returns:
            A list of Snp objects contained in this QueryRegion.
        """
        snps = []
        for region in self.regions:
            snps.extend(region.get_snps())
        return snps

    def keep_snps(self, snps):
        """Updates the Snps contained in this QueryRegion.

        Removes all Snps contained in this QueryRegion, before adding the Snps
        in the input list of Snps back.

        Args:
            snps: A list of Snp objects to keep.
        """
        for region in self.regions:
            region.keep_snps(snps)

    def remove_snps(self, snps):
        """Removes all Snps present in the input snp list.

        Args:
            snps:  A list of Snp objects to remove.
        """
        for region in self.regions:
            region.remove_snps(snps)


class SubRegion(object):

    def __init__(self, chromosome, start_pos, end_pos):
        """Inits SubRegion instance"""
        self.chromosome = chromosome
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.snps = []

    def __contains__(self, pos):
        """Defines the usage of the in operator on SubRegion class."""
        return self.start_pos <= pos <= self.end_pos

    def __str__(self):
        """Defines the string representation of the SubRegion class."""
        ret_str = f"Sub-region in chromosome {self.chromosome}, between positions {self.start_pos}-{self.end_pos}\n " + f"Sub-region contains the following {len(self.snps)} SNPs: \n\n"
        for snp in self.snps:
            ret_str += str(snp)
        return ret_str

    def add_snps(self, snps):
        """Adds the input snps.

        Adds the input snps into this SubRegion, if this SubRegion contains the
        Snp loci.

        Args:
            snps: A list of Snp objects to be added.
        """
        for snp in snps:
            if self.end_pos >= snp.pos >= self.start_pos:
                self.snps.append(snp)

    def clear_snps(self):
        self.snps = []

    def clone(self):
        """Creates a copy of this SubRegion class instance.
        Returns:
            A deep copy of this SubRegion object with the same attributes.
        """
        snp_clones = [snp.clone() for snp in self.snps]
        subregion_clone = SubRegion(self.chromosome, self.start_pos,
                                    self.end_pos)
        subregion_clone.add_snps(snp_clones)
        return subregion_clone

    def get_snps(self):
        return self.snps

    def keep_snps(self, snps):
        """Updates the Snps contained in this SubRegion.

        Removes all Snps contained in this SubRegion, before adding the Snps
        in the input list of Snps back.

        Args:
            snps: A list of Snp objects to keep.
        """
        self.clear_snps()
        self.add_snps(snps)

    def remove_snps(self, snps_to_remove):
        """Removes all Snps present in the input snp list.

        Args:
            snps_to_remove:  A list of Snp objects to remove.
        """
        for snp in self.snps:
            if snp in snps_to_remove:
                self.snps.remove(snp)


class Snp(object):

    def __init__(self, chromosome, pos, rsid, gmaf,
                 alleles="/"):  # alleles retrieval not yet implemented
        """Inits Snp instance"""
        self.chromosome = chromosome
        self.pos = pos
        self.rsid = rsid
        self.gmaf = gmaf
        self.alleles = alleles

    def __str__(self):
        """Defines the string representation of the Snp class."""
        return f"SNP: {self.rsid}, at chr{self.chromosome}:{self.pos}, with GMAF {self.gmaf}\n"

    def __eq__(self, other):
        """Defines the conditions for self==other Snp objects."""
        if self.chromosome == other.chromosome and self.pos == other.pos and self.rsid == other.rsid:
            return True
        return False

    def clone(self):
        """Creates a copy of this Snp class instance.
        Returns:
            A copy of this Snp object with the same attributes.
        """
        return Snp(self.chromosome, self.pos, self.rsid, self.gmaf,
                   self.alleles)
