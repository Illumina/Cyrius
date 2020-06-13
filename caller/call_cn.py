#!/usr/bin/env python3
#
# Cyrius: CYP2D6 genotyper
# Copyright (c) 2019-2020 Illumina, Inc.
#
# Author: Xiao Chen <xchen2@illumina.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
from collections import namedtuple
from scipy.stats import poisson, fisher_exact
import pysam

dir_name = os.path.join(os.path.dirname(os.path.dirname(__file__)), "depth_calling")
if os.path.exists(dir_name):
    sys.path.append(dir_name)
from depth_calling.copy_number_call import (
    call_reg1_cn,
    process_raw_call_gc,
    process_raw_call_denovo,
)
from depth_calling.haplotype import get_haplotypes_from_bam, extract_hap


INTRON1_BP_APPROX = 42130500
EXON9_BP_APPROX = 42126611
P_CUTOFF = 0.05

cn_regions = namedtuple(
    "cn_regions", "total_cn exon9_and_downstream exon9_to_intron1 intron1_upstream"
)
CNVTAG_LOOKUP_TABLE = {
    "star5_star5": cn_regions(2, 0, 0, 0),
    "star13_star13": cn_regions(2, 2, 0, 0),
    "star13intron1_star13intron1": cn_regions(2, 2, 2, 0),
    "star5": cn_regions(3, 1, 1, 1),
    "star13": cn_regions(3, 2, 1, 1),
    "star13intron1": cn_regions(3, 2, 2, 1),
    "star5_star5_star68": cn_regions(3, 0, 0, 1),
    "star5_star68": cn_regions(4, 1, 1, 2),
    "cn2": cn_regions(4, 2, 2, 2),
    #'exon9hyb_star5': cn_regions(4, 1, 2, 2),
    "exon9hyb_star5": cn_regions(4, 2, 2, 2),  # process as cn2
    "dup_star13": cn_regions(4, 3, 2, 2),
    "dup_star13intron1": cn_regions(4, 3, 3, 2),
    "star13_star68": cn_regions(4, 2, 1, 2),
    "cn3": cn_regions(5, 3, 3, 3),
    "exon9hyb": cn_regions(5, 2, 3, 3),
    "star68": cn_regions(5, 2, 2, 3),
    "cn4": cn_regions(6, 4, 4, 4),
    "exon9hyb_exon9hyb": cn_regions(6, 2, 4, 4),
    "star68_star68": cn_regions(6, 2, 2, 4),
    "dup_exon9hyb": cn_regions(6, 3, 4, 4),
    "dup_star68": cn_regions(6, 3, 3, 4),
    "exon9hyb_star68": cn_regions(6, 2, 3, 4),
    "cn5": cn_regions(7, 5, 5, 5),
    "exon9hyb_exon9hyb_exon9hyb": cn_regions(7, 2, 5, 5),
    "star68_star68_star68": cn_regions(7, 2, 2, 5),
    "cn6": cn_regions(8, 6, 6, 6),
    "exon9hyb_exon9hyb_exon9hyb_exon9hyb": cn_regions(8, 2, 6, 6),
    "star68_star68_star68_star68": cn_regions(8, 2, 2, 6),
}

# For these variants the region is clean and we use a less stringent minimum read cutoff
CLEAN_VAR = [
    "g.42129809T>C",
    "g.42129819G>T",
    "g.42128945C>T",
    "g.42126611C>G",
    "g.42130692G>A",
    "g.42127941G>A",
]
# These are noisy (mostly gene conversion) variants that may have misalignments
# resulting in strand bias
NOISY_VAR = [
    "g.42127473C>T",
    "g.42128181A>T",
    "g.42128185C>T",
    "g.42129042T>C",
    "g.42129174C>A",
    "g.42129180A>T",
]


def get_total_cn_per_site(cnvtag, var_db, var_list):
    """
    For variants in non-homology regions, get the total expected CN
    at each site based on the CNV configuration.
    """
    # total number of variant sites in non-homology regions
    num_var_sites = len(var_db.dsnp1)
    variant_names = var_list[:num_var_sites]
    num_var_sites_before_intron1bp = 0
    num_var_sites_after_exon9bp = 0
    for var_name in variant_names:
        var_pos = int(var_name[2:10])
        if var_pos >= INTRON1_BP_APPROX:
            num_var_sites_before_intron1bp += 1
        if var_pos <= EXON9_BP_APPROX:
            num_var_sites_after_exon9bp += 1

    if cnvtag not in CNVTAG_LOOKUP_TABLE:
        return None
    cn_pattern = CNVTAG_LOOKUP_TABLE[cnvtag]
    cn_list = []
    for _ in range(num_var_sites_after_exon9bp):
        cn_list.append(cn_pattern.exon9_and_downstream)
    for _ in range(
        num_var_sites - num_var_sites_before_intron1bp - num_var_sites_after_exon9bp
    ):
        cn_list.append(cn_pattern.exon9_to_intron1)
    for _ in range(num_var_sites_before_intron1bp):
        cn_list.append(cn_pattern.intron1_upstream)
    return cn_list


def call_cn_snp(total_cn, lsnp1, lsnp2, threshold=0.6):
    """
    Call CN for SNP sites between CYP2D6 and CYP2D7.
    Use a loose cutoff as this is for CNV/hybrid group calling.
    """
    cn_prob = []
    for i, count1 in enumerate(lsnp1):
        count2 = lsnp2[i]
        cn_prob.append(call_reg1_cn(total_cn, count1, count2))
    cn_call = process_raw_call_gc(cn_prob, threshold)
    return cn_call


def call_cn_var_homo(total_cn, lsnp1, lsnp2):
    """
    Call CN for variant sites in homology regions.
    """
    cn_prob = []
    for i, count1 in enumerate(lsnp1):
        count2 = lsnp2[i]
        cn_prob.append(call_reg1_cn(total_cn, count1, count2, 4))
    cn_call = []
    for site_call in process_raw_call_denovo(cn_prob, 0.8, 0.65):
        if site_call is None:
            cn_call.append(None)
        else:
            cn_call.append(min(site_call, total_cn - 2))
    return cn_call


def call_cn_var(cnvtag, var_alt, var_ref, alt_forward, alt_reverse, var_list, var_db):
    """
    Call CN for variant sites in non-homology regions.
    Use different minimum read cutoffs for clean variant sites and other sites.
    Total CN at each site is also considered during filtering.
    """
    total_cn = get_total_cn_per_site(cnvtag, var_db, var_list)
    assert total_cn is not None
    cn_prob = []

    for i, forward in enumerate(alt_forward):
        reverse = alt_reverse[i]
        total_ref = var_ref[i]
        total_var = var_alt[i]
        if total_var > 0 and var_list[i] in NOISY_VAR:
            ntotal = forward + reverse
            oddsratio, pvalue = fisher_exact(
                [[forward, reverse], [ntotal / 2, ntotal / 2]]
            )
            if pvalue < P_CUTOFF or forward <= 1 or reverse <= 1:
                total_var = 0

        if var_list[i] in CLEAN_VAR:
            cn_prob.append(call_reg1_cn(total_cn[i], total_var, total_ref, 2))
        elif var_list[i] in NOISY_VAR:
            cn_prob.append(call_reg1_cn(total_cn[i], total_var, total_ref, 7))
        else:
            cn_prob.append(call_reg1_cn(total_cn[i], total_var, total_ref, 4))
    cn_call = process_raw_call_denovo(cn_prob, 0.8, 0.65, total_cn)
    return cn_call


def good_read(read):
    """
    Define read filters
    """
    return read.is_secondary == 0 and read.is_supplementary == 0


def get_allele_counts_42128936(bamfile_handle, genome):
    """
    Search for the inserstions at 42128936 defining
    *30/*40/*58 in read sequences
    """
    long_ins_read = 0
    short_ins_read = 0
    ref_read = 0
    dregion = {
        "19": ("chr22", 42524850, 42524980),
        "37": ("22", 42524850, 42524980),
        "38": ("chr22", 42128848, 42128978),
    }
    region = dregion[genome]
    for read in bamfile_handle.fetch(region[0], region[1], region[2]):
        seq = read.query_sequence
        if good_read(read):
            if "TGGGGCGAAAGGGGCGAAAGGGGCGAAAGGGGCGT" in seq:
                long_ins_read += 1
            elif "TTGGGGCGAAAGGGGCGAAAGGGGCGTC" in seq:
                short_ins_read += 1
            elif "TTGGGGCGAAAGGGGCGTC" in seq:
                ref_read += 1
    return (ref_read, long_ins_read, short_ins_read)


def call_exon9gc(d6_count, d7_count, full_length_cn):
    """
    Call exon 9 conversion
    """
    lsnp1 = [d6_count]
    lsnp2 = [d7_count]

    if full_length_cn is not None:
        full_length_cn = int(full_length_cn)
        cn_prob = []
        for i, count1 in enumerate(lsnp1):
            count2 = lsnp2[i]
            cn_prob.append(call_reg1_cn(full_length_cn, count1, count2, 3))
        cn_prob_processed_stringent = process_raw_call_gc(cn_prob, 0.9)

    return cn_prob_processed_stringent[0]


def call_var42126938(bamfile, cnvtag, site42126938, base_db, target_positions):
    """
    Call variant g.42126938C>T (gene conversion variant in homology region) 
    based on read depth and phased haplotypes
    """
    dcn = {"star5": 3, "cn2": 4}
    assert cnvtag in dcn
    full_length_cn = dcn[cnvtag]
    d6_cn = call_cn_snp(full_length_cn, [site42126938[0]], [site42126938[1]], 0.8)[0]
    var_called = []
    # Whether g.42126938C>T is on the same haplotype as g.42126611C>G
    G_haplotype = False
    if d6_cn is not None and d6_cn < full_length_cn - 2:
        haplotype_per_read = get_haplotypes_from_bam(bamfile, base_db, target_positions)
        recombinant_read_count = extract_hap(haplotype_per_read, [0, 2])
        if "12" in recombinant_read_count and sum(recombinant_read_count["12"]) > 1:
            G_hap_count = extract_hap(haplotype_per_read, [1, 2])
            for _ in range(full_length_cn - 2 - d6_cn):
                var_called.append("g.42126938C>T")
            if "12" in G_hap_count and sum(G_hap_count["12"]) > 1:
                G_haplotype = True
    return var_called, G_haplotype


def get_called_variants(var_list, cn_prob_processed, starting_index=0):
    """
    Return called variants based on called copy number and list of variant names
    """
    total_callset = []
    if starting_index != 0:
        assert len(var_list) == len(cn_prob_processed) + starting_index
    for i, cn_called in enumerate(cn_prob_processed):
        if cn_called is not None and cn_called != 0:
            for _ in range(cn_called):
                total_callset.append(var_list[i + starting_index])
    return total_callset
