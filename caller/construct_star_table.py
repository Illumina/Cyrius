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

from collections import namedtuple


# Exon 9 gene conversion
EXON9GC_ALLELES = ["*36", "*4N", "*57", "*83"]
EXON9GC_PAIR_ALLELES = {"*36": "*10", "*4N": "*4A"}


def make_hap_dic(variant_list, star_set, hap_dic):
    """
    Update variant-to-diplotype dictionary
    """
    while "NA" in variant_list:
        variant_list.remove("NA")
    var_list_joined = "_".join(variant_list)
    if variant_list == []:
        var_list_joined = "NA"
    hap_dic.setdefault(var_list_joined, [])
    if star_set not in hap_dic[var_list_joined]:
        hap_dic[var_list_joined].append(star_set)


def get_hap_table(hap_table):
    """
    Return all the possible variant tables based on the star allele definition file.
    """
    # CN=1, one copy
    dhap = {}
    with open(hap_table) as f:
        for line in f:
            at = line.strip().split()
            star_id = at[0]
            variant_list = sorted(at[1:-2])
            var_list_joined = "_".join(variant_list)
            dhap.setdefault(var_list_joined, star_id)

    # CN=2, two copies
    dhap2 = {}
    for star1 in dhap:
        for star2 in dhap:
            variant_list = sorted(star1.split("_") + star2.split("_"))
            star_set = "_".join(sorted([dhap[star1], dhap[star2]]))
            make_hap_dic(variant_list, star_set, dhap2)

    # CN=3, three copies, for exon9hyb cases
    dhap3 = {}
    for star1 in dhap:
        for star2 in dhap:
            for star3 in dhap:
                variant_list = sorted(
                    star1.split("_") + star2.split("_") + star3.split("_")
                )
                star_set = "_".join(sorted([dhap[star1], dhap[star2], dhap[star3]]))
                make_hap_dic(variant_list, star_set, dhap3)

    # CN=3, three copies, no hybrid, only duplication
    # Assume duplicated copies are identical
    dhap3pair = {}
    for star1 in dhap:
        for star2 in dhap:
            variant_list = sorted(
                star1.split("_") + star2.split("_") + star2.split("_")
            )
            star_set = "_".join(sorted([dhap[star1], dhap[star2], dhap[star2]]))
            make_hap_dic(variant_list, star_set, dhap3pair)

    # CN=4, four copies, no hybrid, only duplication
    # Assume duplicated copies are identical
    # can be star1x2 + star2x2 or star1 + star2x3
    dhap4pair = {}
    for star1 in dhap:
        for star2 in dhap:
            variant_list = sorted(
                star1.split("_")
                + star1.split("_")
                + star2.split("_")
                + star2.split("_")
            )
            star_set = "_".join(
                sorted([dhap[star1], dhap[star1], dhap[star2], dhap[star2]])
            )
            make_hap_dic(variant_list, star_set, dhap4pair)

            variant_list = sorted(
                star1.split("_")
                + star2.split("_")
                + star2.split("_")
                + star2.split("_")
            )
            star_set = "_".join(
                sorted([dhap[star1], dhap[star2], dhap[star2], dhap[star2]])
            )
            make_hap_dic(variant_list, star_set, dhap4pair)

    # exon9hybx2. limit the search to EXON9GC_ALLELES
    dhap_exon9_x2 = {}
    for star1 in dhap:
        for star2 in dhap:
            if dhap[star1] in EXON9GC_ALLELES and dhap[star2] in EXON9GC_ALLELES:
                for star3 in dhap:
                    for star4 in dhap:
                        variant_list = sorted(
                            star1.split("_")
                            + star2.split("_")
                            + star3.split("_")
                            + star4.split("_")
                        )
                        star_set = "_".join(
                            sorted([dhap[star1], dhap[star2], dhap[star3], dhap[star4]])
                        )
                        make_hap_dic(variant_list, star_set, dhap_exon9_x2)

    # exon9hybx3. limit the search to EXON9GC_ALLELES
    dhap_exon9_x3 = {}
    for star1 in dhap:
        for star2 in dhap:
            for star3 in dhap:
                if (
                    dhap[star1] in EXON9GC_ALLELES
                    and dhap[star2] in EXON9GC_ALLELES
                    and dhap[star3] in EXON9GC_ALLELES
                ):
                    for star4 in dhap:
                        for star5 in dhap:
                            variant_list = sorted(
                                star1.split("_")
                                + star2.split("_")
                                + star3.split("_")
                                + star4.split("_")
                                + star5.split("_")
                            )
                            star_set = "_".join(
                                sorted(
                                    [
                                        dhap[star1],
                                        dhap[star2],
                                        dhap[star3],
                                        dhap[star4],
                                        dhap[star5],
                                    ]
                                )
                            )
                            make_hap_dic(variant_list, star_set, dhap_exon9_x3)

    # dup_exon9hyb. For the exon9hyb part, limit the search to EXON9GC_PAIR_ALLELES
    dhap_dup_exon9 = {}
    for star1 in dhap:
        for star2 in dhap:
            if (
                dhap[star1] in EXON9GC_PAIR_ALLELES
                and dhap[star2] == EXON9GC_PAIR_ALLELES[dhap[star1]]
            ):
                for star3 in dhap:
                    for star4 in dhap:
                        variant_list = sorted(
                            star1.split("_")
                            + star2.split("_")
                            + star3.split("_")
                            + star4.split("_")
                        )
                        star_set = "_".join(
                            sorted([dhap[star1], dhap[star2], dhap[star3], dhap[star4]])
                        )
                        make_hap_dic(variant_list, star_set, dhap_dup_exon9)

    star_combinations = namedtuple(
        "star_combinations",
        "dhap dhap2 dhap3 dhap3pair dhap4pair dhap_exon9_x2 dhap_exon9_x3 dhap_dup_exon9",
    )
    return star_combinations(
        dhap,
        dhap2,
        dhap3,
        dhap3pair,
        dhap4pair,
        dhap_exon9_x2,
        dhap_exon9_x3,
        dhap_dup_exon9,
    )
