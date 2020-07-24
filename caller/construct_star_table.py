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
EXON9GC_ALLELES = ["*36", "*4.013", "*57", "*83"]
EXON9GC_PAIR_ALLELES = {"*36": "*10", "*4.013": "*4"}


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
    dstar = {}
    with open(hap_table) as f:
        for line in f:
            at = line.strip().split()
            star_id = at[0]
            variant_list = sorted(at[1:-1])
            var_list_joined = "_".join(variant_list)
            dhap.setdefault(var_list_joined, star_id)
            dstar.setdefault(star_id, var_list_joined)

    # CN=2, two copies
    dhap2 = {}
    for star1 in dstar:
        for star2 in dstar:
            variant_list = sorted(dstar[star1].split("_") + dstar[star2].split("_"))
            star_set = "_".join(sorted([star1, star2]))
            make_hap_dic(variant_list, star_set, dhap2)

    # CN=3, three copies, for exon9hyb cases
    dhap3 = {}
    for star1 in dstar:
        for star2 in dstar:
            for star3 in dstar:
                variant_list = sorted(
                    dstar[star1].split("_")
                    + dstar[star2].split("_")
                    + dstar[star3].split("_")
                )
                star_set = "_".join(sorted([star1, star2, star3]))
                make_hap_dic(variant_list, star_set, dhap3)

    # CN=3, three copies, no hybrid, only duplication
    # Assume duplicated copies are identical
    dhap3pair = {}
    for star1 in dstar:
        for star2 in dstar:
            variant_list = sorted(
                dstar[star1].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
            )
            star_set = "_".join(sorted([star1, star2, star2]))
            make_hap_dic(variant_list, star_set, dhap3pair)

    # CN=4, four copies, no hybrid, only duplication
    # Assume duplicated copies are identical
    # can be star1x2 + star2x2 or star1 + star2x3
    dhap4pair = {}
    for star1 in dstar:
        for star2 in dstar:
            variant_list = sorted(
                dstar[star1].split("_")
                + dstar[star1].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
            )
            star_set = "_".join(sorted([star1, star1, star2, star2]))
            make_hap_dic(variant_list, star_set, dhap4pair)

            variant_list = sorted(
                dstar[star1].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
            )
            star_set = "_".join(sorted([star1, star2, star2, star2]))
            make_hap_dic(variant_list, star_set, dhap4pair)

    # CN=5, five copies, no hybrid, only duplication
    # Assume duplicated copies are identical
    # can be star1x2 + star2x3 or star1 + star2x4
    dhap5pair = {}
    for star1 in dstar:
        for star2 in dstar:
            variant_list = sorted(
                dstar[star1].split("_")
                + dstar[star1].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
            )
            star_set = "_".join(sorted([star1, star1, star2, star2, star2]))
            make_hap_dic(variant_list, star_set, dhap5pair)

            variant_list = sorted(
                dstar[star1].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
            )
            star_set = "_".join(sorted([star1, star2, star2, star2, star2]))
            make_hap_dic(variant_list, star_set, dhap5pair)

    # CN=6, six copies, no hybrid, only duplication
    # Assume duplicated copies are identical
    # can be star1x2 + star2x4 or star1 + star2x5 or star1x3 + star2x3
    dhap6pair = {}
    for star1 in dstar:
        for star2 in dstar:
            variant_list = sorted(
                dstar[star1].split("_")
                + dstar[star1].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
            )
            star_set = "_".join(sorted([star1, star1, star2, star2, star2, star2]))
            make_hap_dic(variant_list, star_set, dhap6pair)

            variant_list = sorted(
                dstar[star1].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
            )
            star_set = "_".join(sorted([star1, star2, star2, star2, star2, star2]))
            make_hap_dic(variant_list, star_set, dhap6pair)

            variant_list = sorted(
                dstar[star1].split("_")
                + dstar[star1].split("_")
                + dstar[star1].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
                + dstar[star2].split("_")
            )
            star_set = "_".join(sorted([star1, star1, star1, star2, star2, star2]))
            make_hap_dic(variant_list, star_set, dhap6pair)

    # exon9hybx2. limit the search to EXON9GC_ALLELES
    dhap_exon9_x2 = {}
    for star1 in EXON9GC_ALLELES:
        for star2 in EXON9GC_ALLELES:
            if star1 in dstar and star2 in dstar:
                for star3 in dstar:
                    for star4 in dstar:
                        variant_list = sorted(
                            dstar[star1].split("_")
                            + dstar[star2].split("_")
                            + dstar[star3].split("_")
                            + dstar[star4].split("_")
                        )
                        star_set = "_".join(sorted([star1, star2, star3, star4]))
                        make_hap_dic(variant_list, star_set, dhap_exon9_x2)

    # exon9hybx3. limit the search to EXON9GC_ALLELES
    dhap_exon9_x3 = {}
    for star1 in EXON9GC_ALLELES:
        for star2 in EXON9GC_ALLELES:
            for star3 in EXON9GC_ALLELES:
                if star1 in dstar and star2 in dstar and star3 in dstar:
                    for star4 in dstar:
                        for star5 in dstar:
                            variant_list = sorted(
                                dstar[star1].split("_")
                                + dstar[star2].split("_")
                                + dstar[star3].split("_")
                                + dstar[star4].split("_")
                                + dstar[star5].split("_")
                            )
                            star_set = "_".join(
                                sorted([star1, star2, star3, star4, star5])
                            )
                            make_hap_dic(variant_list, star_set, dhap_exon9_x3)

    # exon9hybx4. limit the search to *10 and *36
    dhap_exon9_x4 = {}
    variant_list = sorted(
        dstar["*10"].split("_")
        + dstar["*10"].split("_")
        + dstar["*36"].split("_")
        + dstar["*36"].split("_")
        + dstar["*36"].split("_")
        + dstar["*36"].split("_")
    )
    star_set = "_".join(sorted(["*10", "*10", "*36", "*36", "*36", "*36"]))
    make_hap_dic(variant_list, star_set, dhap_exon9_x4)

    # dup_exon9hyb. For the exon9hyb part, limit the search to EXON9GC_PAIR_ALLELES
    dhap_dup_exon9 = {}
    for star1 in EXON9GC_PAIR_ALLELES:
        star2 = EXON9GC_PAIR_ALLELES[star1]
        for star3 in dstar:
            for star4 in dstar:
                variant_list = sorted(
                    dstar[star1].split("_")
                    + dstar[star2].split("_")
                    + dstar[star3].split("_")
                    + dstar[star4].split("_")
                )
                star_set = "_".join(sorted([star1, star2, star3, star4]))
                make_hap_dic(variant_list, star_set, dhap_dup_exon9)

    star_combinations = namedtuple(
        "star_combinations",
        "dhap dhap2 dhap3 dhap3pair dhap4pair dhap5pair dhap6pair dhap_exon9_x2 dhap_exon9_x3 dhap_exon9_x4 dhap_dup_exon9",
    )
    return star_combinations(
        dhap,
        dhap2,
        dhap3,
        dhap3pair,
        dhap4pair,
        dhap5pair,
        dhap6pair,
        dhap_exon9_x2,
        dhap_exon9_x3,
        dhap_exon9_x4,
        dhap_dup_exon9,
    )
