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
from collections import namedtuple
import re

CNVTAG_TO_GENOTYPE = {
    "star5_star5": "*5/*5",
    "star13_star13": "*13/*13",
    "star13intron1_star13intron1": "*13/*13",
    "star5_star5_star68": "*5/*68",
}
# These suballeles below are not converted to main alleles as
# they reflect SVs.
KEPT_SUBALLELES = ["*4.013"]
# Rare alleles lead to nonunique diplotypes. Select against these
# when there are nonunique calls.
RARE_ALLELES = ["*34", "*39", "*4.009", "*139"]
raw_star = namedtuple("raw_star", "call_info candidate star_call")


def get_var_list(var_observed):
    """
    Merge called variants
    """
    if var_observed == []:
        return "NA"
    return "_".join(var_observed)


def convert_to_main_allele(list_of_star):
    """
    Convert suballeles to main alleles
    """
    converted_list = set()
    for stars in list_of_star:
        star_split = stars.split("_")
        converted_star = []
        for star in star_split:
            if star not in KEPT_SUBALLELES:
                converted_star.append(star.split(".")[0])
            else:
                converted_star.append(star)
        converted_list.add("_".join(sorted(converted_star)))
    return list(converted_list)


def get_star(var_observed, dic):
    """
    Return the star allele call and assign tags like unique and nonunique.
    Output: information about the match, raw star allele calls,
    parsed (selected) star allele calls
    """
    var_observed = sorted(var_observed)
    var_list = get_var_list(var_observed)
    match_tag = None
    raw_stars = []
    processed_stars = []

    # Compare called variants against the star combination table
    if var_list not in dic:
        # No match
        match_tag = "no_match"
    else:
        # One copy. Only one haplotype to match
        if "*" in dic[var_list]:
            match_tag = "unique_match"
            raw_stars = [dic[var_list]]
            processed_stars = convert_to_main_allele(raw_stars)
        # More than one match
        elif len(dic[var_list]) > 1:
            raw_stars = dic[var_list]
            processed_stars = convert_to_main_allele(raw_stars)
            if len(processed_stars) == 1:
                match_tag = "unique_star"
            else:
                rare_stars_found = []
                for haplotype in raw_stars:
                    for rare_allele in RARE_ALLELES:
                        if rare_allele in haplotype:
                            rare_stars_found.append(haplotype)
                processed_stars = convert_to_main_allele(
                    [a for a in raw_stars if a not in rare_stars_found]
                )
                if len(processed_stars) == 1:
                    match_tag = "pick_common_allele"
                else:
                    match_tag = "more_than_one_match"
        else:
            # Unique match
            match_tag = "unique_match"
            raw_stars = dic[var_list]
            processed_stars = convert_to_main_allele([raw_stars[0]])

    return raw_star(match_tag, raw_stars, processed_stars)


def call_star68(var_observed, cnvcall, dic):
    """
    Return the star allele call for star68-related cases
    Output: information about the match, raw star allele calls, parsed (selected) star allele calls
    """
    matchtag = []
    cn = cnvcall.split("_").count("star68")
    # g.42130692G>A is almost always in the D6 part of the hybrid gene
    # try removing it/them and then match against the star table
    count0692 = var_observed.count("g.42130692G>A")
    matchtag.append(get_star(var_observed, dic))
    for _ in range(min(cn, count0692)):
        var_observed.remove("g.42130692G>A")
        matchtag.append(get_star(var_observed, dic))
    # Unique match
    if len(matchtag) == 1:
        return matchtag[0]
    num_match = 0
    for tag in matchtag:
        if tag[0] != "no_match":
            num_match += 1
    # no match
    if num_match == 0:
        return raw_star("no_match", [], [])
    hap_list = []
    for tag in matchtag:
        hap_list += tag[-1]
    if len(set(hap_list)) == 1:
        return raw_star("unique_match", hap_list, [hap_list[0]])

    rare_stars_found = []
    for haplotype in hap_list:
        for rare_allele in RARE_ALLELES:
            if rare_allele in haplotype:
                rare_stars_found.append(haplotype)
    processed_stars = [a for a in hap_list if a not in rare_stars_found]
    if len(set(processed_stars)) == 1:
        return raw_star("pick_common_allele", hap_list, [processed_stars[0]])

    # The one with g.42130692G>A removed is the most likely case
    return raw_star("pick_common_allele", hap_list, [hap_list[-1]])


def get_dic(cnvcall, star_combinations):
    """
    Return the variant table to use for each CNV/hybrid group
    """
    if cnvcall in ["star5", "star13", "star5_star68", "star13_star68"]:
        return star_combinations.dhap
    elif cnvcall in [
        "cn2",
        "star13intron1",
        "dup_star13",
        "exon9hyb_star5",
    ] or cnvcall.startswith("star68"):
        return star_combinations.dhap2
    elif cnvcall in ["exon9hyb", "exon9hyb_star68"]:
        return star_combinations.dhap3
    elif cnvcall in ["cn3", "dup_star13intron1", "dup_star68"]:
        return star_combinations.dhap3pair
    elif cnvcall == "cn4":
        return star_combinations.dhap4pair
    elif cnvcall == "cn5":
        return star_combinations.dhap5pair
    elif cnvcall == "cn6":
        return star_combinations.dhap6pair
    elif cnvcall == "exon9hyb_exon9hyb":
        return star_combinations.dhap_exon9_x2
    elif cnvcall == "exon9hyb_exon9hyb_exon9hyb":
        return star_combinations.dhap_exon9_x3
    elif cnvcall == "exon9hyb_exon9hyb_exon9hyb_exon9hyb":
        return star_combinations.dhap_exon9_x4
    elif cnvcall == "dup_exon9hyb":
        return star_combinations.dhap_dup_exon9
    return None


def get_final_call_clean(final_call, cnvcall, spacer_cn):
    """
    Clean up final call to report diplotypes in *#/*# format whenever possible.
    """
    # zero or more than one set
    final_call = sorted(final_call)
    if len(final_call) == 2 and cnvcall == "cn2":
        diplotype1 = final_call[0].split("_")
        diplotype2 = final_call[1].split("_")
        return "/".join(diplotype1) + ";" + "/".join(diplotype2)
    if final_call == [] or len(final_call) > 1:
        if final_call == ["*10_*10_*4.013", "*10_*36_*4"]:
            return "*4/*36+*10"
        return ";".join(final_call)

    called_stars = final_call[0]
    if cnvcall == "star5_star68":
        if called_stars == "*4":
            return "*5/*68+*4"
        return "*68/" + called_stars

    if cnvcall == "star13_star68":
        if called_stars == "*4":
            return "*13/*68+*4"
        return "*13_" + called_stars + "_*68"

    split_call = called_stars.split("_")
    if cnvcall == "cn2":
        return "/".join(split_call)

    if cnvcall == "star13intron1":
        if split_call[0] == "*2":
            return "*13/" + split_call[1]
        if split_call[1] == "*2":
            return "*13/" + split_call[0]
        return None
    if cnvcall == "dup_star13intron1":
        unique_stars = list(set(split_call))
        if len(unique_stars) == 2:
            if split_call.count(unique_stars[0]) == 2:
                return "*13+" + unique_stars[0] + "/" + unique_stars[1]
            elif split_call.count(unique_stars[1]) == 2:
                return "*13+" + unique_stars[1] + "/" + unique_stars[0]
        elif len(unique_stars) == 1:
            return "*13+" + unique_stars[0] + "/" + unique_stars[0]
        return None
    if cnvcall == "dup_star13":
        if spacer_cn is not None and spacer_cn == 1:
            split_call.append("*13")
            return "_".join(split_call)
        # default to cn2
        return "/".join(split_call)

    if cnvcall in ["star5", "star13"]:
        return called_stars + "/" + "*" + cnvcall[4:]

    if cnvcall in ["cn3", "cn4", "cn5", "cn6"]:
        dup_allele = []
        call_set = []
        for star_allele in split_call:
            if split_call.count(star_allele) > 1:
                if star_allele not in dup_allele:
                    call_set.append(
                        star_allele + "x" + str(split_call.count(star_allele))
                    )
                    dup_allele.append(star_allele)
            else:
                call_set.append(star_allele)
        if len(call_set) == 2:
            return "/".join(call_set)
        if len(call_set) == 1:
            var = split_call[0]
            if cnvcall == "cn3":
                return var + "/" + var + "x2"
            if cnvcall == "cn4":
                return var + "x2/" + var + "x2"
            if cnvcall == "cn5":
                return var + "x2/" + var + "x3"
            if cnvcall == "cn6":
                return var + "x3/" + var + "x3"
        return "_".join(call_set)

    if cnvcall == "exon9hyb_star5":
        # for these two cases, check spacer CN to determine if they are on the same chromosome
        if split_call == ["*10", "*10"]:
            if spacer_cn is not None and spacer_cn > 2:
                return "*5/*36+*10"
            else:
                return None
        if split_call == ["*4", "*4"]:
            if spacer_cn is not None and spacer_cn > 2:
                return "*5/*4.013+*4"
        return "/".join(split_call)

    if cnvcall == "exon9hyb":
        if "*4" in split_call and "*4.013" in split_call:
            remain_index = [
                n
                for n in range(3)
                if n not in [split_call.index("*4"), split_call.index("*4.013")]
            ]
            assert len(remain_index) == 1
            return split_call[remain_index[0]] + "/*4.013+*4"
        if "*10" in split_call and "*36" in split_call:
            remain_index = [
                n
                for n in range(3)
                if n not in [split_call.index("*10"), split_call.index("*36")]
            ]
            assert len(remain_index) == 1
            return split_call[remain_index[0]] + "/*36+*10"
        if split_call.count("*36") == 2:
            remain_star = [a for a in split_call if a != "*36"]
            return "*36+*36/" + remain_star[0]

    if cnvcall == "dup_exon9hyb":
        var = []
        if "*36" in split_call:
            var = [a for a in split_call if a not in ["*10", "*36"]]
        elif "*4.013" in split_call:
            var = [a for a in split_call if a not in ["*4", "*4.013"]]
        if len(var) == 2:
            split_call.remove(var[0])
            split_call.remove(var[1])
            if len(set(var)) == 1:
                return var[0] + "x2/" + "+".join(sorted(split_call, reverse=True))
            else:
                return "+".join(var) + "/" + "+".join(sorted(split_call, reverse=True))
        if var == []:
            if called_stars == "*10_*10_*10_*36":
                return "*10x2/*36+*10"
            if called_stars == "*4_*4_*4_*4.013":
                return "*4x2/*4.013+*4"

    if (
        cnvcall == "exon9hyb_exon9hyb"
        or cnvcall == "exon9hyb_exon9hyb_exon9hyb"
        or cnvcall == "exon9hyb_exon9hyb_exon9hyb_exon9hyb"
    ):
        if called_stars == "*4_*4_*4.013_*4.013":
            return "*4.013+*4/*4.013+*4"
        if called_stars == "*10_*10_*36_*36":
            return "*36+*10/*36+*10"
        if called_stars == "*10_*36_*36_*36":
            return "*36+*10/*36+*36"
        if called_stars == "*10_*10_*36_*36_*36":
            return "*36+*10/*36+*36+*10"
        if called_stars == "*10_*10_*36_*36_*36_*36":
            return "*36+*36+*10/*36+*36+*10"
        if (
            cnvcall == "exon9hyb_exon9hyb_exon9hyb"
            and "*10" in split_call
            and "*83" in split_call
            and split_call.count("*36") == 2
        ):
            split_call.remove("*10")
            split_call.remove("*36")
            split_call.remove("*36")
            split_call.remove("*83")
            return split_call[0] + "/*36+*36+*83+*10"
        var = [a for a in split_call if a not in ["*10", "*36", "*83"]]
        if len(var) == 1:
            split_call.remove(var[0])
            return var[0] + "/" + "+".join(sorted(split_call, reverse=True))

    if "star68" in cnvcall:
        cn = cnvcall.split("_").count("star68")
        if len(set(cnvcall.split("_"))) == 1:
            if "*4" in split_call:
                var = [a for a in split_call if a != "*4"]
                if len(var) == 1:
                    genotype = var[0] + "/"
                    for _ in range(cn):
                        genotype += "*68+"
                    genotype += "*4"
                    return genotype
                elif split_call == ["*4", "*4"]:
                    if cn == 1:
                        return "*4/*68+*4"
                    if cn == 2:
                        return "*68+*4/*68+*4"
                    elif cn == 3:
                        return "*68+*4/*68+*68+*4"
                    elif cn == 4:
                        return "*68+*68+*4/*68+*68+*4"
            # *45 is found in tandem with *68 in Africans
            if "*45" in split_call:
                var = [a for a in split_call if a != "*45"]
                if len(var) == 1:
                    genotype = var[0] + "/"
                    for _ in range(cn):
                        genotype += "*68+"
                    genotype += "*45"
                    return genotype

        if cnvcall == "dup_star68":
            var = [a for a in split_call if a not in ["*4", "*68"]]
            if len(var) == 2 and len(set(var)) == 1:
                return var[0] + "x2/*68+*4"
            if var == [] and called_stars == "*4_*4_*4":
                return "*4x2/*68+*4"
        if cnvcall == "exon9hyb_star68":
            if called_stars == "*4_*4_*4.013":
                return "*4.013+*4/*68+*4"
        for _ in range(cn):
            split_call.append("*68")
        return "_".join(split_call)

    return called_stars


def update_variants(var_observed, cnvcall, exon9):
    """
    Update variants based on called CNV.
    """
    # g.42129809T>C and g.42129819G>T should have the same CN.
    if "g.42129809T>C" in var_observed:
        for _ in range(
            var_observed.count("g.42129809T>C") - var_observed.count("g.42129819G>T")
        ):
            var_observed.append("g.42129819G>T")

    # g.42127556T>C is included in g.42127565T>C definition for *108.
    if "g.42127565T>C" in var_observed:
        for _ in range(
            var_observed.count("g.42127565T>C") - var_observed.count("g.42127556T>C")
        ):
            var_observed.append("g.42127556T>C")

    # g.42126611C>G is in the D6 part of the hybrid gene.
    if "star13" in cnvcall and "intron1" not in cnvcall:
        if "g.42126611C>G" in var_observed:
            var_observed.remove("g.42126611C>G")

    if "exon9hyb" in cnvcall and cnvcall != "exon9hyb_star5":
        cn = cnvcall.split("_").count("exon9hyb")
        if "exon9gc" not in var_observed:
            for _ in range(cn):
                # exon9gc always comes with g.42126611C>G
                var_observed.append("exon9gc")
                var_observed.append("g.42126611C>G")
        # Add these variants if they are not called to the sufficient copy number.
        # These variants belong to *10 or *4
        for var_to_add in ["g.42130692G>A"]:
            if var_to_add in var_observed and var_observed.count(var_to_add) <= cn:
                var_observed.append(var_to_add)
        for var_to_add in ["g.42129754G>A"]:
            if (
                var_to_add in var_observed
                and var_observed.count(var_to_add) <= cn
                and "g.42128945C>T" not in var_observed
            ):
                var_observed.append(var_to_add)
        for var_to_add in ["g.42128945C>T", "g.42129809T>C", "g.42129819G>T"]:
            if (
                var_to_add in var_observed
                and var_observed.count(var_to_add) <= cn
                and "g.42129754G>A" not in var_observed
            ):
                var_observed.append(var_to_add)

    exon9_values = namedtuple(
        "exon9_values", "exon9_cn exon9cn_in_consensus exon9_raw_site1 exon9_raw_site2"
    )

    # exon 9 gene conversion by itself, without the fusion
    if "exon9hyb" in cnvcall or cnvcall == "cn2":
        if (
            exon9.exon9cn_in_consensus is not None
            and exon9.exon9_cn is not None
            and exon9.exon9cn_in_consensus > exon9.exon9_cn
            and exon9.exon9_cn <= 1
        ):
            if "g.42126611C>G" not in var_observed or (
                min(exon9.exon9_raw_site1, exon9.exon9_raw_site2) < 1.15
                and max(exon9.exon9_raw_site1, exon9.exon9_raw_site2) < 1.2
            ):
                for _ in range(exon9.exon9cn_in_consensus - exon9.exon9_cn):
                    var_observed.append("exon9gc")
                if "g.42126611C>G" not in var_observed:
                    for _ in range(exon9.exon9cn_in_consensus - exon9.exon9_cn):
                        var_observed.append("g.42126611C>G")

    while "NA" in var_observed:
        var_observed.remove("NA")

    return var_observed


def match_star(
    var_observed,
    cnvcall,
    spacer_cn,
    star_combinations,
    exon9,
    var42126938_G_haplotype,
    var42127803_diff_haplotype,
):
    """
    Return the star allele call based on the called cnv/hybrid group and small variants
    """
    star_call = namedtuple("star_call", "call_info variants_called raw_call clean_call")

    if cnvcall in CNVTAG_TO_GENOTYPE:
        called_genotype = CNVTAG_TO_GENOTYPE[cnvcall]
        return star_call("unique_match", "", called_genotype, called_genotype)

    dic = get_dic(cnvcall, star_combinations)
    if dic is None:
        return star_call(None, None, None, None)
    # print(dic)

    var_observed = update_variants(var_observed, cnvcall, exon9)

    if "star68" not in cnvcall:
        matchtag = get_star(var_observed, dic)
        if cnvcall in ["cn3", "cn4", "cn5", "cn6"] and "no_match" in matchtag.call_info:
            # for cn3 and cn4, try adding a variant and match again
            variant_tried = []
            matched_calls = []
            for variant_to_try in var_observed:
                if (
                    var_observed.count(variant_to_try) == int(cnvcall[-1]) - 2
                    and variant_to_try not in variant_tried
                ):
                    variant_tried.append(variant_to_try)
                    var_observed_new = var_observed + [variant_to_try]
                    matchtag_new = get_star(var_observed_new, dic)
                    if matchtag_new.call_info == "unique_match":
                        matched_calls.append(matchtag_new)
            if len(matched_calls) == 1:
                matchtag_new = matched_calls[0]
                final_call = matchtag_new.star_call
                final_call_clean = get_final_call_clean(final_call, cnvcall, spacer_cn)
                call_info = matchtag_new.call_info
                raw_call = matchtag_new.candidate
                return star_call(
                    call_info, " ".join(var_observed), raw_call, final_call_clean
                )

        final_call = matchtag.star_call
        final_call_clean = get_final_call_clean(final_call, cnvcall, spacer_cn)
        call_info = matchtag.call_info
        if call_info == "more_than_one_match" and cnvcall == "cn2":
            if sorted(re.split(r"[;/]+", final_call_clean)) == [
                "*1",
                "*27",
                "*32",
                "*41",
            ]:
                if var42126938_G_haplotype:
                    final_call_clean = "*1/*32"
                else:
                    final_call_clean = "*27/*41"
            if sorted(re.split(r"[;/]+", final_call_clean)) == [
                "*1",
                "*119",
                "*2",
                "*41",
            ]:
                if var42127803_diff_haplotype:
                    final_call_clean = "*119/*2"
                else:
                    final_call_clean = "*1/*41"
        raw_call = matchtag.candidate
        return star_call(call_info, " ".join(var_observed), raw_call, final_call_clean)

    else:
        var_observed_68 = var_observed.copy()
        matchtag = call_star68(var_observed, cnvcall, dic)
        final_call = matchtag.star_call
        final_call_clean = get_final_call_clean(final_call, cnvcall, spacer_cn)
        call_info = matchtag.call_info
        raw_call = matchtag.candidate
        return star_call(
            call_info, " ".join(var_observed_68), raw_call, final_call_clean
        )
