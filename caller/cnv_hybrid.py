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

import numpy as np
from collections import Counter, namedtuple


REP_END_POSITION = 3  # end of REP sites
EXON9_END_POSITION = 9  # exon9 site
INTRON4_BP_POSITION = 40  # intron4 bp site
INTRON1_BP_POSITION = 74  # intron1 bp site
REP_SITES_CN = 2  # CN in the REP sites
EXON9_TO_INTRON4_SITES_MIN = 18  # minimum number of sites to make a consensus
INTRON4_TO_INTRON1_SITES_MIN = 19  # minimum number of sites to make a consensus
INTRON1_UPSTREAM_SITES_MIN = 25  # minimum number of sites to make a consensus
EXON9REGION_SITES_MIN = 4  # minimum number of sites to make a consensus
INTRON1_UPSTREAM_SITES_MIN_LOOSE = 15


def get_cnvtag(total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn):
    """
    Return a tag for the called CNV/hybrid group based on detected CN switching point
    at SNP sites between CYP2D6 and CYP2D7.
    The four regions we are looking at are (from left to right, reverse of the direction
    of the gene) REP sites far downstream of the gene, exon9 region sites (including
    sites downstream of but close to the gene), exon9 to intron1 sites and sites
    upstream of intron1.
    """
    exon9region_sites_consensus = None
    exon9_intron4_sites_consensus = None
    intron4_intron1_sites_consensus = None
    intron1_upstream_sites_consensus = None

    exon9_intron4_sites = [
        a
        for a in cn_call_per_site[EXON9_END_POSITION:INTRON4_BP_POSITION]
        if a is not None
    ]
    exon9_intron4_sites_counter = sorted(
        Counter(exon9_intron4_sites).items(), key=lambda kv: kv[1], reverse=True
    )
    if exon9_intron4_sites_counter != []:
        exon9_intron4_sites_consensus = (
            exon9_intron4_sites_counter[0][0]
            if exon9_intron4_sites_counter[0][1] >= EXON9_TO_INTRON4_SITES_MIN
            else None
        )

    intron4_intron1_sites = [
        a
        for a in cn_call_per_site[INTRON4_BP_POSITION:INTRON1_BP_POSITION]
        if a is not None
    ]
    intron4_intron1_sites_counter = sorted(
        Counter(intron4_intron1_sites).items(), key=lambda kv: kv[1], reverse=True
    )

    if intron4_intron1_sites_counter != []:
        intron4_intron1_sites_consensus = (
            intron4_intron1_sites_counter[0][0]
            if intron4_intron1_sites_counter[0][1] >= INTRON4_TO_INTRON1_SITES_MIN
            else None
        )

    if (
        exon9_intron4_sites_consensus is None
        and intron4_intron1_sites_consensus is not None
    ):
        exon9_intron4_sites_consensus = intron4_intron1_sites_consensus
    elif (
        intron4_intron1_sites_consensus is None
        and exon9_intron4_sites_consensus is not None
    ):
        intron4_intron1_sites_consensus = exon9_intron4_sites_consensus

    intron1_upstream_sites = [
        a for a in cn_call_per_site[INTRON1_BP_POSITION:] if a is not None
    ]
    intron1_upstream_sites_counter = sorted(
        Counter(intron1_upstream_sites).items(), key=lambda kv: kv[1], reverse=True
    )
    if intron1_upstream_sites_counter != []:
        intron1_upstream_sites_consensus = (
            intron1_upstream_sites_counter[0][0]
            if intron1_upstream_sites_counter[0][1] >= INTRON1_UPSTREAM_SITES_MIN
            else None
        )
    if intron1_upstream_sites_consensus is None:
        if (
            intron1_upstream_sites.count(total_cn - 2)
            >= INTRON1_UPSTREAM_SITES_MIN_LOOSE
        ):
            intron1_upstream_sites_consensus = total_cn - 2

    if spacer_cn is not None:
        # spacer CN indicates the number of copies starting with CYP2D7
        exon9region_sites_consensus = total_cn - spacer_cn
        if (
            exon9gc_call_stringent is not None
            and exon9_intron4_sites_consensus is not None
        ):
            # Use exon9gc_call_stringent as the CN of CYP2D6 in the rare case of *10D,
            # when fusion breakpoint is before sites downstream of CYP2D6 but after exon9.
            if (
                exon9region_sites_consensus < exon9gc_call_stringent
                and exon9gc_call_stringent <= exon9_intron4_sites_consensus
            ):
                exon9region_sites_consensus = exon9gc_call_stringent
            elif (
                exon9region_sites_consensus > exon9gc_call_stringent
                and exon9gc_call_stringent >= exon9_intron4_sites_consensus
            ):
                exon9region_sites_consensus = exon9gc_call_stringent
    else:
        exon9region_sites = [
            a
            for a in cn_call_per_site[REP_END_POSITION:EXON9_END_POSITION]
            if a is not None
        ]
        exon9region_sites_counter = sorted(
            Counter(exon9region_sites).items(), key=lambda kv: kv[1], reverse=True
        )
        if exon9region_sites_counter != []:
            exon9region_sites_consensus = (
                exon9region_sites_counter[0][0]
                if exon9region_sites_counter[0][1] >= EXON9REGION_SITES_MIN
                else None
            )
    if exon9region_sites_consensus is None and exon9gc_call_stringent is not None:
        exon9region_sites_consensus = exon9gc_call_stringent

    cn_regions = namedtuple(
        "cn_regions",
        "rep exon9_and_downstream exon9_to_intron4 intron4_to_intron1 intron1_upstream",
    )
    consensus = cn_regions(
        REP_SITES_CN,
        exon9region_sites_consensus,
        exon9_intron4_sites_consensus,
        intron4_intron1_sites_consensus,
        intron1_upstream_sites_consensus,
    )

    # CNVs that result in an increase in CYP2D6 CN.
    cn_increase = ["dup", "exon9hyb", "star68"]
    # CNVs that result in a decrease in CYP2D6 CN.
    cn_decrease = ["star13intron1", "star13", "star5"]

    change_point = []
    sv_call = None

    # There are only two copies of complete CYP2D7. Assuming no SV in CYP2D7.
    if consensus.intron1_upstream is None or consensus.intron1_upstream != total_cn - 2:
        return (sv_call, consensus)

    # Assign CNV events based on CNs of the different regions.
    if None not in [consensus.intron4_to_intron1, consensus.intron1_upstream]:
        for _ in range(consensus.intron4_to_intron1 - consensus.intron1_upstream):
            change_point.append("star13intron1")
        for _ in range(consensus.intron1_upstream - consensus.intron4_to_intron1):
            change_point.append("star68")
    if None not in [consensus.exon9_to_intron4, consensus.intron4_to_intron1]:
        for _ in range(consensus.exon9_to_intron4 - consensus.intron4_to_intron1):
            change_point.append("star13intron1")
    if None not in [consensus.exon9_and_downstream, consensus.exon9_to_intron4]:
        for _ in range(consensus.exon9_and_downstream - consensus.exon9_to_intron4):
            change_point.append("star13")
        for _ in range(consensus.exon9_to_intron4 - consensus.exon9_and_downstream):
            change_point.append("exon9hyb")
    if None not in [consensus.rep, consensus.exon9_and_downstream]:
        for _ in range(consensus.rep - consensus.exon9_and_downstream):
            change_point.append("star5")
        for _ in range(consensus.exon9_and_downstream - consensus.rep):
            change_point.append("dup")

    if check_cn_match(
        change_point, cn_increase, cn_decrease, consensus.intron1_upstream
    ):
        sv_call = transform_cnvtag("_".join(sorted(change_point)))
        return (sv_call, consensus)

    # no CNV
    if [
        consensus.exon9_to_intron4,
        consensus.intron4_to_intron1,
        consensus.intron1_upstream,
    ] == [2, 2, 2]:
        return ("cn2", consensus)

    return (sv_call, consensus)


def transform_cnvtag(cnvtag):
    """
    Rename some cnv tags for downstream processing.
    """
    split_call = cnvtag.split("_")
    # exon9hyb_star5 and dup_star13 are unlikely to occur together with yet another sv.
    if cnvtag != "exon9hyb_star5":
        while "exon9hyb" in split_call and "star5" in split_call:
            split_call.remove("exon9hyb")
            split_call.remove("star5")
    if cnvtag != "dup_star13":
        while "dup" in split_call and "star13" in split_call:
            split_call.remove("dup")
            split_call.remove("star13")
    if split_call.count("dup") == len(split_call):
        return "cn" + str(len(split_call) + 2)
    if cnvtag == "dup_dup_exon9hyb_star13intron1":
        return "cn4"
    return "_".join(split_call)


def check_cn_match(sv_list, cn_increase, cn_decrease, final_cn):
    """
    Check that the CNV combination produces the right final copy number.
    """
    if sv_list == []:
        return False
    initial_cn = 2
    for sv in sv_list:
        if sv in cn_increase:
            initial_cn += 1
        if sv in cn_decrease:
            initial_cn -= 1
    if initial_cn == final_cn:
        return True
    return False
