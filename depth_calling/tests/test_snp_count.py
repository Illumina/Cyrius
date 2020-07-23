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


import sys
import os
import pytest

from ..snp_count import (
    get_snp_position,
    get_supporting_reads,
    get_supporting_reads_single_region,
    get_fraction,
)
from ..utilities import open_alignment_file

TOTAL_NUM_SITES = 16
test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestParseSNPFile(object):
    def test_parse_snp_file(self):
        snp_file = os.path.join(test_data_dir, "SMN_SNP_37.txt")
        dsnp1, dsnp2, nchr, dindex = get_snp_position(snp_file)
        assert len(dsnp1) == TOTAL_NUM_SITES
        assert len(dsnp2) == TOTAL_NUM_SITES
        assert dsnp1["70247773_12"] == "C_T"
        assert dsnp2["69372353_12"] == "C_T"
        assert dsnp1["70247724_11"] == "G_A"
        assert dsnp2["69372304_11"] == "G_A"
        assert nchr == "5"

        snp_file = os.path.join(test_data_dir, "SMN_SNP_19.txt")
        dsnp1, dsnp2, nchr, dindex = get_snp_position(snp_file)
        assert len(dsnp1) == TOTAL_NUM_SITES
        assert len(dsnp2) == TOTAL_NUM_SITES
        assert dsnp1["70247773_12"] == "C_T"
        assert dsnp2["69372353_12"] == "C_T"
        assert dsnp1["70247724_11"] == "G_A"
        assert dsnp2["69372304_11"] == "G_A"
        assert nchr == "chr5"

        snp_file = os.path.join(test_data_dir, "SMN_SNP_38.txt")
        dsnp1, dsnp2, nchr, dindex = get_snp_position(snp_file)
        assert len(dsnp1) == TOTAL_NUM_SITES
        assert len(dsnp2) == TOTAL_NUM_SITES
        assert dsnp1["70951946_12"] == "C_T"
        assert dsnp2["70076526_12"] == "C_T"
        assert dsnp1["70951463_10"] == "T_C"
        assert dsnp2["70076043_10"] == "T_C"
        assert nchr == "chr5"

        # test indels and reverse complement
        snp_file = os.path.join(test_data_dir, "SMN_SNP_37_test.txt")
        dsnp1, dsnp2, nchr, dindex = get_snp_position(snp_file)
        assert len(dsnp1) == TOTAL_NUM_SITES
        assert len(dsnp2) == TOTAL_NUM_SITES
        assert dsnp1["70245876_1"] == "T_G"
        assert dsnp2["69370451_1"] == "A_C"
        assert dsnp1["70246016_2"] == "G_T"
        assert dsnp2["69370591_2"] == "C_A"
        assert dsnp1["70248108_15"] == "CAC_CC"
        assert dsnp2["69372688_15"] == "CAC_CC"


class TestReadCount(object):
    def test_get_snp_count(self):
        snp_file = os.path.join(test_data_dir, "SMN_SNP_37.txt")
        dsnp1, dsnp2, nchr, dindex = get_snp_position(snp_file)

        bam1 = os.path.join(test_data_dir, "NA12878.bam")
        bamfile1 = open_alignment_file(bam1)
        lsnp1, lsnp2 = get_supporting_reads(bamfile1, dsnp1, dsnp2, nchr, dindex)
        assert lsnp1 == [0, 0, 0, 0, 0, 0, 29, 35, 26, 39, 29, 35, 32, 37, 39, 39]
        assert lsnp2 == [0, 0, 0, 0, 0, 0, 12, 39, 39, 32, 26, 55, 45, 33, 42, 18]

        bam2 = os.path.join(test_data_dir, "NA12885.bam")
        bamfile2 = open_alignment_file(bam2)
        lsnp1, lsnp2 = get_supporting_reads(bamfile2, dsnp1, dsnp2, nchr, dindex)
        assert lsnp1 == [46, 32, 45, 36, 34, 14, 36, 54, 38, 34, 41, 41, 40, 51, 40, 37]
        assert lsnp2 == [35, 35, 32, 29, 35, 59, 22, 28, 32, 24, 34, 32, 33, 28, 38, 21]

        lsnp1, lsnp2, forward, reverse = get_supporting_reads_single_region(
            bamfile2, dsnp1, nchr, dindex
        )
        assert lsnp1 == [46, 32, 45, 36, 26, 14, 36, 54, 38, 34, 41, 41, 40, 51, 40, 34]
        assert lsnp2 == [0, 1, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        # test indels and reverse complement
        snp_file = os.path.join(test_data_dir, "SMN_SNP_37_test.txt")
        dsnp1, dsnp2, nchr, dindex = get_snp_position(snp_file)
        lsnp1, lsnp2, forward, reverse = get_supporting_reads_single_region(
            bamfile2, dsnp1, nchr, dindex
        )
        assert lsnp1 == [46, 32, 45, 36, 26, 14, 36, 54, 38, 34, 41, 41, 40, 51, 40, 19]
        assert lsnp2 == [0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16]

    def test_get_fraction(self):
        lsnp1 = [16, 15, 32, 25, 28, 0]
        lsnp2 = [40, 45, 31, 30, 27, 0]
        smn1_fraction = get_fraction(lsnp1, lsnp2)
        assert smn1_fraction == [16 / 56, 15 / 60, 32 / 63, 25 / 55, 28 / 55, 0]
