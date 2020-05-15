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

from ..haplotype import get_base1_base2, get_hap_counts, extract_hap
from ..snp_count import get_snp_position

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestUtilities(object):
    def test_get_bases(self):
        snp_file = os.path.join(test_data_dir, "SMN_SNP_37.txt")
        base_db = get_snp_position(snp_file)
        target_positions = [11, 12, 13]
        base1, base2 = get_base1_base2(base_db, target_positions)
        assert base1 == ["G", "C", "A"]
        assert base2 == ["A", "T", "G"]

    def test_get_haplotype(self):
        snp_file = os.path.join(test_data_dir, "SMN_SNP_37.txt")
        base_db = get_snp_position(snp_file)
        target_positions = [11, 12]
        base1, base2 = get_base1_base2(base_db, target_positions)
        assert base1 == ["G", "C"]
        assert base2 == ["A", "T"]
        dread = {
            "read1": {0: "G", 1: "C"},
            "read2": {0: "A", 1: "T"},
            "read3": {0: "G", 1: "T"},
            "read4": {0: "G", 1: "T"},
        }
        dhap = get_hap_counts(dread, base1, base2)
        hap_count = extract_hap(dhap, [0, 1])
        assert hap_count["11"] == [1]
        assert hap_count["22"] == [1]
        assert hap_count["12"] == [1, 1]

        hap_count = extract_hap(dhap, [1])
        assert hap_count["1"] == [1]
        assert hap_count["2"] == [1, 1, 1]
