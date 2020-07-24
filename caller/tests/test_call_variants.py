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
import pysam


from ..call_variants import (
    process_raw_call_gc,
    process_raw_call_denovo,
    get_allele_counts_var42128936,
    call_exon9gc,
    get_called_variants,
)


TOTAL_SITE = 118
test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestCallCN(object):
    def test_call_42128936(self):
        bam = pysam.AlignmentFile(os.path.join(test_data_dir, "NA23275.bam"), "rb")
        ref_read, long_ins_read, short_ins_read = get_allele_counts_var42128936(
            bam, "37"
        )
        assert long_ins_read == 6

    def test_get_called_variants(self):
        var_list = ["var1", "var2", "var3", "var4"]
        cn_called = [0, None, 1, 2]
        var_called = get_called_variants(var_list, cn_called)
        assert var_called == ["var3", "var4", "var4"]

        var_list = ["var1", "var2", "var3", "var4", "var5"]
        cn_called = [0, None, 1, 2]
        var_called = get_called_variants(var_list, cn_called, 1)
        assert var_called == ["var4", "var5", "var5"]
