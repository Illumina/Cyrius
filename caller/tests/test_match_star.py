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


from ..construct_star_table import get_hap_table
from ..match_star_allele import (
    convert_to_main_allele,
    get_final_call_clean,
    CNVTAG_TO_GENOTYPE,
    get_dic,
)

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from star_caller import CNV_ACCEPTED


class TestMatchStar(object):
    def test_accepted_cnv(self):
        star_table = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
            "data",
            "star_table.txt",
        )
        star_combinations = get_hap_table(star_table)
        for cnvtag in CNV_ACCEPTED:
            if cnvtag not in CNVTAG_TO_GENOTYPE:
                dic = get_dic(cnvtag, star_combinations)
                assert dic is not None

    def test_check_name(self):
        var_called = ["*1_*2"]
        main_allele = convert_to_main_allele(var_called)
        assert main_allele == ["*1_*2"]

        var_called = ["*1_*4"]
        main_allele = convert_to_main_allele(var_called)
        assert main_allele == ["*1_*4"]

        var_called = ["*1_*4", "*1_*2"]
        main_allele = convert_to_main_allele(var_called)
        assert len(main_allele) == 2
        assert "*1_*2" in main_allele
        assert "*1_*4" in main_allele

        var_called = ["*1_*4", "*1_*4.009"]
        main_allele = convert_to_main_allele(var_called)
        assert main_allele == ["*1_*4"]

    def test_clean_call(self):
        cnvcall = "star5"
        spacer_cn = None

        final_call = ["*1"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*5"

        cnvcall = "cn2"
        spacer_cn = None

        final_call = []
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == ""

        final_call = ["*1_*2", "*4_*17"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*2;*4/*17"

        final_call = ["*1_*2"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*2"

        cnvcall = "exon9hyb_star5"
        final_call = ["*10_*10"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call is None

        spacer_cn = 2
        final_call = ["*10_*10"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call is None

        spacer_cn = 3
        final_call = ["*10_*10"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*5/*36+*10"

        spacer_cn = 4
        final_call = ["*10_*10"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*5/*36+*10"

        cnvcall = "star13"
        spacer_cn = None
        final_call = ["*2"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*2/*13"

        cnvcall = "star13intron1"
        final_call = ["*1_*2"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*13/*1"

        final_call = ["*1_*1"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call is None

        cnvcall = "dup_star13intron1"
        final_call = ["*1_*2_*2"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*13+*2/*1"

        cnvcall = "dup_star13intron1"
        final_call = ["*9_*4_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*13+*4/*9"

        cnvcall = "dup_star13intron1"
        final_call = ["*4_*4_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*13+*4/*4"

        cnvcall = "dup_star13"
        final_call = ["*1_*2"]
        spacer_cn = None
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*2"
        spacer_cn = 1
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1_*2_*13"

        cnvcall = "cn3"
        final_call = ["*2_*2_*2"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*2/*2x2"

        final_call = ["*2_*2_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*2x2/*4"

        final_call = ["*1_*2_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1_*2_*4"

        cnvcall = "cn4"
        final_call = ["*2_*2_*2_*2"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*2x2/*2x2"

        final_call = ["*2_*2_*4_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*2x2/*4x2"

        final_call = ["*2_*2_*2_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*2x3/*4"

        cnvcall = "exon9hyb"
        final_call = ["*10_*10_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*10/*36+*10"

        final_call = ["*1_*4_*4.013"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*4.013+*4"

        cnvcall = "exon9hyb_exon9hyb"
        final_call = ["*10_*10_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*36+*10/*36+*10"

        final_call = ["*1_*10_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*36+*36+*10"

        cnvcall = "exon9hyb_exon9hyb_exon9hyb"
        final_call = ["*10_*10_*36_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*36+*10/*36+*36+*10"

        final_call = ["*1_*10_*36_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*36+*36+*36+*10"

        cnvcall = "exon9hyb_exon9hyb_exon9hyb_exon9hyb"
        final_call = ["*10_*10_*36_*36_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*36+*36+*10/*36+*36+*10"

        final_call = ["*1_*10_*36_*36_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*36+*36+*36+*36+*10"

        cnvcall = "star5_star68"
        final_call = ["*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*5/*68+*4"

        final_call = ["*10"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*68/*10"

        cnvcall = "star68"
        final_call = ["*4_*40"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*40/*68+*4"

        final_call = ["*10_*40"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*10_*40_*68"

        cnvcall = "star68_star68"
        final_call = ["*4_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*68+*4/*68+*4"

        final_call = ["*10_*40"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*10_*40_*68_*68"

        cnvcall = "star68_star68_star68"
        final_call = ["*4_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*68+*4/*68+*68+*4"

        final_call = ["*2_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*2/*68+*68+*68+*4"

        final_call = ["*10_*40"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*10_*40_*68_*68_*68"
