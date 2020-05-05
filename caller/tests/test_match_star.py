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


from ..match_star_allele import check_name, get_final_call_clean


class TestMatchStar(object):
    def test_check_name(self):
        var_called = ["*1_*2"]
        main_allele = check_name(var_called)
        assert list(main_allele) == ["*1_*2"]

        var_called = ["*1_*4A"]
        main_allele = check_name(var_called)
        assert list(main_allele) == ["*1_*4"]

        var_called = ["*1_*4A", "*1_*2"]
        main_allele = check_name(var_called)
        assert len(main_allele) == 2
        assert "*1_*2" in main_allele
        assert "*1_*4" in main_allele

        var_called = ["*1_*4A", "*1_*4D"]
        main_allele = check_name(var_called)
        assert list(main_allele) == ["*1_*4"]

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
        assert clean_call == "*5/*10+*36"

        spacer_cn = 4
        final_call = ["*10_*10"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*5/*10+*36"

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
        assert clean_call == "*2+*13/*1"

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
        assert clean_call == "*2/*2x3;*2x2/*2x2"

        final_call = ["*2_*2_*4_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*2x2/*4x2"

        final_call = ["*2_*2_*2_*4"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*2x3/*4"

        cnvcall = "exon9hyb"
        final_call = ["*10_*10_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*10/*10+*36"

        final_call = ["*1_*4A_*4N"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*4A+*4N"

        cnvcall = "exon9hyb_exon9hyb"
        final_call = ["*10_*10_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*10+*36/*10+*36"

        final_call = ["*1_*10_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*10+*36+*36"

        cnvcall = "exon9hyb_exon9hyb_exon9hyb"
        final_call = ["*10_*10_*36_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*10+*36/*10+*36+*36"

        final_call = ["*1_*10_*36_*36_*36"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*1/*10+*36+*36+*36"

        cnvcall = "star5_star68"
        final_call = ["*4A"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*5/*4A+*68"

        final_call = ["*10"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*5/*10+*68;*68/*10"

        cnvcall = "star68"
        final_call = ["*4A_*40"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*40/*4A+*68"

        final_call = ["*10_*40"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*10_*40_*68"

        cnvcall = "star68_star68"
        final_call = ["*4A_*4A"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*4A+*68/*4A+*68"

        final_call = ["*10_*40"]
        clean_call = get_final_call_clean(final_call, cnvcall, spacer_cn)
        assert clean_call == "*10_*40_*68_*68"
