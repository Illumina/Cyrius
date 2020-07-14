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
from collections import namedtuple
from ..cnv_hybrid import get_cnvtag


TOTAL_SITE = 117
cn_regions = namedtuple(
    "cn_regions",
    "rep exon9_and_downstream exon9_to_intron4 intron4_to_intron1 intron1_upstream",
)


class TestCallCNVGroup(object):
    def test_get_cnvtag(self):
        total_cn = 5
        cn_call_per_site = [2] * 3 + [3] * (TOTAL_SITE - 3)
        rawv = cn_call_per_site
        exon9gc_call_stringent = 3
        spacer_cn = 2
        cnvtag = get_cnvtag(
            total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn
        )
        assert cnvtag[0] == "cn3"

        cn_call_per_site = [2] * 9 + [3] * (TOTAL_SITE - 9)
        rawv = cn_call_per_site
        exon9gc_call_stringent = 2
        spacer_cn = 3
        cnvtag = get_cnvtag(
            total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn
        )
        assert cnvtag[0] == "exon9hyb"
        assert cnvtag[1] == cn_regions(2, 2, 3, 3, 3)

        # *10D
        cn_call_per_site = [2] * 7 + [3] * (TOTAL_SITE - 7)
        rawv = cn_call_per_site
        exon9gc_call_stringent = 3
        spacer_cn = 3
        cnvtag = get_cnvtag(
            total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn
        )
        assert cnvtag[0] == "cn3"
        assert cnvtag[1] == cn_regions(2, 3, 3, 3, 3)

        # conversion to CYP2D7 downstream of the gene
        cn_call_per_site = [2] * 3 + [1] * 4 + [2, 2] + [3] * (TOTAL_SITE - 9)
        rawv = cn_call_per_site
        exon9gc_call_stringent = 2
        spacer_cn = 4
        cnvtag = get_cnvtag(
            total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn
        )
        assert cnvtag[0] == "exon9hyb"
        assert cnvtag[1] == cn_regions(2, 2, 3, 3, 3)

        cn_call_per_site = [2] * 74 + [3] * (TOTAL_SITE - 74)
        rawv = cn_call_per_site
        exon9gc_call_stringent = 2
        spacer_cn = 3
        cnvtag = get_cnvtag(
            total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn
        )
        assert cnvtag[0] == "star68"

        total_cn = 6
        cn_call_per_site = [2] * 3 + [3] * 71 + [4] * (TOTAL_SITE - 74)
        rawv = cn_call_per_site
        exon9gc_call_stringent = 3
        spacer_cn = 3
        cnvtag = get_cnvtag(
            total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn
        )
        assert cnvtag[0] == "dup_star68"

        total_cn = 4
        cn_call_per_site = [2] * 3 + [1] * 6 + [2] * (TOTAL_SITE - 9)
        rawv = cn_call_per_site
        exon9gc_call_stringent = 1
        spacer_cn = 3
        cnvtag = get_cnvtag(
            total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn
        )
        assert cnvtag[0] == "exon9hyb_star5"

        # exon 9 gene conversion
        total_cn = 4
        cn_call_per_site = [2] * 7 + [1] * 2 + [2] * (TOTAL_SITE - 9)
        rawv = cn_call_per_site
        exon9gc_call_stringent = 1
        spacer_cn = 2
        cnvtag = get_cnvtag(
            total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn
        )
        assert cnvtag[0] == "cn2"

        # fusion deletion with breakpoint in intron4
        total_cn = 3
        cn_call_per_site = [2] * 40 + [1] * (TOTAL_SITE - 40)
        rawv = cn_call_per_site
        exon9gc_call_stringent = 2
        spacer_cn = 1
        cnvtag = get_cnvtag(
            total_cn, rawv, cn_call_per_site, exon9gc_call_stringent, spacer_cn
        )
        assert cnvtag[0] == "star13intron1"
        assert cnvtag[1] == cn_regions(2, 2, 2, 1, 1)
