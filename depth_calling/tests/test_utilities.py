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

from ..utilities import parse_region_file

test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")


class TestUtilities(object):
    def test_parse_reigon_file(self):
        region_file = os.path.join(test_data_dir, "SMN_region_19_short.bed")
        region_dic = parse_region_file(region_file)
        assert len(region_dic["norm"]) == 500
        assert len(region_dic["exon16"]) == 2
        assert len(region_dic["exon78"]) == 2
