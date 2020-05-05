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

from ..copy_number_call import (
    call_reg1_cn,
    process_raw_call_gc,
    process_raw_call_denovo,
)


class TestCallCN(object):
    def test_call_reg1_cn(self):
        call = call_reg1_cn(4, 30, 30)
        assert call == [2]
        call = call_reg1_cn(4, 10, 30)
        assert call == [1]
        call = call_reg1_cn(4, 30, 10)
        assert call == [3]
        call = call_reg1_cn(3, 30, 30)
        assert call == [2, 0.689, 1, 0.311]

    def test_process_raw_call_gc(self):
        lcn = [[1], [1, 0.8, 2, 0.2]]
        filtered_call = process_raw_call_gc(lcn, 0.7)
        assert filtered_call == [1, 1]
        filtered_call = process_raw_call_gc(lcn, 0.9)
        assert filtered_call == [1, None]
        filtered_call = process_raw_call_gc(lcn, 0.9, keep_none=False)
        assert filtered_call == [1]

    def test_process_raw_call_denovo(self):
        lcn = [[1], [0, 0.8, 1, 0.2], [2, 0.6, 1, 0.4]]
        filtered_call = process_raw_call_denovo(lcn, 0.9, 0.7)
        assert filtered_call == [1, None, 1]
        filtered_call = process_raw_call_denovo(lcn, 0.9, 0.55)
        assert filtered_call == [1, None, 2]

        lcn = [[1], [1, 0.65, 0, 0.35], [2, 0.6, 1, 0.4]]
        filtered_call = process_raw_call_denovo(lcn, 0.9, 0.62, [1, 1, 1])
        assert filtered_call == [1, 1, 1]
        filtered_call = process_raw_call_denovo(lcn, 0.9, 0.7, [1, 1, 1])
        assert filtered_call == [1, 0, 1]
