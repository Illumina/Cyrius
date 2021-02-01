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


from collections import namedtuple
import pysam


def parse_region_file(region_file):
    """Return the set of regions for counting from a bed file."""
    region_dic = {}
    with open(region_file) as read_region:
        for line in read_region:
            nchr, region_start, region_end, region_name, region_type, region_gc = (
                line.strip().split()
            )
            region_start = int(region_start)
            region_end = int(region_end)
            region = (nchr, region_start, region_end, region_name)
            region_dic.setdefault(region_type, []).append((region, region_gc))
    return region_dic


def parse_gmm_file(gmm_file):
    """Return the gmm parameters stored in input file."""
    dpar_tmp = {}
    with open(gmm_file) as read_gmm:
        for line in read_gmm:
            split_line = line.strip().split()
            dpar_tmp.setdefault(split_line[0], {})
            list_value = [a.split(":")[-1] for a in split_line[2:]]
            dpar_tmp[split_line[0]].setdefault(split_line[1], list_value)
    return dpar_tmp


def open_alignment_file(alignment_file, reference_fasta=None, index_filename=None):
    if alignment_file.endswith("cram"):
        return pysam.AlignmentFile(
            alignment_file, "rc", reference_filename=reference_fasta, index_filename=index_filename
        )
    return pysam.AlignmentFile(alignment_file, "rb", index_filename=index_filename)
