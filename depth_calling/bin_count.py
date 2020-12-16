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


from collections import namedtuple, OrderedDict
import multiprocessing as mp
from functools import partial
import pysam
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from .utilities import open_alignment_file


MAD_CONSTANT = 1.4826


def get_normed_depth(bamf, region_dic, nCores=1, reference=None, gc_correct=True):
    """
    Return the normalized depth values and coverage stats for a sample
    given a bam file
    """
    counts_for_normalization, gc_for_normalization, region_type_cn, read_length = count_reads_and_prepare_for_normalization(
        bamf, region_dic, nCores, reference
    )
    normed_depth = normalize(
        counts_for_normalization,
        gc_for_normalization,
        region_type_cn,
        read_length,
        gc_correct,
    )

    return normed_depth


def get_normed_depth_from_count(count_file, region_dic, read_length, gc_correct=True):
    """
    Return the normalized depth values and coverage stats for a sample from
    a count file.
    """
    count_dic = get_count_from_file(count_file)
    counts_for_normalization, gc_for_normalization, region_type_cn = process_counts_and_prepare_for_normalization(
        count_dic, region_dic
    )
    normed_depth = normalize(
        counts_for_normalization,
        gc_for_normalization,
        region_type_cn,
        read_length,
        gc_correct,
    )

    return normed_depth


def median_correction(counts):
    """Return values corrected by median."""
    y_counts = np.array(counts)
    y_counts = y_counts / np.median(y_counts)
    return y_counts


def gc_correction(counts, gc, scale_coefficient=0.9):
    """Return values corrected by GC content."""
    y_counts = np.array(counts)
    x_gc = np.array(gc)
    y_counts = y_counts / np.median(y_counts)
    value_lowess = lowess(y_counts, x_gc, return_sorted=False)
    sample_median = np.median(y_counts)
    gc_corrected = []
    for i in range(len(y_counts)):
        scale_factor = scale_coefficient * min(y_counts[i], 2)
        gc_corrected.append(
            y_counts[i] + scale_factor * (sample_median - value_lowess[i])
        )

    return gc_corrected


def get_read_count(bamfile, region, mapq_cutoff=0):
    """
    Return the number of reads that align to a region.
    Keep duplicate reads.
    Keep unmapped reads with mapped mates.
    """
    reads = bamfile.fetch(region[0], region[1], region[2])
    nreads = 0
    for read in reads:
        if (
            read.mapq >= mapq_cutoff
            and read.is_secondary == 0
            and read.is_supplementary == 0
            and read.reference_start >= region[1]
            and read.reference_start < region[2]
        ):
            nreads += 1
    return nreads


def mad(list_of_number):
    """Return the median of absolute deviation."""
    med = np.median(list_of_number)
    return MAD_CONSTANT * np.median([abs(a - med) for a in list_of_number])


def normalize(
    counts_for_normalization,
    gc_for_normalization,
    region_type_cn,
    read_length,
    gc_correct,
):
    """
    Return the normalized depth values for a sample.
    Median normalization and/or GC normalization
    """
    gc_corrected_depth = gc_correction(counts_for_normalization, gc_for_normalization)
    if gc_correct is True:
        # GC normalization
        corrected_depth = gc_corrected_depth
    else:
        # Median normalization
        corrected_depth = median_correction(counts_for_normalization)

    vmedian = np.median(counts_for_normalization) * read_length
    vmad = round(mad(gc_corrected_depth) / np.median(gc_corrected_depth), 3)

    # Also return median depth of this sample and the coverage MAD.
    norm_count = {}
    for i, (region_type, hap_cn) in enumerate(region_type_cn.items()):
        if vmedian == 0:
            norm_count.setdefault(region_type, None)
        else:
            norm_count.setdefault(region_type, 2 * hap_cn * corrected_depth[i])
    depth_value = namedtuple("depth_value", "normalized mediandepth mad")
    normalized_bin = depth_value(norm_count, vmedian, vmad)
    return normalized_bin


def get_read_length(reads, number_to_count=2000):
    """Return the median read length from a set of reads."""
    counter = 0
    read_length = []
    for read in reads:
        counter += 1
        read_length.append(read.query_length)
        if counter > number_to_count:
            break
    return np.median(read_length)


def count_reads_and_prepare_for_normalization(
    bamf, region_dic, nCores=1, reference=None
):
    """
    Return the normalized depth values and coverage stats for a sample
    given a bam file
    """
    bamfile = open_alignment_file(bamf, reference)
    # Store read counts in each interval
    counts_for_normalization = []
    gc_for_normalization = []
    region_type_cn = OrderedDict()
    for region_type in region_dic:
        if region_type != "norm":
            lcount = []
            region_length = None
            hap_cn = None
            for (region, gc) in region_dic[region_type]:
                region_reads = get_read_count(bamfile, region)
                lcount.append(region_reads)
                if "_hapcn" in region[3]:
                    region_length = region[2] - region[1]
                    gc_for_normalization.append(float(gc))
                    hap_cn = int(region[3].split("hapcn")[1])
            if region_length is None or hap_cn is None:
                raise Exception("Problem with region definition. Length not specified.")
            count_sum = sum(lcount)
            counts_for_normalization.append(count_sum / (hap_cn * region_length))
            region_type_cn.setdefault(region_type, hap_cn)

    # Get read length from the last region
    reads = bamfile.fetch(region[0], region[1], region[2])
    read_length = get_read_length(reads)

    lregion = [
        (region[0], region[1], region[2], gc) for (region, gc) in region_dic["norm"]
    ]

    if nCores == 1:
        for region_out in get_normalization_region_values(lregion, bamf, reference):
            counts_for_normalization.append(region_out[0])
            gc_for_normalization.append(region_out[1])
    else:
        get_normed_depth_bam = partial(
            get_normalization_region_values, bam=bamf, reference=reference
        )
        region_groups = partition(lregion, nCores)
        pool = mp.Pool(nCores)
        result = pool.map(get_normed_depth_bam, region_groups)
        pool.close()
        for result_group in result:
            for region_out in result_group:
                counts_for_normalization.append(region_out[0])
                gc_for_normalization.append(region_out[1])

    bamfile.close()

    return counts_for_normalization, gc_for_normalization, region_type_cn, read_length


def partition(lst, n):
    """Partion a list"""
    return [lst[i::n] for i in range(n)]


def get_normalization_region_values(l, bam, reference=None):
    """Perform read counting in a list of regions."""
    bamfile = open_alignment_file(bam, reference)
    lcount = []
    for region in l:
        num_reads = get_read_count(bamfile, region)
        region_length = int(region[2]) - int(region[1])
        norm_depth = num_reads / region_length
        region_gc = float(region[-1])
        lcount.append((norm_depth, region_gc))
    bamfile.close()
    return lcount


def get_count_from_file(count_file):
    """Parse count file"""
    count_dic = {}
    with open(count_file) as f:
        for line in f:
            at = line.strip().split()
            count_dic.setdefault(at[3], int(at[-1]))
    return count_dic


def process_counts_and_prepare_for_normalization(count_dic, region_dic):
    """
    Return the normalized depth values and coverage stats for a sample from
    a count dictionary.
    """
    counts_for_normalization = []
    gc_for_normalization = []
    region_type_cn = OrderedDict()
    for region_type in region_dic:
        if region_type != "norm":
            lcount = []
            region_length = None
            hap_cn = None
            for (region, gc) in region_dic[region_type]:
                lcount.append(count_dic[region[3]])
                if "_hapcn" in region[3]:
                    region_length = region[2] - region[1]
                    gc_for_normalization.append(float(gc))
                    hap_cn = float(region[3].split("hapcn")[1])
            if region_length is None or hap_cn is None:
                raise Exception("Problem with region definition. Length not specified.")
            count_sum = sum(lcount)
            counts_for_normalization.append(count_sum / (hap_cn * region_length))
            region_type_cn.setdefault(region_type, hap_cn)

    for (region, gc) in region_dic["norm"]:
        region_length = int(region[2]) - int(region[1])
        counts_for_normalization.append(count_dic[region[3]] / region_length)
        gc_for_normalization.append(gc)

    return counts_for_normalization, gc_for_normalization, region_type_cn
