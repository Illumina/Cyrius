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


from scipy.stats import poisson

POSTERIOR_CUTOFF_STRINGENT = 0.9
ERROR_RATE = 0.01


def call_reg1_cn(full_cn, count_reg1, count_reg2, min_read=0):
    """
    Return the reg1 copy number call at each site based on Poisson likelihood,
    with a minimum read support cutoff
    """
    if full_cn is None:
        return [None]
    if full_cn == 0:
        return [0]
    if full_cn == 1 and count_reg1 > min_read:
        return [1]
    prob = []
    nsum = count_reg1 + count_reg2
    if nsum == 0:
        return [None]
    for i in range(full_cn + 1):
        depthexpected = float(nsum) * float(i) / float(full_cn)
        if i == 0:
            depthexpected = (ERROR_RATE / 3) * float(nsum)
        if i == full_cn:
            depthexpected = float(nsum) - ERROR_RATE * float(nsum)
        if count_reg1 <= count_reg2:
            prob.append(poisson.pmf(int(count_reg1), depthexpected))
        else:
            prob.append(poisson.pmf(int(count_reg2), depthexpected))
    sum_prob = sum(prob)
    if sum_prob == 0:
        return [None]
    post_prob = [float(a) / float(sum_prob) for a in prob]
    if count_reg2 < count_reg1:
        post_prob = post_prob[::-1]
    post_prob_sorted = sorted(post_prob, reverse=True)
    if (
        post_prob.index(post_prob_sorted[0]) != 0
        and count_reg1 <= min_read
        and count_reg2 >= min_read
    ):
        return [0]
    if post_prob_sorted[0] >= POSTERIOR_CUTOFF_STRINGENT:
        return [post_prob.index(post_prob_sorted[0])]
    # output the two most likely scenarios
    cn_prob_filtered = [
        post_prob.index(post_prob_sorted[0]),
        round(post_prob_sorted[0], 3),
        post_prob.index(post_prob_sorted[1]),
        round(post_prob_sorted[1], 3),
    ]
    return cn_prob_filtered


def process_raw_call_gc(cn_prob, post_cutoff, keep_none=True):
    """
    Filter raw CN calls based on posterior probablity cutoff.
    For gene conversion cases, i.e. SNVs between paralogs.
    """
    cn_prob_filtered = []
    for cn_call in cn_prob:
        if len(cn_call) == 1:
            call_value = cn_call[0]
            if call_value is not None or keep_none:
                cn_prob_filtered.append(call_value)
        elif cn_call[1] > post_cutoff:
            cn_prob_filtered.append(cn_call[0])
        elif keep_none:
            cn_prob_filtered.append(None)
    return cn_prob_filtered


def process_raw_call_denovo(
    cn_prob, post_cutoff1, post_cutoff2, list_total_cn=None, keep_none=True
):
    """
    Filter raw CN calls based on posterior probablity cutoff.
    For de novel variant calling, i.e. non-gene-conversion cases.
    For less confident calls that are not copy number zero,
    return the smaller CN call.
    Also keep the variant if called CN is equal to total CN at the site.
    This makes sure a variant is always kept if CN>=1
    but we can't distinguish 1 vs 2.
    """
    cn_prob_filtered = []
    for i, cn_call in enumerate(cn_prob):
        if len(cn_call) == 1:
            call_value = cn_call[0]
            if call_value is not None or keep_none:
                cn_prob_filtered.append(call_value)
        else:
            if list_total_cn is not None:
                total_cn = list_total_cn[i]
                keep_var = (cn_call[0] > 0 and cn_call[2] > 0) or (
                    cn_call[0] == total_cn or cn_call[2] == total_cn
                )
            else:
                keep_var = cn_call[0] > 0 and cn_call[2] > 0
            if cn_call[1] > post_cutoff1:
                cn_prob_filtered.append(cn_call[0])
            elif keep_var:
                if cn_call[1] > post_cutoff2:
                    cn_prob_filtered.append(cn_call[0])
                else:
                    cn_prob_filtered.append(min(cn_call[0], cn_call[2]))
            elif keep_none:
                cn_prob_filtered.append(None)
    return cn_prob_filtered
