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


from .snp_count import passing_read


def get_haplotypes_from_bam(bamfile_handle, base_db, target_positions):
    """
    Get haplotypes for each read at a set of target positions
    """
    dread = {}
    dread = get_bases_per_read(bamfile_handle, base_db, target_positions)
    base1, base2 = get_base1_base2(base_db, target_positions)
    dhaplotype = get_hap_counts(dread, base1, base2)
    return dhaplotype


def get_haplotypes_from_bam_single_region(bamfile_handle, base_db, target_positions):
    dread = {}
    dread = get_bases_per_read(
        bamfile_handle, base_db, target_positions, region=0, min_mapq=10
    )
    base1, base2 = get_base1_base2(base_db, target_positions)
    dhaplotype = get_hap_counts(dread, base1, base2)
    return dhaplotype


def get_bases_per_read(
    bamfile_handle, base_db, target_positions, region=None, min_mapq=0
):
    dread = {}
    nchr = base_db.nchr
    dindex = base_db.dindex
    dsnps = [base_db.dsnp1, base_db.dsnp2]
    if region is not None:
        if region == 0:
            dsnps = [base_db.dsnp1]
        elif region == 1:
            dsnps = [base_db.dsnp2]
    for dsnp in dsnps:
        for snp_position_ori in dsnp:
            dsnp_index = dindex[snp_position_ori]
            snp_position = int(snp_position_ori.split("_")[0])
            if dsnp_index in target_positions:
                for pileupcolumn in bamfile_handle.pileup(
                    nchr,
                    snp_position - 1,
                    snp_position + 1,
                    truncate=True,
                    stepper="nofilter",
                    ignore_overlaps=False,
                    ignore_orphan=False,
                ):
                    site_position = pileupcolumn.pos + 1
                    if site_position == snp_position:
                        reg1_allele, reg2_allele = dsnp[snp_position_ori].split("_")
                        for read in pileupcolumn.pileups:
                            if (
                                passing_read(read)
                                and read.alignment.mapping_quality >= min_mapq
                            ):
                                read_name = read.alignment.query_name
                                read_seq = read.alignment.query_sequence
                                start_pos = read.query_position
                                end_pos = start_pos + min(
                                    len(reg1_allele), len(reg2_allele)
                                )
                                if end_pos < len(read_seq):
                                    hap = read_seq[start_pos:end_pos]
                                    if read_name not in dread:
                                        dread.setdefault(read_name, {})
                                        for pos in target_positions:
                                            dread[read_name].setdefault(pos, None)
                                    if dread[read_name][dsnp_index] not in [None, hap]:
                                        dread[read_name][dsnp_index] = None
                                    dread[read_name][dsnp_index] = hap
    return dread


def get_hap_counts(dread, base1, base2):
    """
    Translate bases into haplotypes
    """
    dread_hap = {}
    for read in dread:
        read_bases = dread[read]
        assert len(read_bases) == len(base1) == len(base2)
        pos_list = ["x"] * len(read_bases)
        for i, pos in enumerate(read_bases):
            base = read_bases[pos]
            if base is not None:
                for allele in base1[i].split(","):
                    if base == allele.upper():
                        pos_list[i] = "1"
                for allele in base2[i].split(","):
                    if base == allele.upper():
                        pos_list[i] = "2"
        dread_hap.setdefault(read, "".join(pos_list))
    return dread_hap


def get_base1_base2(base_db, target_positions):
    """
    Get expected bases corresponding to different alleles/paralogs
    at target positions 
    """
    base1 = [None] * len(base_db.dsnp1)
    base2 = [None] * len(base_db.dsnp1)
    dindex = base_db.dindex
    for pos in base_db.dsnp1:
        dsnp_index = dindex[pos]
        if dsnp_index in target_positions:
            index = int(pos.split("_")[1])
            allele1, allele2 = base_db.dsnp1[pos].split("_")
            base1[index] = allele1
            base2[index] = allele2
    return [base1[i] for i in target_positions], [base2[i] for i in target_positions]


def extract_hap(dhaplotype, positions_to_extract):
    """
    Extract haplotypes at certain positions
    """
    hap_count = {}
    for read_name in dhaplotype:
        hap = dhaplotype[read_name]
        hap_base = [hap[pos] for pos in positions_to_extract]
        if "x" not in hap_base:
            hap_count.setdefault("".join(hap_base), []).append(1)
    return hap_count
