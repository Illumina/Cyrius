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


import os
import sys
import argparse
import json
import logging
import datetime
from collections import namedtuple, OrderedDict
import pysam


from depth_calling.snp_count import (
    get_supporting_reads,
    get_supporting_reads_single_region,
    get_fraction,
    get_snp_position,
)
from depth_calling.gmm import Gmm
from depth_calling.utilities import (
    parse_gmm_file,
    parse_region_file,
    open_alignment_file,
)
from depth_calling.bin_count import (
    get_normed_depth,
    get_normed_depth_from_count,
    get_read_length,
)
from caller.call_variants import (
    NOISY_VAR,
    call_cn_snp,
    call_cn_var,
    call_cn_var_homo,
    get_allele_counts_var42128936,
    update_var42128936,
    get_called_variants,
    call_exon9gc,
    call_var42126938,
    call_var42127526_var42127556,
    call_var42127803hap,
)
from caller.cnv_hybrid import get_cnvtag
from caller.construct_star_table import get_hap_table
from caller.match_star_allele import match_star

MAD_THRESHOLD = 0.11
EXON9_SITE1 = 7
EXON9_SITE2 = 8
HIGH_CN_DEPTH_THRESHOLD = 7.5
HAPLOTYPE_VAR = ["g.42126938C>T", "g.42127803C>T", "g.42127526C>T_g.42127556T>C"]
resource_info = namedtuple(
    "resource_info",
    "genome gmm_parameter region_dic snp_db var_db var_homo_db haplotype_db var_list star_combinations",
)
exon9_values = namedtuple(
    "exon9_values", "exon9_cn exon9cn_in_consensus exon9_raw_site1 exon9_raw_site2"
)
# Below are the SV configurations that the caller is able to call
CNV_ACCEPTED = [
    "star5_star5",
    "star13_star13",
    "star13intron1_star13intron1",
    "star5",
    "star13",
    "star13intron1",
    "star5_star5_star68",
    "star5_star68",
    "cn2",
    "exon9hyb_star5",
    "dup_star13",
    "dup_star13intron1",
    "star13_star68",
    "cn3",
    "exon9hyb",
    "star68",
    "cn4",
    "exon9hyb_exon9hyb",
    "star68_star68",
    "dup_exon9hyb",
    "dup_star68",
    "exon9hyb_star68",
    "cn5",
    "exon9hyb_exon9hyb_exon9hyb",
    "star68_star68_star68",
    "cn6",
    "exon9hyb_exon9hyb_exon9hyb_exon9hyb",
    "star68_star68_star68_star68",
]


def load_parameters():
    """Return parameters."""
    parser = argparse.ArgumentParser(
        description="Call CYP2D6 genotypes from a WGS BAM file."
    )
    parser.add_argument(
        "-m",
        "--manifest",
        help="Manifest listing absolute paths to BAM/CRAM files",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--genome",
        help="Reference genome, select from 19, 37, or 38",
        required=True,
    )
    parser.add_argument("-o", "--outDir", help="Output directory", required=True)
    parser.add_argument("-p", "--prefix", help="Prefix to output file", required=True)
    parser.add_argument(
        "-t",
        "--threads",
        help="Optional, number of threads to use. Default is 1",
        type=int,
        required=False,
        default=1,
    )
    parser.add_argument(
        "--countFilePath", help="Optional path to count files", required=False
    )
    parser.add_argument(
        "-r",
        "--reference",
        help="Optional path to reference fasta file for CRAM",
        required=False,
    )

    args = parser.parse_args()
    if args.genome not in ["19", "37", "38"]:
        raise Exception("Genome not recognized. Select from 19, 37, or 38")

    return args


def d6_star_caller(
    bam, call_parameters, threads, count_file=None, reference_fasta=None, index_name=None
):
    """Return CYP2D6 star allele diplotype calls for each sample."""
    d6_call = namedtuple(
        "d6_call",
        "Coverage_MAD Median_depth Total_CN Spacer_CN Total_CN_raw \
        Spacer_CN_raw Variants_called CNV_group Genotype Filter Raw_star_allele \
        Call_info Exon9_CN CNV_consensus d67_snp_call d67_snp_raw \
        Variant_raw_count",
    )
    # 1. Read counting and normalization
    bamfile = open_alignment_file(bam, reference_fasta, index_filename=index_name)
    if count_file is not None:
        reads = bamfile.fetch()
        read_length = get_read_length(reads)
        normalized_depth = get_normed_depth_from_count(
            count_file, call_parameters.region_dic, read_length
        )
    else:
        normalized_depth = get_normed_depth(
            bam, call_parameters.region_dic, threads, reference=reference_fasta
        )

    # no-call after normalizaton
    if normalized_depth.normalized["d67"] is None:
        sample_call = d6_call(
            normalized_depth.mad,
            normalized_depth.mediandepth,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        return sample_call

    # 2. GMM and CN call
    # There are two regions to call CN based on depth: total CYP2D6+CYP2D7, and CYP2D7 spacer region
    cn_call = namedtuple("cn_call", "d67_cn d67_depth spacer_cn spacer_depth")
    gmm_d67 = Gmm()
    gmm_d67.set_gmm_par(call_parameters.gmm_parameter, "d67")
    gcall_d67 = gmm_d67.gmm_call(normalized_depth.normalized["d67"])
    gmm_spacer = Gmm()
    gmm_spacer.set_gmm_par(call_parameters.gmm_parameter, "spacer")
    gcall_spacer = gmm_spacer.gmm_call(normalized_depth.normalized["spacer"])
    high_cn_low_confidence = False
    if gcall_d67.cn is None and gcall_d67.depth_value > HIGH_CN_DEPTH_THRESHOLD:
        high_cn_low_confidence = True
        raw_cn_call = cn_call(
            int(round(gcall_d67.depth_value)),
            gcall_d67.depth_value,
            gcall_spacer.cn,
            gcall_spacer.depth_value,
        )
    else:
        raw_cn_call = cn_call(
            gcall_d67.cn,
            gcall_d67.depth_value,
            gcall_spacer.cn,
            gcall_spacer.depth_value,
        )

    # 3. Get allele counts at D6/D7 SNP (base difference) sites and target variant sites
    # D6/D7 base difference sites. Get read counts at both D6/D7 positions.
    snp_db = call_parameters.snp_db
    snp_d6, snp_d7 = get_supporting_reads(
        bamfile, snp_db.dsnp1, snp_db.dsnp2, snp_db.nchr, snp_db.dindex
    )

    # Variants not in homology regions. Get read counts only at D6 positions.
    var_db = call_parameters.var_db
    var_alt, var_ref, var_alt_forward, var_alt_reverse = get_supporting_reads_single_region(
        bamfile, var_db.dsnp1, var_db.nchr, var_db.dindex
    )
    # Look more carefully for insertions at 42128936 from reads
    var_list = call_parameters.var_list
    ref_read, long_ins_read, short_ins_read = get_allele_counts_var42128936(
        bamfile, call_parameters.genome
    )
    var_alt, var_ref = update_var42128936(
        var_list, var_alt, var_ref, ref_read, long_ins_read, short_ins_read
    )
    # Variants in homology regions. Get read counts at both D6/D7 positions.
    var_homo_db = call_parameters.var_homo_db
    var_homo_alt, var_homo_ref = get_supporting_reads(
        bamfile,
        var_homo_db.dsnp1,
        var_homo_db.dsnp2,
        var_homo_db.nchr,
        var_homo_db.dindex,
    )
    # This ordered dictionary is for final reporting.
    raw_count = OrderedDict()
    non_homology_variant_count = len(var_alt)
    for i in range(len(call_parameters.var_list)):
        if i < non_homology_variant_count:
            if var_list[i] in NOISY_VAR:
                raw_count.setdefault(
                    var_list[i],
                    "%i(%i:%i),%i"
                    % (var_alt[i], var_alt_forward[i], var_alt_reverse[i], var_ref[i]),
                )
            else:
                raw_count.setdefault(var_list[i], "%i,%i" % (var_alt[i], var_ref[i]))
        else:
            raw_count.setdefault(
                var_list[i],
                "%i,%i"
                % (
                    var_homo_alt[i - non_homology_variant_count],
                    var_homo_ref[i - non_homology_variant_count],
                ),
            )

    # no-call due to total copy number calling
    if raw_cn_call.d67_cn is None:
        sample_call = d6_call(
            normalized_depth.mad,
            normalized_depth.mediandepth,
            raw_cn_call.d67_cn,
            raw_cn_call.spacer_cn,
            raw_cn_call.d67_depth,
            raw_cn_call.spacer_depth,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            raw_count,
        )
        return sample_call

    # 4. Call CNV and hybrids
    d6_fraction = get_fraction(snp_d6, snp_d7)
    raw_d6_cn = [round(raw_cn_call.d67_cn * a, 3) for a in d6_fraction]
    cn_call_snp = call_cn_snp(raw_cn_call.d67_cn, snp_d6, snp_d7)

    # exon9gc
    exon9gc_call_stringent = call_exon9gc(
        snp_d6[EXON9_SITE1 : EXON9_SITE2 + 1],
        snp_d7[EXON9_SITE1 : EXON9_SITE2 + 1],
        raw_cn_call.d67_cn,
    )
    cnvtag, consensus = get_cnvtag(
        raw_cn_call.d67_cn,
        raw_d6_cn,
        cn_call_snp,
        exon9gc_call_stringent,
        raw_cn_call.spacer_cn,
    )

    # no-call due to CNV group calling
    if cnvtag is None or cnvtag not in CNV_ACCEPTED:
        sample_call = d6_call(
            normalized_depth.mad,
            normalized_depth.mediandepth,
            raw_cn_call.d67_cn,
            raw_cn_call.spacer_cn,
            raw_cn_call.d67_depth,
            raw_cn_call.spacer_depth,
            None,
            cnvtag,
            None,
            None,
            None,
            None,
            exon9gc_call_stringent,
            ",".join(str(a) for a in consensus),
            ",".join(str(a) for a in cn_call_snp),
            ",".join(str(a) for a in raw_d6_cn),
            raw_count,
        )
        return sample_call

    # 5. Call variants
    # homology region
    cn_call_var_homo = call_cn_var_homo(raw_cn_call.d67_cn, var_homo_alt, var_homo_ref)
    # non-homology region
    cn_call_var = call_cn_var(
        cnvtag, var_alt, var_ref, var_alt_forward, var_alt_reverse, var_list, var_db
    )
    # call haplotypes
    haplotype_db = call_parameters.haplotype_db
    site42126938_count, var42126938, var42126938_G_haplotype = call_var42126938(
        bamfile, raw_cn_call.d67_cn, haplotype_db["g.42126938C>T"]
    )
    raw_count.setdefault(
        "g.42126938C>T", "%i,%i" % (site42126938_count[1], site42126938_count[0])
    )

    site42127526_count, site42127556_count, var42127526 = call_var42127526_var42127556(
        bamfile, cnvtag, haplotype_db["g.42127526C>T_g.42127556T>C"]
    )
    raw_count.setdefault(
        "g.42127526C>T", "%i,%i" % (site42127526_count[1], site42127526_count[0])
    )
    raw_count.setdefault(
        "g.42127556T>C", "%i,%i" % (site42127556_count[1], site42127556_count[0])
    )

    var42127803_diff_haplotype = call_var42127803hap(
        bamfile, cnvtag, haplotype_db["g.42127803C>T"]
    )

    # 6. Call star allele
    total_callset = get_called_variants(var_list, cn_call_var)
    called_var_homo = get_called_variants(var_list, cn_call_var_homo, len(cn_call_var))
    total_callset += called_var_homo
    total_callset += var42126938
    total_callset += var42127526

    star_called = match_star(
        total_callset,
        cnvtag,
        raw_cn_call.spacer_cn,
        call_parameters.star_combinations,
        exon9_values(
            exon9gc_call_stringent,
            consensus.exon9_and_downstream,
            raw_d6_cn[EXON9_SITE1],
            raw_d6_cn[EXON9_SITE2],
        ),
        var42126938_G_haplotype,
        var42127803_diff_haplotype,
    )

    genotype_filter = None
    final_star_allele_call = None
    # no-call due to star allele matching
    if star_called.call_info and star_called.call_info != "no_match":
        final_star_allele_call = star_called.clean_call
        if final_star_allele_call:
            if ";" in final_star_allele_call:
                genotype_filter = "More_than_one_possible_genotype"
            elif "/" not in final_star_allele_call:
                genotype_filter = "Not_assigned_to_haplotypes"
            elif high_cn_low_confidence:
                genotype_filter = "LowQ_high_CN"
            else:
                genotype_filter = "PASS"

    sample_call = d6_call(
        normalized_depth.mad,
        normalized_depth.mediandepth,
        raw_cn_call.d67_cn,
        raw_cn_call.spacer_cn,
        raw_cn_call.d67_depth,
        raw_cn_call.spacer_depth,
        star_called.variants_called.split(),
        cnvtag,
        final_star_allele_call,
        genotype_filter,
        star_called.raw_call,
        star_called.call_info,
        exon9gc_call_stringent,
        ",".join(str(a) for a in consensus),
        ",".join(str(a) for a in cn_call_snp),
        ",".join(str(a) for a in raw_d6_cn),
        raw_count,
    )
    bamfile.close()
    return sample_call


def prepare_resource(datadir, parameters):
    genome = parameters.genome
    region_file = os.path.join(datadir, "CYP2D6_region_%s.bed" % genome)
    snp_file = os.path.join(datadir, "CYP2D6_SNP_%s.txt" % genome)
    gmm_file = os.path.join(datadir, "CYP2D6_gmm.txt")
    star_table = os.path.join(datadir, "star_table.txt")
    variant_file = os.path.join(datadir, "CYP2D6_target_variant_%s.txt" % genome)
    variant_homology_file = os.path.join(
        datadir, "CYP2D6_target_variant_homology_region_%s.txt" % genome
    )
    haplotype_file = os.path.join(datadir, "CYP2D6_haplotype_%s.txt" % genome)
    star_combinations = get_hap_table(star_table)

    for required_file in [
        region_file,
        snp_file,
        variant_file,
        variant_homology_file,
        haplotype_file,
        gmm_file,
    ]:
        if os.path.exists(required_file) == 0:
            raise Exception("File %s not found." % required_file)

    snp_db = get_snp_position(snp_file)
    var_db = get_snp_position(variant_file)
    var_homo_db = get_snp_position(variant_homology_file)
    haplotype_db = {}
    for variant in HAPLOTYPE_VAR:
        haplotype_db.setdefault(variant, get_snp_position(haplotype_file, variant))
    var_list = []
    with open(variant_file) as f:
        for line in f:
            if line[0] != "#":
                var_name = line.split()[-1]
                var_list.append(var_name)
    with open(variant_homology_file) as f:
        for line in f:
            if line[0] != "#":
                var_name = line.split()[-1]
                var_list.append(var_name)
    gmm_parameter = parse_gmm_file(gmm_file)
    region_dic = parse_region_file(region_file)
    call_parameters = resource_info(
        genome,
        gmm_parameter,
        region_dic,
        snp_db,
        var_db,
        var_homo_db,
        haplotype_db,
        var_list,
        star_combinations,
    )
    return call_parameters


def main():
    parameters = load_parameters()
    manifest = parameters.manifest
    outdir = parameters.outDir
    prefix = parameters.prefix
    reference_fasta = parameters.reference
    threads = parameters.threads
    path_count_file = parameters.countFilePath
    logging.basicConfig(level=logging.DEBUG)

    if os.path.exists(outdir) == 0:
        os.makedirs(outdir)

    # Prepare data files
    datadir = os.path.join(os.path.dirname(__file__), "data")
    call_parameters = prepare_resource(datadir, parameters)

    out_json = os.path.join(outdir, prefix + ".json")
    out_tsv = os.path.join(outdir, prefix + ".tsv")
    final_output = {}
    with open(manifest) as read_manifest:
        for line in read_manifest:
            bam_name = line.strip()
            index_name = None
            if '##idx##' in bam_name:
                bam_name, index_name = bam_name.split('##idx##')

            sample_id = os.path.splitext(os.path.basename(bam_name))[0]
            count_file = None
            if path_count_file is not None:
                count_file = os.path.join(path_count_file, sample_id + "_count.txt")
            if "://" not in bam_name and os.path.exists(bam_name) == 0:
                logging.warning("Input file for sample %s does not exist.", sample_id)
            else:
                logging.info(
                    "Processing sample %s at %s", sample_id, datetime.datetime.now()
                )
                cyp2d6_call = d6_star_caller(
                    bam_name, call_parameters, threads, count_file, reference_fasta, index_name=index_name
                )._asdict()
                # Use normalized coverage MAD across stable regions
                # as a sample QC measure.
                if cyp2d6_call["Coverage_MAD"] > MAD_THRESHOLD:
                    logging.warning(
                        "Sample %s has uneven coverage. CN calls may be unreliable.",
                        sample_id,
                    )
                final_output.setdefault(sample_id, cyp2d6_call)

    # Write to json
    logging.info("Writing to json at %s", datetime.datetime.now())
    with open(out_json, "w") as json_output:
        json.dump(final_output, json_output)

    # Write to tsv
    logging.info("Writing to tsv at %s", datetime.datetime.now())
    header = ["Sample", "Genotype", "Filter"]
    with open(out_tsv, "w") as tsv_output:
        tsv_output.write("\t".join(header) + "\n")
        for sample_id in final_output:
            final_call = final_output[sample_id]
            output_per_sample = [
                sample_id,
                final_call["Genotype"],
                final_call["Filter"],
            ]
            tsv_output.write("\t".join(str(a) for a in output_per_sample) + "\n")


if __name__ == "__main__":
    main()
