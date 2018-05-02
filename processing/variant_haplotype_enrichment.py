__author__ = 'konradk'

import argparse
from utils import *
from slack_utils import *
from pprint import pprint

finland_vds_path = 'gs://konradk/finland/FINRISK_QC_hard.vds'
finland_temp_vds_path = 'gs://konradk/finland/FINRISK_QC_hard_10000.vds'
haplotypes_path = 'gs://konradk/finland/haplotypes/Engage*bgz'
sample_mapping_path = 'gs://konradk/finland/all_birth_wgs_exome_combined.csv'
haplotypes_kt_path = 'gs://konradk/finland/haplotypes/Engagex_haplotypes.kt'
output_path = 'gs://konradk/finland/shared_variants.tsv.bgz'

# Ran this to get more partitions


def main(args):
    hc = HailContext()

    if args.write_haplotypes:
        kt_columns = ['sample1', 'sample1chrom', 'sample2', 'sample2chrom', 'chrom', 'start', 'stop', 'rsid1', 'rsid2', 'something', 'size', 'unit', 'int1', 'int2', 'int3']
        haplotypes_kt = (hc.import_table(haplotypes_path, impute=True, no_header=True)
                         .rename(kt_columns))
        haplotypes_kt = haplotypes_kt.annotate('interval = Interval(str(chrom), start, stop)').key_by('interval')
        haplotypes_kt.write(haplotypes_kt_path, overwrite=args.overwrite)

    # vds = hc.read(finland_vds_path)
    #
    # # Sample renaming and annotation
    # # vds = vds.rename_samples({x: x.replace(' ', '_') for x in vds.sample_ids})  # Translate spaces into underscores (apparently we don't need this for now)
    # samples_kt = hc.import_table(sample_mapping_path, impute=True, delimiter=',').key_by('exome_id')
    # vds = vds.annotate_samples_table(samples_kt, root='sa.meta')
    #
    # vds = vds.filter_samples_expr('isDefined(sa.meta.ID2)')
    # vds = vds.annotate_variants_expr('va.callstats = gs.callStats(g => v)')
    # vds = vds.filter_variants_expr('va.callstats.AC[1] > 1')
    #
    # vds.repartition(10000).write(finland_temp_vds_path, overwrite=True)
    vds = hc.read(finland_temp_vds_path)

    # Check how many samples got a mapping
    pprint(vds.sample_schema)
    vds.query_samples('samples.map(s => isDefined(sa.meta.ID2)).counter()')
    # for gs://konradk/finland/FINRISK_QC_hard.vds, this is: {False: 4673L, True: 4734L}

    haplotypes_kt = (hc.read_table(haplotypes_kt_path)
                     .filter('sample1 != sample2')
                     .annotate('sample_pair = [sample1, sample2].toSet()')
                     .select(['sample_pair', 'interval'])
    )
    # pprint(haplotypes_kt.schema)

    # Annotate VDS with haplotypes (each variant now has a va.haps, which is a set of pairs (stored as a set of strings (ID2's)))
    vds = vds.annotate_variants_table(haplotypes_kt, root='va.haps', product=True)
    vds = vds.annotate_variants_expr(['va.haps = va.haps.toSet()',
                                      'va.hom_refs = gs.filter(g => g.isHomRef()).map(g => sa.meta.ID2).collect().toSet()',
                                      'va.hets = gs.filter(g => g.isHet()).map(g => sa.meta.ID2).collect().toSet()',
                                      'va.hom_alts = gs.filter(g => g.isHomVar()).map(g => sa.meta.ID2).collect().toSet()'
    ])
    # pprint(vds.variant_schema)

    # vds = vds.split_multi()  # Don't need to split since Andrea's VDS was already split

    vds = vds.annotate_variants_expr(
        ['va.shared_het = va.haps.filter(x => va.hets.contains(x.toArray()[0]) && va.hets.contains(x.toArray()[1])).size()',
         'va.shared_ref = va.haps.filter(x => va.hom_refs.contains(x.toArray()[0]) && va.hom_refs.contains(x.toArray()[1])).size()',
         'va.n_hom_ref = va.hom_refs.size()',
         'va.n_het = va.hets.size()',
         'va.n_hom_alt = va.hom_alts.size()'
         ])
    vds = vds.annotate_variants_expr(
        ['va.non_shared_het = va.n_het*(va.n_het - 1)//2 - va.shared_het',
         'va.non_shared_ref = va.n_hom_ref*(va.n_hom_ref - 1)//2 - va.shared_ref']
    )
    vds = process_consequences(vds)
    vds.export_variants(output_path, 'v.contig, v.start, v.ref, v.alt(), va.shared_het, va.shared_ref, va.n_hom_ref, va.n_het, va.n_hom_alt, va.non_shared_het, va.non_shared_ref')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--write_haplotypes', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    try_slack('@konradjk', main, args)