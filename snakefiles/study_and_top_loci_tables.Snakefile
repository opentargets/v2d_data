'''
Contains Snakemake rules to make the study and top loci tables.

The rules for each table are interspersed with each other due to
inter-dependencies.

'''

import pandas as pd
import numpy as np
import json

assert(pd.__version__ >= '0.24')

#
# Top loci table --------------------------------------------------------------
#

rule get_gwas_cat_assoc:
    ''' Download GWAS catalog association file
    '''
    input:
        HTTPRemoteProvider().remote(
            'https://www.ebi.ac.uk/gwas/api/search/downloads/alternative')
        # FTPRemoteProvider().remote(
        #     'ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated.tsv')
    output:
        tmpdir + '/{version}/gwas-catalog-associations_ontology-annotated.tsv'
    shell:
        'cp {input} {output}'

rule get_variant_index:
    ''' Download variant index site list
    '''
    input:
        GSRemoteProvider().remote(
            config['var_index_sitelist'], keep_local=KEEP_LOCAL)
    output:
        tmpdir + '/variant-annotation.sitelist.tsv.gz'
    shell:
        'cp {input} {output}'

rule extract_gwascat_rsids_from_variant_index:
    ''' Makes set of GWAS Catalog rsids and chrom:pos strings. Then reads
        these from the variant index. Takes ~2 mins.
    '''
    input:
        gwascat = rules.get_gwas_cat_assoc.output,
        invar = rules.get_variant_index.output
    output:
        tmpdir + '/{version}/variant-annotation.sitelist.GWAScat.tsv.gz'
    shell:
        'pypy3 scripts/extract_from_variant-index.py '
        '--gwas {input.gwascat} '
        '--vcf {input.invar} '
        '--out {output}'

rule annotate_gwas_cat_with_variant_ids:
    ''' Annotates rows in the gwas catalog assoc file with variant IDs from
        a VCF
    '''
    input:
        gwascat = rules.get_gwas_cat_assoc.output,
        invar = rules.extract_gwascat_rsids_from_variant_index.output
    output:
        tmpdir + \
            '/{version}/gwas-catalog-associations_ontology_variantID-annotated.tsv'
    shell:
        'python scripts/annotate_gwascat_variantids.py '
        '--gwas {input.gwascat} '
        '--invar {input.invar} '
        '--out {output}'

# "Study table" rule that needs to be above `convert_gwas_catalog_to_standard`
rule make_gwas_cat_studies_table:
    ''' Make GWAS Catalog table.

        We need to join GWAS Catalog's study table with the top loci
        (association) table in order to get subphenotype in the
        "P-VALUE (TEXT)" field.

        We can't use the top loci table by itself as this would only include
        studies with >0 associations (may not be valid for phewas studies).
    '''
    input:
        toploci = rules.annotate_gwas_cat_with_variant_ids.output,
        gwas_study = HTTPRemoteProvider().remote(
            'https://www.ebi.ac.uk/gwas/api/search/downloads/studies_alternative', keep_local=KEEP_LOCAL),
        ancestries = HTTPRemoteProvider().remote(
            'https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry', keep_local=KEEP_LOCAL)
    output:
        main = tmpdir + '/{version}/gwas-catalog_study_table.json',
        lut = tmpdir + '/{version}/gwas-catalog_study_id_lut.tsv'
    shell:
        'python scripts/make_gwas_cat_study_table.py '
        '--in_gwascat_study {input.gwas_study} '
        '--in_toploci {input.toploci} '
        '--in_ancestries {input.ancestries} '
        '--outf {output.main} '
        '--out_id_lut {output.lut}'

rule convert_gwas_catalog_to_standard:
    ''' Outputs the GWAS Catalog association data into a standard format
    '''
    input:
        gwas = rules.annotate_gwas_cat_with_variant_ids.output,
        id_lut = rules.make_gwas_cat_studies_table.output.lut
    output:
        out_assoc = tmpdir + \
            '/{version}/gwas-catalog-associations_ot-format.tsv',
        log = 'logs/{version}/gwas-cat-assocs.log'
    shell:
        'python scripts/format_gwas_assoc.py '
        '--inf {input.gwas} '
        '--id_lut {input.id_lut} '
        '--outf {output.out_assoc} '
        '--log {output.log}'

rule cluster_gwas_catalog:
    ''' GWAS Catalog contains none independent loci from curation error result-
        ing in summary statistics being bulk imported from supplementary tables.
        To correct this, an additional clustering step has been added.

        In addition to clustering, this script also filters by p-value.
    '''
    input:
        rules.convert_gwas_catalog_to_standard.output.out_assoc
    output:
        out_assoc = tmpdir + \
            '/{version}/gwas-catalog-associations_ot-format.clustered.tsv',
        log = 'logs/{version}/gwas-cat-assocs_clustering.log'
    params:
        min_p = config['gwas_cat_min_pvalue'],
        dist = config['gwas_cat_cluster_dist_kb'],
        min_loci = config['gwas_cat_cluster_min_loci'],
        multi_prop = config['gwas_cat_cluster_multi_proportion']
    shell:
        'python scripts/cluster_gwas_catalog_associations.py '
        '--inf {input} '
        '--outf {output.out_assoc} '
        '--min_p {params.min_p} '
        '--cluster_dist_kb {params.dist} '
        '--cluster_min_loci {params.min_loci} '
        '--cluster_multi_prop {params.multi_prop} '
        '--log {output.log}'

# "Study table" rule that need to be above `merge_study_tables`
rule make_UKB_studies_table:
    ''' Makes study table for UKB summary statistics (Neale v2 and SAIGE)
    '''
    input:
        manifest = GSRemoteProvider().remote(
            config['ukb_manifest'], keep_local=KEEP_LOCAL),
        efos = GSRemoteProvider().remote(
            config['ukb_efo_curation'], keep_local=KEEP_LOCAL)
    output:
        study_table = tmpdir + '/{version}/UKB_study_table.json',
        prefix_counts = tmpdir + '/{version}/UKB_prefix_counts.tsv'
    shell:
        'python scripts/make_UKB_study_table.py '
        '--in_manifest {input.manifest} '
        '--in_efos {input.efos} '
        '--prefix_counts {output.prefix_counts} '
        '--outf {output.study_table}'

# "Study table" rule that need to be above `make_summarystat_toploci_table`
rule merge_study_tables:
    ''' Merges the GWAS Catalog and Neale UK Biobank study tables together.
    '''
    input:
        gwas = rules.make_gwas_cat_studies_table.output.main,
        ukb = rules.make_UKB_studies_table.output.study_table
    output:
        tmpdir + '/{version}/merged_study_table.old.json'
    shell:
        'python scripts/merge_study_tables.py '
        '--in_gwascat {input.gwas} '
        '--in_ukb {input.ukb} '
        '--output {output}'

rule merge_FINNGEN_study_tables:
    input:
        old_table=rules.merge_study_tables.output,
        finn_manifest=config['FINNGEN_manifest'],
        finn_efo=config['FINNGEN_efo_curation']
    output:
        tmpdir + '/{version}/merged_study_table.json'
    shell:
        'python scripts/Make_FINNGEN_entries.py '
        '--in_manifest {input.finn_manifest} '
        '--in_EFO {input.finn_efo} '
        '--in_study_table {input.old_table}'
        '--outf {output} '

rule make_summarystat_toploci_table:
    ''' Converts the toploci table produce from the finemapping pipeline to
        standardised format. Study table is need to know if a study is
        case-control or not.
    '''
    input:
        toploci = GSRemoteProvider().remote(
            config['toploci'], keep_local=KEEP_LOCAL),
        study_info = rules.merge_FINNGEN_study_tables.output
    output:
        tmpdir + '/{version}/sumstat-associations_ot-format.tsv'
    shell:
        'python scripts/format_sumstat_toploci_assoc.py '
        '--inf {input.toploci} '
        '--study_info {input.study_info} '
        '--outf {output}'

rule merge_gwascat_and_sumstat_toploci:
    ''' Merges associations from gwas_cat and sumstat derived
    '''
    input:
        gwascat = rules.cluster_gwas_catalog.output.out_assoc,
        sumstat = rules.make_summarystat_toploci_table.output
    output:
        'output/{version}/toploci.parquet'
    shell:
        'python scripts/merge_top_loci_tables.py '
        '--in_gwascat {input.gwascat} '
        '--in_sumstat {input.sumstat} '
        '--output {output}'

#
# Study table -----------------------------------------------------------------
#

rule get_efo_categories:
    ''' Uses OLS API to get "therapeutic area" / category for each EFO
    '''
    input:
        study_table = rules.merge_FINNGEN_study_tables.output,
        therapeutic_areas = config['efo_therapeutic_areas']
    output:
        tmpdir + '/{version}/efo_categories.json'
    shell:
        'python scripts/get_therapeutic_areas.py '
        '--in_study {input.study_table} '
        '--in_ta {input.therapeutic_areas} '
        '--output {output}'

rule list_studies_with_sumstats:
    ''' Makes a list of files with sumstats
    '''
    params:
        url=config['gcs_sumstat_pattern'],
        prev_url=config['prev_gcs_sumstat_pattern']
    output:
        tmpdir + '/{version}/studies_with_sumstats.tsv'
    run:
        shell('gsutil -m ls -d "{params.url}" > {output}')
	shell('gsutil -m ls -d "{params.prev_url}" >> {output}')

rule study_table_to_parquet:
    ''' Converts study table to final parquet.
        
        study_table: merged study table
        sumstat_studies: list of study IDs with sumstats
        top_loci: top loci table (needed to add number of associated loci)
        efo_anno: EFO to therapeautic area mapping json
        therapeutic_areas: list of "therapeutic areas", needed to sort the
                           efo_anno
        
        This task will fail if there are any study IDs with summary stats
        that are not found in the merged study table, this is required
        because of https://github.com/opentargets/genetics/issues/354 
    '''
    input:
        study_table = rules.merge_FINNGEN_study_tables.output,
        sumstat_studies = rules.list_studies_with_sumstats.output,
        efo_anno = rules.get_efo_categories.output,
        top_loci = rules.merge_gwascat_and_sumstat_toploci.output,
        therapeutic_areas = config['efo_therapeutic_areas']
    output:
        'output/{version}/studies.parquet'
    params:
        unknown_label = config['therapeutic_area_unknown_label']
    shell:
        'python scripts/study_table_to_parquet.py '
        '--in_study_table {input.study_table} '
        '--in_toploci {input.top_loci} '
        '--sumstat_studies {input.sumstat_studies} '
        '--efo_categories {input.efo_anno} '
        '--in_ta {input.therapeutic_areas} '
        '--unknown_label {params.unknown_label} '
        '--output {output}'
