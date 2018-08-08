#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import argparse
import pandas as pd
import numpy as np
import re

def main():

    # Parse args
    args = parse_args()

    #
    # Load GWAS Catalog data
    #

    # Load gwas catalog data
    gwas_raw = pd.read_csv(args.gwas, sep='\t', header=0, dtype={'P-VALUE':str})

    # Set a row (association) id
    gwas_raw['assoc_id'] = list(range(1, gwas_raw.shape[0] + 1))
    gwas = gwas_raw.copy()

    # Split rows by multiple variants
    gwas['CHR_ID'] = gwas['CHR_ID'].astype(str).str.split(';')
    gwas['CHR_POS'] = gwas['CHR_POS'].astype(str).str.split(';')
    gwas['SNPS'] = gwas['SNPS'].astype(str).str.split('; ')
    gwas['STRONGEST SNP-RISK ALLELE'] = gwas['STRONGEST SNP-RISK ALLELE'].astype(str).str.split('; ')

    # Assert that same number of variants exist
    assert(((gwas['CHR_ID'].apply(len)) == (gwas['CHR_POS'].apply(len))).all())
    assert(((gwas['SNPS'].apply(len)) == (gwas['STRONGEST SNP-RISK ALLELE'].apply(len))).all())
    # Filter rows (2 exist) whether CHR_ID != SNPS length
    to_keep = ((gwas['CHR_ID'].apply(len)) == (gwas['SNPS'].apply(len)))
    gwas = gwas.loc[to_keep, :]
    # Add column to flag that multiple variants were reported
    gwas['multi_SNPs_reported'] = np.where(gwas['CHR_ID'].apply(len) > 1, 1, 0)

    # Explode multivars into different rows
    gwas = explode(gwas, ['CHR_ID', 'CHR_POS', 'SNPS', 'STRONGEST SNP-RISK ALLELE'])

    #
    # Load HRC VCF data
    #

    # Load
    vcf_hrc = pd.read_csv(args.vcf_hrc, sep='\t', header=None).iloc[:, 0:5]
    vcf_hrc.columns = ['chrom', 'pos', 'rsid', 'ref', 'alt']

    # Assert only 1 rsid per row
    assert((vcf_hrc.rsid.str.split(',').apply(len) == 1).all())

    # Split multiallelic
    vcf_hrc['alt'] = vcf_hrc['alt'].str.split(',')
    vcf_hrc['is_multiallelic'] = np.where(vcf_hrc['alt'].apply(len) > 1, 1, 0)
    vcf_hrc = explode(vcf_hrc, ['alt'])

    #
    # Merge to HRC. First on rsids, then on chr:pos
    #

    # Make sure types match
    vcf_hrc['chrom'] = vcf_hrc['chrom'].astype(str)
    vcf_hrc['pos'] = vcf_hrc['pos'].astype(int)

    # Make new rsid column based on either SNP_ID_CURRENT or SNPS
    gwas['best_rsid'] = gwas.apply(get_best_rsid, axis=1)

    # Merge on rsids
    gwas_rs = gwas.loc[gwas.best_rsid.str.startswith('rs'), :]
    merged_rs = pd.merge(gwas_rs, vcf_hrc,
                         how='left',
                         left_on='best_rsid', right_on='rsid')

    # Merge on chr:pos
    gwas_chr = gwas.loc[gwas.best_rsid.str.startswith('chr'), :]
    gwas_chr[['chrom_37', 'pos_37']] = ( gwas_chr.best_rsid.apply(str_to_chrompos)
                                                           .apply(pd.Series) )
    merged_chr = pd.merge(gwas_chr, vcf_hrc,
                          how='left',
                          left_on=['chrom_37', 'pos_37'],
                          right_on=['chrom', 'pos'])

    # Combine
    merged = pd.concat([merged_rs, merged_chr])


    # Split these into those with HRC mapping and those without
    has_hrc = ~pd.isnull(merged.rsid)
    merged_w_hrc = merged.loc[has_hrc, :]
    merged_wo_hrc = merged.loc[~has_hrc, gwas.columns]

    #
    # Load Ensembl VCF data
    #

    # Load
    vcf_ensembl = pd.read_csv(args.vcf_ensembl, sep='\t', header=None).iloc[:, 0:5]
    vcf_ensembl.columns = ['chrom', 'pos', 'rsid', 'ref', 'alt']

    # Assert only 1 rsid per row
    assert((vcf_ensembl.rsid.str.split(',').apply(len) == 1).all())

    # Split multiallelic
    vcf_ensembl['alt'] = vcf_ensembl['alt'].str.split(',')
    vcf_ensembl['is_multiallelic'] = np.where(vcf_ensembl['alt'].apply(len) > 1, 1, 0)
    vcf_ensembl = explode(vcf_ensembl, ['alt'])

    #
    # Merge to Ensembl. First on rsids, then on chr:pos
    #

    # Make sure types match
    vcf_ensembl['chrom'] = vcf_ensembl['chrom'].astype(str)
    vcf_ensembl['pos'] = vcf_ensembl['pos'].astype(int)

    # Make new rsid column based on either SNP_ID_CURRENT or SNPS
    merged_wo_hrc['best_rsid'] = merged_wo_hrc.apply(get_best_rsid, axis=1)

    # Merge on rsids
    merged_wo_hrc_rs = merged_wo_hrc.loc[merged_wo_hrc.best_rsid.str.startswith('rs'), :]
    merged_rs = pd.merge(merged_wo_hrc_rs, vcf_ensembl,
                         how='left',
                         left_on='best_rsid', right_on='rsid')

    # Merge on chr:pos
    merged_wo_hrc_chr = merged_wo_hrc.loc[merged_wo_hrc.best_rsid.str.startswith('chr'), :]
    merged_wo_hrc_chr[['chrom_37', 'pos_37']] = ( merged_wo_hrc_chr.best_rsid.apply(str_to_chrompos)
                                                           .apply(pd.Series) )
    merged_chr = pd.merge(merged_wo_hrc_chr, vcf_ensembl,
                          how='left',
                          left_on=['chrom_37', 'pos_37'],
                          right_on=['chrom', 'pos'])

    # Combine
    ens_merged = pd.concat([merged_rs, merged_chr])

    #
    # Combine HRC and Ensembl merges
    #

    hrc_ens_merged = pd.concat([merged_w_hrc, ens_merged])


    # Drop rows without mapping in VCF
    hrc_ens_merged = hrc_ens_merged.dropna(subset=['rsid'])
    hrc_ens_merged['pos'] = hrc_ens_merged['pos'].astype(int)

    #
    # Create variant IDs
    #

    # Make a variant_id
    hrc_ens_merged['variant_id_b37'] = hrc_ens_merged.apply(make_var_id, axis=1)

    #
    # For rows with a risk allele reported, see if concordant with ref or alt
    # alleles on forward or reverse strands
    #

    # Extract risk allele
    hrc_ens_merged['risk_allele'] = hrc_ens_merged['STRONGEST SNP-RISK ALLELE'].apply(extract_risk_allele)
    # Check concordancy
    hrc_ens_merged['risk_allele_concordant'] = hrc_ens_merged[['risk_allele', 'ref', 'alt']].apply(check_concordancy, axis=1)
    # Remove rows where not concordant
    hrc_ens_merged = hrc_ens_merged.loc[hrc_ens_merged['risk_allele_concordant'] != 'no', :]

    #
    # Tidy
    #

    # Drop unneeded columns
    hrc_ens_merged = hrc_ens_merged.drop(labels=[
        'alt',
        'best_rsid',
        'chrom',
        'chrom_37',
        'pos',
        'pos_37',
        'risk_allele',
        'risk_allele_concordant',
        'ref'], axis=1)
    hrc_ens_merged = hrc_ens_merged.drop_duplicates()


    #
    # Collapse multiple variants into one row per locus
    #

    # Collapse multiple variants into a single row separated by ';'
    multivar = hrc_ens_merged.loc[:, ['assoc_id', 'rsid', 'variant_id_b37']]
    multivar_collapsed = (
        multivar.groupby('assoc_id')
                .agg({'rsid':combine_rows, 'variant_id_b37':combine_rows})
                .reset_index() )

    # Merge collapsed variant IDs back to the raw gwas data
    out_df = pd.merge(gwas_raw, multivar_collapsed,
                      how='left',
                      on='assoc_id')

    out_df.to_csv(args.out, sep='\t', index=None)

    return 0

def check_concordancy(row):
    ''' Check whether risk allele is concordant with ref or alt on wither forward or reverse.
    Args:
        row (pd.series)
    Returns:
        str of yes, no or ambiguous
    '''
    if pd.isnull(row['risk_allele']):
        ret = 'ambiguous'
    elif any( [row['risk_allele'] == row['alt'],
               row['risk_allele'] == row['ref'],
               row['risk_allele'] == revcomp(row['alt']),
               row['risk_allele'] == revcomp(row['ref'])] ):
        ret = 'yes'
    else:
        ret = 'no'
    return ret

def comp(s):
    """ Complementary sequence """
    com_dict = {"A":"T", "T":"A", "G":"C", "C":"G"}
    com_seq = "".join([com_dict.get(nucl, "N") for nucl in s])
    return com_seq

def revcomp(s):
    """ Reverse complement sequence """
    return comp(s[::-1])

def extract_risk_allele(s):
    ''' Takes a string from GWAS Catalog STRONGEST SNP-RISK ALLELE field and
        returns the risk allele, if any.
    Args:
        row (pd.series)
    Returns:
        str
    '''
    mtch = re.match(r'.*-([A|T|G|C]+)', s)
    if mtch:
        allele = mtch.group(1)
    else:
        allele = None
    return allele

def combine_rows(items):
    return ';'.join(items)

def make_var_id(row):
    ''' Make variant ID
    Args:
        row (pd.series)
    Returns:
        str
    '''
    # print(row)
    return '{chrom}_{pos}_{ref}_{alt}'.format(**row)

def get_best_rsid(row):
    ''' Returns the best rsid from:
            1. SNP_ID_CURRENT
            2. SNPS
            3. None
    Args:
        row (pd.series)
    Returns:
        str
    '''
    if not pd.isnull(row['SNP_ID_CURRENT']):
        try:
            return 'rs{0}'.format(int(row['SNP_ID_CURRENT']))
        except ValueError:
            return 'rs{0}'.format(int(row['SNP_ID_CURRENT'].replace('b', '')))
    elif not pd.isnull(row['SNPS']):
        return row['SNPS']
    return None

def str_to_chrompos(s):
    ''' Converts gwas catalog snp str to chrom and pos
    Args:
        s (str)
    Returns:
        str, int
    '''

    chrom, pos = ( s.replace('chr', '')
                    .replace('_', ':')
                    .replace('.', ':')
                    .replace('-', ':')
                    .split(':')[:2] )
    return str(chrom), int(pos)

def explode(df, columns):
    ''' Explodes multiple columns
    '''
    idx = np.repeat(df.index, df[columns[0]].str.len())
    a = df.T.reindex_axis(columns).values
    concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
    p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
    return pd.concat([df.drop(columns, axis=1), p], axis=1).reset_index(drop=True)

def parse_args():
    """ Load command line args """
    parser = argparse.ArgumentParser()
    parser.add_argument('--gwas', metavar="<file>", help=('GWAS Catalog input'), type=str, required=True)
    parser.add_argument('--vcf_hrc', metavar="<file>", help=("HRC VCF input"), type=str, required=True)
    parser.add_argument('--vcf_ensembl', metavar="<file>", help=("Ensembl VCF input"), type=str, required=True)
    parser.add_argument('--out', metavar="<str>", help=("Output"), type=str, required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
