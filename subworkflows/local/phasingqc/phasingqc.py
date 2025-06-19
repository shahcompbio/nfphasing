import cyvcf2

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import click

vcf_filename = "phased_variants.vcf.gz"
patient_id = "patient_27976"

def read_vcf(vcf_path:str) -> pd.DataFrame:
    """
    Generator that yields DataFrame chunks of size `chunk_size`.
    """
    vcf = cyvcf2.VCF(vcf_path)
    records = []
    
    for variant in vcf:

        ps = None
        if variant.format('PS') is not None:
            assert variant.format('PS').shape == (1, 1)
            ps = variant.format('PS')[0, 0]

        gt_qual = None
        if variant.gt_quals.shape is not None:
            assert variant.gt_quals.shape == (1,)
            gt_qual = variant.gt_quals[0]

        assert len(variant.genotypes) == 1
        assert len(variant.genotypes[0]) == 3
        allele_a = variant.genotypes[0][0]
        allele_b = variant.genotypes[0][1]
        is_phased = variant.genotypes[0][2]

        records.append({
            'CHROM': variant.CHROM,
            'start': variant.start,
            'end': variant.end,
            'is_indel': variant.is_indel,
            'QUAL': variant.QUAL,
            'FILTER': variant.FILTER,
            'num_called': variant.num_called,
            'num_het': variant.num_het,
            'num_hom_alt': variant.num_hom_alt,
            'num_hom_ref': variant.num_hom_ref,
            'ps': ps,
            'gt_qual': gt_qual,
            'allele_a': allele_a,
            'allele_b': allele_b,
            'is_phased': is_phased,
        })

    vcf_data = pd.DataFrame(records)

    vcf_data['FILTER'] = vcf_data['FILTER'].fillna('PASS')
    vcf_data['ps'] = vcf_data['ps'].fillna(-1).astype(int)

    return vcf_data
    

def prepare_vcf_data(vcf_path: str) -> pd.DataFrame:
    """
    Prepare VCF data for analysis.
    """

    vcf_data:pd.DataFrame = read_vcf(vcf_path)
    # Filtering for het positions
    vcf_data = vcf_data.loc[(vcf_data['num_called'] == 1) &(vcf_data['num_het'] == 1),] #type:ignore
    return vcf_data


# Visualize genotype quality for het positions
def plot_genotype_quality(vcf_data: pd.DataFrame, patient_id: str):
    """
    Visualize genotype quality for het positions.
    """
    plt.figure(figsize=(3, 3))
    sns.histplot(data=vcf_data, x='gt_qual', bins=10)
    sns.despine()
    plt.title(patient_id)
    plt.savefig(f"{patient_id}_genotype_quality.png", dpi=300, bbox_inches='tight')


# Block stats
def calculate_block_stats(vcf_data: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate statistics for each phase block.
    """
  
    block_stats = vcf_data.query('ps >= 0').groupby('ps').agg(
        n_variants=('ps', 'size'),
        start=('start', 'min'),
        end=('end', 'max'),
        fraction_allele_a_1=('allele_a', 'mean'),
        fraction_indels=('is_indel', 'mean'),
    ).reset_index()
    block_stats['length'] = block_stats['end'] - block_stats['start'] + 1
    return block_stats


# Visualize number of variants per block
def plot_variants_per_block(block_stats: pd.DataFrame, patient_id: str):
    """
    Visualize number of variants per block.
    """
    plt.figure(figsize=(3, 3))
    sns.histplot(x='n_variants', data=block_stats, log_scale=True, bins=20)
    sns.despine()
    plt.title(patient_id)
    plt.savefig(f"{patient_id}_variants_per_block.png", dpi=300, bbox_inches='tight')
    plt.close()

# Visualize distribution of block lengths
def plot_block_length_distribution(block_stats: pd.DataFrame, patient_id: str):
    """
    Visualize distribution of block lengths.
    """
    plt.figure(figsize=(3, 3))
    sns.histplot(x='length', data=block_stats, log_scale=True, bins=20)
    sns.despine()
    plt.title(patient_id)
    plt.savefig(f"{patient_id}_block_length_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()

# Visualize fraction of allele a that are the alternate allele
def plot_fraction_allele_a(block_stats: pd.DataFrame, patient_id: str):
    """
    Visualize fraction of allele a that are the alternate allele.
    """
    plt.figure(figsize=(3, 3))
    sns.histplot(x='fraction_allele_a_1', data=block_stats.query('n_variants > 10'), bins=20)
    sns.despine()
    plt.title(patient_id)
    plt.savefig(f"{patient_id}_fraction_allele_a.png", dpi=300, bbox_inches='tight')
    plt.close()
# Compute n50
# Phasing N50 is the minimum phase block length, where the sum of its phase blocks with all larger phase blocks spans ≥50% of the total phase length.

# Code stolen from https://github.com/sinamajidian/phaseme
def compute_n50(block_stats: pd.DataFrame) -> int:
    """
    Compute the N50 of the phase blocks.
    N50 is the minimum phase block length, where the sum of its phase blocks with all larger phase blocks spans ≥50% of the total phase length.
    """    
    # Blocks sorted by length
    values_sorted = sorted(block_stats['length'].values, reverse=True)

    # Compute total length of blocks / 2
    n2 = int(sum(values_sorted)/2)

    # Compute the minimum cumulative value that is greater than equal to half
    # the total lengths of blocks, find its index
    csum = np.cumsum(values_sorted)
    csumn2 = min(csum[csum >= n2])
    ind = np.where(csum == csumn2)

    # Compute the n50: the length of the largest 
    n50 = values_sorted[int(ind[0])]

    return n50

def export_n50(block_stats: pd.DataFrame, patient_id: str):
    """
    Export N50 to a file.
    """
    try:
        n50 = compute_n50(block_stats)
    except ValueError as e:
        print(f"Error computing N50: {e}")
        n50 = -1
    with open(f"{patient_id}_n50.txt", 'w') as f:
        f.write(f"{n50}\n")


@click.command()
@click.option(
    "--vcf",
    type=click.Path(exists=True),
    required=True,
    help="Path to the VCF file to analyze.",
)
@click.option(
    "--sample",
    required=True,
    help="Sample ID for the analysis.",
)
def main(vcf, sample):
    """
    Main function to run the phasing QC analysis.
    """
    # Read and prepare VCF data
    vcf_data = prepare_vcf_data(vcf)

    # Plot genotype quality
    plot_genotype_quality(vcf_data, sample)

    # Calculate block stats
    block_stats = calculate_block_stats(vcf_data)

    # Plot variants per block
    plot_variants_per_block(block_stats, sample)

    # Plot block length distribution
    plot_block_length_distribution(block_stats, sample)

    # Plot fraction of allele a
    plot_fraction_allele_a(block_stats, sample)

    # Export N50
    export_n50(block_stats, sample)

if __name__ == "__main__":
    main()
