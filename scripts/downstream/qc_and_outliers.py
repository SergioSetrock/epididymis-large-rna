import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import LocalOutlierFactor
import argparse
import os

def log2_cpm(counts):
    """Calculate log2(Counts Per Million + 1) for normalization."""
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
    return np.log2(cpm + 1)

def plot_pca(pca_df, variance, out_path, title, palette, markers):
    """Generate and save a PCA plot."""
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x='PC1', y='PC2', hue='Treatment', style='Generation',
        palette=palette, markers=markers, data=pca_df, s=150, edgecolor='k'
    )
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel(f"PC1 ({variance[0]*100:.1f}%)", fontweight='bold')
    plt.ylabel(f"PC2 ({variance[1]*100:.1f}%)", fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

def plot_correlation(log_cpm, out_path, title):
    """Generate and save a sample correlation heatmap."""
    corr_matrix = log_cpm.corr(method='pearson')
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_matrix, cmap='viridis', annot=False, xticklabels=True, yticklabels=True)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="QC, PCA, and KNN Outlier Detection.")
    parser.add_argument("--counts", required=True, help="Cleaned counts matrix")
    parser.add_argument("--meta", required=True, help="Metadata CSV")
    parser.add_argument("--outdir", required=True, help="Output directory for plots and tables")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Define colors and markers
    palette = {'Control': 'green', 'iAs': 'salmon'}
    markers = {'F1': 'o', 'F2': '^'}

    # 1. Load Data
    counts = pd.read_csv(args.counts, sep='\t', index_col=0)
    meta = pd.read_csv(args.meta, index_col='Sample')
    
    # Ensure metadata matches counts columns
    shared_samples = [s for s in counts.columns if s in meta.index]
    counts = counts[shared_samples]
    meta = meta.loc[shared_samples]

    # Filter lowly expressed genes for faster/cleaner PCA (sum > 10 across all samples)
    counts_filtered = counts[counts.sum(axis=1) > 10]

    # 2. Normalize (log2 CPM)
    log_cpm = log2_cpm(counts_filtered)

    # 3. Pre-QC Correlation and PCA
    plot_correlation(log_cpm, f"{args.outdir}/01_Pre_QC_Correlation.png", "Pre-QC Sample Correlation")
    
    # Run PCA
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(log_cpm.T)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(scaled_data)
    
    pca_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'], index=log_cpm.columns)
    pca_df = pca_df.join(meta)
    
    plot_pca(pca_df, pca.explained_variance_ratio_, f"{args.outdir}/02_Pre_QC_PCA.png", 
             "Pre-QC PCA", palette, markers)

    # 4. KNN Outlier Detection (Local Outlier Factor)
    # n_neighbors=5 is a good default for ~20 samples. 
    lof = LocalOutlierFactor(n_neighbors=5, contamination=0.10) # Assumes max ~10% outliers
    outlier_labels = lof.fit_predict(scaled_data) # 1 for inliers, -1 for outliers
    pca_df['Outlier_Status'] = ['Outlier' if x == -1 else 'Inlier' for x in outlier_labels]

    # Plot Outlier Identification
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x='PC1', y='PC2', hue='Outlier_Status', style='Generation',
        palette={'Inlier': 'blue', 'Outlier': 'red'}, markers=markers, data=pca_df, s=150, edgecolor='k'
    )
    plt.title("KNN Outlier Identification (LOF)", fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    for i, row in pca_df[pca_df['Outlier_Status'] == 'Outlier'].iterrows():
        plt.text(row['PC1']+1, row['PC2']+1, i, color='red', weight='bold') # Label the outliers
    plt.tight_layout()
    plt.savefig(f"{args.outdir}/03_KNN_Outlier_Identification.png", dpi=300)
    plt.close()

    # 5. Remove Outliers and Re-plot
    inlier_samples = pca_df[pca_df['Outlier_Status'] == 'Inlier'].index.tolist()
    counts_clean = counts[inlier_samples]
    meta_clean = meta.loc[inlier_samples]
    log_cpm_clean = log_cpm[inlier_samples]

    print(f"Removed {len(counts.columns) - len(inlier_samples)} outliers.")

    # Post-QC Correlation
    plot_correlation(log_cpm_clean, f"{args.outdir}/04_Post_QC_Correlation.png", "Post-QC Sample Correlation")

    # Post-QC PCA
    scaled_data_clean = scaler.fit_transform(log_cpm_clean.T)
    pca_clean = PCA(n_components=2)
    pcs_clean = pca_clean.fit_transform(scaled_data_clean)
    
    pca_df_clean = pd.DataFrame(pcs_clean, columns=['PC1', 'PC2'], index=log_cpm_clean.columns)
    pca_df_clean = pca_df_clean.join(meta_clean)
    
    plot_pca(pca_df_clean, pca_clean.explained_variance_ratio_, f"{args.outdir}/05_Post_QC_PCA.png", 
             "Post-QC PCA (Outliers Removed)", palette, markers)

    # 6. Save final cleaned counts and metadata for DESeq2
    counts_clean.to_csv(f"{args.outdir}/counts_no_outliers.tsv", sep='\t')
    meta_clean.to_csv(f"{args.outdir}/metadata_no_outliers.csv")

if __name__ == "__main__":
    main()
