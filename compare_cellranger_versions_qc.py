

[ df.head():
          Sample Version  Mean_Reads_Per_Cell  Median_Genes_Per_Cell  Fraction_Reads_In_Cells  Estimated_Cells
0  MDS001-09-203   7.2.0              16426.0                 1110.0                     96.6          15156.0
1  MDS001-09-203   9.0.0              16574.0                 1118.0                     96.3          15020.0
2  MDS005-09-247   7.2.0              77711.0                 1824.0                     69.4           3486.0
3  MDS005-09-247   9.0.0              80386.0                 1842.0                     68.6           3370.0
4  MDS006-08-249   7.2.0              52182.0                 1418.0                     75.5           3649.0

 df.columns:
Index(['Sample', 'Version', 'Mean_Reads_Per_Cell', 'Median_Genes_Per_Cell',
       'Fraction_Reads_In_Cells', 'Estimated_Cells'],
      dtype='object')

 Paired t-test result on Mean Reads per Cell:
T-statistic: -2.987 p-value: 0.0083
Significant difference between versions (p < 0.05)

 df_melted.head():
          Sample Version               Metric    Value
0  MDS001-09-203   7.2.0  Mean_Reads_Per_Cell  16426.0
1  MDS001-09-203   9.0.0  Mean_Reads_Per_Cell  16574.0
2  MDS005-09-247   7.2.0  Mean_Reads_Per_Cell  77711.0
3  MDS005-09-247   9.0.0  Mean_Reads_Per_Cell  80386.0
4  MDS006-08-249   7.2.0  Mean_Reads_Per_Cell  52182.0

 df_melted.columns:
Index(['Sample', 'Version', 'Metric', 'Value'], dtype='object')
/net/beegfs/scratch/mafechkar/MDS_Data/compare_cellranger_versions_qc.py:135: FutureWarning:

Passing `palette` without assigning `hue` is deprecated and will be removed in v0.14.0. Assign the `x` variable to `hue` and set `legend=False` for the same effect.





import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Real sample names
samples = [
    "MDS001-09-203", "MDS005-09-247", "MDS006-08-249", "MDS010-09-299",
    "MDS016-09-478", "MDS023-10-053", "MDS029-10-118", "MDS038-10-241",
    "MDS059-10-531", "MDS065-10-609", "MDS154-13-486", "MDS155-13-606",
    "MDS167-13-913", "MDS169-13-919", "MDS180-14-164", "MDS189-14-527",
    "MDS201-15-093", "MDS212-15-463"
]

# Output directories for different Cell Ranger versions
versions = {
    "7.2.0": "MDS_Data/MDS_OUTS_CellRangerCount",
    "9.0.0": "MDS_Data/MDS_OUTS_CellRangerCount_9.0.0"
}

# Standard metrics to extract (exact match)
metrics_to_extract = {
    "Mean Reads per Cell": "Mean_Reads_Per_Cell",
    "Median Genes per Cell": "Median_Genes_Per_Cell",
    "Fraction Reads in Cells": "Fraction_Reads_In_Cells"
}

# Metrics with fallback names
fallback_metrics = {
    "Estimated_Cells": ["Estimated Number of Cells", "Total Genes Detected"]
}

# Gather data
data = []

for sample in samples:
    for version, base_dir in versions.items():
        csv_path = os.path.join(
            base_dir,
            f"{sample}_count_output",
            f"{sample}_count",
            "outs",
            "metrics_summary.csv",
        )
        if os.path.exists(csv_path):
            df_metrics = pd.read_csv(csv_path, index_col=0).T
            row = {"Sample": sample, "Version": version}

            # Standard metrics
            for original_label, short_label in metrics_to_extract.items():
                if original_label in df_metrics.columns:
                    value = df_metrics.at[df_metrics.index[0], original_label]
                    row[short_label] = float(str(value).replace(",", ""))
                else:
                    print(f"⚠️  Missing metric '{original_label}' in file: {csv_path}")
                    row[short_label] = None

            # Fallback metrics
            for short_label, options in fallback_metrics.items():
                found = False
                for label in options:
                    if label in df_metrics.columns:
                        value = df_metrics.at[df_metrics.index[0], label]
                        row[short_label] = float(str(value).replace(",", ""))
                        found = True
                        break
                if not found:
                    print(f"⚠️  Missing fallback metrics {options} for '{short_label}' in file: {csv_path}")
                    row[short_label] = None

            data.append(row)
        else:
            print(f"❌ Missing file: {csv_path}")

# Build dataframe
df = pd.DataFrame(data)

# Export to CSV
df.to_csv("cellranger_comparison_metrics.csv", index=False)

# Bar plots for all metrics
df_melted = df.melt(id_vars=["Sample", "Version"], var_name="Metric", value_name="Value")
sns.set(style="whitegrid")
g = sns.catplot(
    data=df_melted,
    kind="bar",
    x="Sample",
    y="Value",
    hue="Version",
    col="Metric",
    col_wrap=2,
    sharey=False,
    height=4,
    aspect=1.5
)
g.set_titles("{col_name}")
g.fig.subplots_adjust(top=0.9)
g.fig.suptitle("Cell Ranger v7.2.0 vs v9.0.0 - Per Sample Bar Plots", fontsize=16)
plt.xticks(rotation=90)

#  Boxplot for Mean Reads per Cell
plt.figure(figsize=(6, 5))
sns.boxplot(data=df, x="Version", y="Mean_Reads_Per_Cell", palette="pastel")
plt.title("Distribution of Mean Reads per Cell by Version")
plt.ylabel("Mean Reads per Cell")
plt.xlabel("Cell Ranger Version")
plt.tight_layout()
plt.show()
