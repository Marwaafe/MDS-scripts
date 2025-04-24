[mafechkar@ares MDS_Data]$ python compare_cellranger_versions_qc.py
Missing file: /trinity/home/mafechkar/MDS_OUTS_CellRangerCount/MDS001-09-203_count_output/MDS001-09-203_count/outs/metrics_summary.csv
Missing file: /trinity/home/mafechkiar/MDS_OUTS_CellRangerCount_9.0/MDS001-09-203_count_output/MDS001-09-203_count/outs/metrics_summary.csv
Traceback (most recent call last):
  File "/trinity/home/mafechkar/miniconda3/lib/python3.12/site-packages/pandas/core/indexes/base.py", line 3805, in get_loc
    return self._engine.get_loc(casted_key)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "index.pyx", line 167, in pandas._libs.index.IndexEngine.get_loc
  File "index.pyx", line 196, in pandas._libs.index.IndexEngine.get_loc
  File "pandas/_libs/hashtable_class_helper.pxi", line 7081, in pandas._libs.hashtable.PyObjectHashTable.get_item
  File "pandas/_libs/hashtable_class_helper.pxi", line 7089, in pandas._libs.hashtable.PyObjectHashTable.get_item
KeyError: 'Estimated Number of Cells'

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/net/beegfs/scratch/mafechkar/MDS_Data/compare_cellranger_versions_qc.py", line 45, in <module>
    value = df_metrics.at[df_metrics.index[0], original_label]
            ~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/trinity/home/mafechkar/miniconda3/lib/python3.12/site-packages/pandas/core/indexing.py", line 2575, in __getitem__
    return super().__getitem__(key)
           ^^^^^^^^^^^^^^^^^^^^^^^^
  File "/trinity/home/mafechkar/miniconda3/lib/python3.12/site-packages/pandas/core/indexing.py", line 2527, in __getitem__
    return self.obj._get_value(*key, takeable=self._takeable)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/trinity/home/mafechkar/miniconda3/lib/python3.12/site-packages/pandas/core/frame.py", line 4214, in _get_value
    series = self._get_item_cache(col)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/trinity/home/mafechkar/miniconda3/lib/python3.12/site-packages/pandas/core/frame.py", line 4638, in _get_item_cache
    loc = self.columns.get_loc(item)
          ^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/trinity/home/mafechkar/miniconda3/lib/python3.12/site-packages/pandas/core/indexes/base.py", line 3812, in get_loc
    raise KeyError(key) from err
KeyError: 'Estimated Number of Cells'
[mafechkar@ares MDS_Data]$



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

# Version-specific output directories
versions = {
    "7.2.0": "MDS_Data/MDS_OUTS_CellRangerCount",
    "9.0.0": "MDS_Data/MDS_OUTS_CellRangerCount_9.0.0"
}

# Metrics to extract
metrics_to_extract = {
    "Estimated Number of Cells": "Estimated_Cells",
    "Mean Reads per Cell": "Mean_Reads_Per_Cell",
    "Median Genes per Cell": "Median_Genes_Per_Cell",
    "Fraction Reads in Cells": "Fraction_Reads_In_Cells"
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
            for original_label, short_label in metrics_to_extract.items():
                value = df_metrics.at[df_metrics.index[0], original_label]
                row[short_label] = float(str(value).replace(",", ""))
            data.append(row)
        else:
            print(f"Missing file: {csv_path}")

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

# Boxplot for Mean Reads per Cell
plt.figure(figsize=(6, 5))
sns.boxplot(data=df, x="Version", y="Mean_Reads_Per_Cell", palette="pastel")
plt.title("Distribution of Mean Reads per Cell by Version")
plt.ylabel("Mean Reads per Cell")
plt.xlabel("Cell Ranger Version")
plt.tight_layout()
plt.show()
