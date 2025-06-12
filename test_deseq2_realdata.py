import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from preprocessing import preprocessing

if __name__ == "__main__":
    OUTPUT_PATH = "./output_files/"
    os.makedirs(OUTPUT_PATH, exist_ok=True) 

    counts_df, metadata = preprocessing()
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()

    conditions = metadata["condition"].value_counts().index[:2].tolist()
    print("Using contrast between:", conditions)

    ds = DeseqStats(dds, contrast=["condition", conditions[0], conditions[1]], inference=inference)
    ds.summary()

    ds.results_df.to_csv(os.path.join(OUTPUT_PATH, "results.csv"))