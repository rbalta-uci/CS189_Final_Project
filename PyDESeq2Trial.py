from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
import os

def load_data():
    files = os.listdir("data")
    
    for f in files:
        print(f)
        f = load_example_data(modality="raw_counts",dataset=f,debug=False,)
        print(f.head())
load_data()