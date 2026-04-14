import pandas as pd, matplotlib.pyplot as plt, pickle as pkl, numpy as np, glob, os
import argparse

def parse():
    return

if __name__=="__main__":
    files=glob.glob("*/probe*csv")
    for i in range(len(files)):
        if i==0:
            z=pd.read_csv(files[i],sep=" ")
        else:
            z=pd.concat([z,pd.read_csv(files[i],sep=" ")])
    z.to_csv("benchmark.csv",sep=" ",index=False,header=True)