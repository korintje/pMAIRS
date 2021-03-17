import pandas as pd
import numpy as np
import mairs

angles = [44.0, 39.0, 34.0, 29.0, 24.0, 19.0, 14.0, 9.0]
fnames = ["5S1", "5S2", "5S3", "5S4", "5S5", "5S6", "5S7", "5S8"]
fnames_bkg = ["5B1", "5B2", "5B3", "5B4", "5B5", "5B6", "5B7", "5B8"]

sample = mairs.MultiBeams()
for fname, angle in zip(fnames, angles):
    df = pd.read_csv("test/{}.CSV".format(fname), header=None)
    spectrum = np.array([df[0].values, df[1].values])
    sample.load_beam(angle, spectrum)

bkg = mairs.MultiBeams()
for fname, angle in zip(fnames_bkg, angles):
    df = pd.read_csv("test/{}.CSV".format(fname), header=None)
    spectrum = np.array([df[0].values, df[1].values])
    bkg.load_beam(angle, spectrum)

dataset = mairs.MultiBeamsSet()
dataset.load_background(bkg)
dataset.load_sample(sample)

result_df = dataset.get_op_ip()
result_df.to_csv('out.csv')

thetas_df = dataset.get_thetas()
thetas_df.dropna().to_csv('theta.csv')