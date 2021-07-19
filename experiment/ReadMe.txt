This is the code of all experiments. For the Experiment 2-4, data are not included.

The subfolders represent different experiments.

Files description:

Experiment 1:
PGmodel_fig.m  --------- A script to generate the clipped Poissonian-Gaussian models.
dis1_i.mat     --------- The probability distribution of model i.
ConvexOpt.m    --------- Our algorithm. For the input dis1_i, it will stabilize the variance and save in result1_i.mat.
result1_i.mat  --------- The stabilized-variance and VST function.
exp1.m         --------- The main script of Experiment 1. It includes the implementation of peer methods and the comparison of all methods.

Experiment 2:
dataSimulation.m ------- A script to generate the synthetic data.
HistogramCount.m ------- The distribution estimation method with truncated Gaussian model fitting.
exp2.m         --------- The main script of Experiment 2.

Experiment 3:
exp3.m         --------- The main script of Experiment 3.

Experiment 4:
exp4.m         --------- The main script of Experiment 4.
PSNR.mat       --------- The PNSR and MSE vary of FOV, method or subsets.
exp4_plot.m    --------- Plot the box plot of PNSR based on PSNR.mat
bm3d           --------- BM3D algorithm


exp1.m can be run directly. For other experiments, please install the MOSEK optimizer and get the full data by this link:

https://drive.google.com/file/d/1P_wG8j2i7mBtF64GJytK-F_modcsxbaf/view?usp=sharing


