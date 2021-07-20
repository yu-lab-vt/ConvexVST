# Welcome to ConvexVST
ConvexVST (A **Convex** Optimization Approach to **V**ariance-**S**tabilizing **T**ransformation) is an approach to solve the variance-stabilizing transformation (VST) problem, transforming heteroscedastic data to homoscedastic data so that they are more tractable for subsequent analysis. ConvexVST can cast the VST problem into a convex optimization problem, which can always be efficiently solved, identified the specific structure of the convex problem, which further improved the efficiency of the proposed algorithm, and showed that any finite discrete distributions and the discretized version of any continuous distributions from real data can be variance-stabilized in an easy and nonparametric way. 

If you have any feedback or issue, you are welcome to either post them or send email to yug@vt.edu or mengfanw@vt.edu (Guoqiang Yu and Mengfan Wang at Virginia Tech).

## Overview
### Variance-Stabilizing Transformation
#### Definition
VST problem is to transform heteroscedastic distributions to homoscedastic. That is, we want to find a transform function $f$ so that $Var[f(X_/theta)]$ is constant.


