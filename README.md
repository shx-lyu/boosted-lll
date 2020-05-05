# boosted-lll
lattice reduction for communications. The project contains the implementations of LLL, boosted LLL, deep LLL and SR-SIC.

# How to use
run main_Fig1.m or main_Fig2.m

# on main_Fig1.m
Comparing the performance of LLL, boosted LLL, deep LLL and SR-SIC. The lattice bases are random correlated matrices.
The performance metrics are: {'Orthogonality defect','Basis length','Shortest Vector','Average reduction time/s'}.
Set the $M$ value to indicate the desired metric. By default, $M=1$ (Orthogonality defect).

# on main_Fig2.m
Comparing the performance of LLL, boosted LLL, deep LLL and SR-SIC. The lattice bases are taken from Integer-Forcing linear receivers.
The performance metrics are: {'Ergodic rate $R_E$/bpcu','Basis length','Shortest Vector','Average reduction time/s'}.
Set the $M$ value to indicate the desired metric. By default, $M=2$ (Basis length).

# Maintainer
The project is currently maintained by:
Shanxiang Lyu, shanxianglyu@gmail.com

# Contributing
boosted LLL welcomes contributions.
