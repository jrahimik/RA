#!/bin/bash
#SBATCH --job-name="RA_lasso_linear"
#SBATCH --mem=64GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=Javad@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --time=6-00:00:00







## This script is the last part of the STEMI data
Rscript LASSO_linear_regression_compare_auc_dist_matched_javad_v2.R
