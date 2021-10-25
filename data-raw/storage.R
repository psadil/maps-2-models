
n <- 40000

# T1 (2 TB)
t1p <- 50 * n

# fMRI (2.68e7; 26.8 TB)
fp <- 670 * n

# total (TB: 28.8)
storage <- t1p / 1e6 + fp / 1e6

# price
# https://www.oracle.com/cloud/price-list.html#storage
sp <- 30 * 1000 * 0.0255 * 12

# compute E2 => free (https://www.oracle.com/cloud/price-list.html#compute-vm)

# GPU (https://www.oracle.com/cloud/price-list.html#compute-gpu)
# 300 hrs / month
gpu <- 300 * 12 * 2.95 * 2

gpu + sp
