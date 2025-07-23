import numpy as np
import matplotlib.pyplot as plt

# Set random seed for reproducibility
np.random.seed(42)

def generate_mean_std_pairs(low, high, n, num_groups):
	total_samples = n * num_groups
	data = np.random.randint(low, high, size=total_samples).reshape(-1, n)
	means = np.mean(data, axis=1)
	stds = np.std(data, axis=1, ddof=1)
	return means, stds

# Parameters
sample_sizes = [3, 5, 7, 10, 20, 50, 100]
value_ranges = [5, 10, 50, 100, 500]
num_groups = 100  # how many mean-std points per condition

# Create plot grid
fig, axs = plt.subplots(len(value_ranges), len(sample_sizes), figsize=(18, 10))
fig.suptitle("Mean vs Standard Deviation Across Sample Sizes and Value Ranges", fontsize=16)

for i, vr in enumerate(value_ranges):
	for j, n in enumerate(sample_sizes):
		means, stds = generate_mean_std_pairs(1, 1+vr, n, num_groups)

		ax = axs[i, j]
		ax.scatter(means, stds, s=10, edgecolors='none')
		ax.set_title(f"Range={vr}, n={n}", fontsize=10)
		ax.set_xlabel("Mean")
		ax.set_ylabel("Standard Deviation")

plt.tight_layout()
plt.savefig("check.png", dpi=100)
