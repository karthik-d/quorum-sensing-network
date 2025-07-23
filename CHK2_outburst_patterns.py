import numpy as np
import matplotlib.pyplot as plt

np.random.seed(0)

# Define sample sizes and value ranges
sample_sizes = [3, 5, 10, 20]
value_ranges = [(10, 12), (10, 15), (10, 20), (10, 30)]
num_groups = 1000  # number of samples per setting

# Create subplot grid
fig, axs = plt.subplots(len(value_ranges), len(sample_sizes), figsize=(16, 12), sharex=True, sharey=True)
fig.suptitle("Discretization of Mean vs. Standard Deviation\nEffect of Sample Size and Value Range", fontsize=16)

for i, (low, high) in enumerate(value_ranges):
    for j, n in enumerate(sample_sizes):
        means = []
        stds = []

        for _ in range(num_groups):
            group = np.random.choice(np.arange(low, high + 1), size=n)
            means.append(np.mean(group))
            stds.append(np.std(group, ddof=1))

        ax = axs[i, j]
        ax.scatter(means, stds, s=10, alpha=0.6, color='firebrick')
        ax.set_title(f"Range=({low},{high}), n={n}", fontsize=10)
        ax.grid(True)

        if i == len(value_ranges) - 1:
            ax.set_xlabel("Mean")
        if j == 0:
            ax.set_ylabel("Standard Deviation")

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("check_radial.png")
