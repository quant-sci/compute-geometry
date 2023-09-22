import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

plt.rc('text', usetex=True)
plt.rc('font', size=12)

sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 2.5})

circle_color = "lightslategray"
points_color = "navy"

import time

from minimum_circle import ComputeMinimumCircle
sizes = [20, 200, 2000, 200000, 2000000]

import sys
sys.path.append("../")

if __name__ == "__main__":
    circle_computer = ComputeMinimumCircle()
    sizes = [10, 20, 30, 40, 50]  # Example list of input sizes
    circle_computer.plot_min_circle_random(sizes)