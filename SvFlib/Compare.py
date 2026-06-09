import COMMON as co
from  os  import listdir, getcwd
import re
import numpy as np
import matplotlib.pyplot as plt

def find_res_files():
    """
    Find all .res files in the current directory

    Returns:
        list[str]: names of all .res files
    """
    res_files = []
    files = listdir(getcwd())
    for f in files:
        p  = f.rfind('.res')
        if p + 4 == len(f):
                res_files.append(f)
    return res_files


def get_points(filename):
    """
    Extract information about results of optimization from .res file

    Args:
        filename (str) - full name (with .res) of results file

    Returns:
        ndarray[int] - array of number of iterations of optimization
        ndarray[float] - array of results on each iteration of optimization
        ndarray[list[float]] - array of points evaluated in optimization
    """

    with open(filename) as f:
        text = f.read()
    
    matches = re.findall(r'Num\s+(\d+)\s+Val\s+([\d\.]+)\s+Arg\s+\[([\d\.]+)\]', text)
    points = sorted([(int(n), float(v), float(a)) for n, v, a in matches], key=lambda x: x[0])
    
    iterations = np.array([p[0] for p in points])
    vals = np.array([p[1] for p in points])
    args = np.array([p[2] for p in points])
    
    return iterations, vals, args


def compare_results():
    if co.ResToCompare.lower() == "All".lower():
        target_res_files = find_res_files()
    else:
        target_res_files = co.ResToCompare

    print(f"Target res files = {target_res_files}")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
   # plt.figure(figsize=(10,6))
    for res_file in target_res_files:
        iterations, vals, args = get_points(res_file)
        print(f"In {res_file:}")
        print(f"Vals: {vals}")
        print(f"Args: {args}")

        # Plot best-so-far curve starting after initial design
        history = vals

        # Best so far across all evaluations
        best_so_far = np.minimum.accumulate(history)
        # Start the red line after initial design
        x_best = iterations
        y_best = best_so_far
        ax1.plot(
            x_best,
            y_best,
            linewidth=2,
            label=res_file,
        )
        ax2.plot(
            x_best[12:],
            y_best[12:],
            linewidth=2,
            label=res_file,
        )

    ax1.set_xlabel("Iteration", fontsize=11)
    ax1.set_ylabel("CVError", fontsize=11)

    title = "Optimization Progress"
    ax1.set_title(title, fontsize=12)

    ax1.legend(fontsize=10)
    # ax1.yscale("log")
    ax1.grid(True, alpha=0.3)
    plt.savefig("Comparison_of_optimization_strategies.png")

    #plt.show()
