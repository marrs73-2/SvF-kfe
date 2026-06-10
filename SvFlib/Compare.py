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
        int - number of initial points
    """

    with open(filename) as f:
        lines = f.readlines()

    iterations = []
    vals = []
    args = []
    initial_count = 0
    
    # Read .res file line by line
    for line in lines:
            # Find lines with optimization point data in them using regular expression
            match = re.search(r'Num\s+(\d+)\s+Val\s+([\d\.]+)\s+Arg\s+\[([^\]]+)\]', line)

            # Extract data for found points
            if match:
                iterations.append(int(match.group(1)))
                vals.append(float(match.group(2)))
                args.append(match.group(3))
                # Проверяем наличие INITIAL в той же строке
                if 'INITIAL' in line: initial_count += 1
    
    # Sort points by iteration number (they're sometimes mixed in process)
    sorted_indices = np.argsort(iterations)

    return (np.array(iterations)[sorted_indices],
            np.array(vals)[sorted_indices],
            np.array(args)[sorted_indices],
            initial_count)



def compare_results():
    """
    Draw graphs of optimization progress of different tasks solved beforehand for comparison.
    
    Initial points are scattered in the left area and optimization progress 
    is shown as the line of best result over the iterations for each task.
    """
        
    # Get names of processed .res files
    if co.ResFilesToCompare.lower() == "All".lower():
        # If All option was used find all .res files in the directory
        target_res_files = find_res_files()
    else:
        # Otherways, take the inputted array of names
        target_res_files = co.ResFilesToCompare.strip().strip('[]').split(',')
        target_res_files = [file.strip().strip('"\'') for file in target_res_files]

    print(f"Target res files = {target_res_files}")
    print(f"co.Tail_start = {co.Tail_start}")
    tail_start = int(co.Tail_start)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    did_plot_initial = False
   
    # Draw data from each .res file
    for res_file in target_res_files:
        iterations, vals, args, initial_count = get_points(res_file)
        print(f"In {res_file:}")
        print(f"Vals: {vals}")
        print(f"Args: {args}")
        print(f"Iterations: {iterations}")
        print(f"Initial count: {initial_count}")

        if did_plot_initial == False:
            did_plot_initial = True
            ax1.axvspan(0, initial_count, alpha=0.15, color="gray", zorder=0)

            ax1.scatter(
                iterations[0:initial_count],
                vals[0:initial_count],
                alpha=0.6,
                s=50,
                label=f"Initial design (n={initial_count})",
                color="gray",
                edgecolors="black",
                linewidth=0.5,
                zorder=2,
            )
        # Plot best-so-far curve starting after initial design
        history = vals

        # Best so far across all evaluations
        best_so_far = np.minimum.accumulate(history)
        # Start the red line after initial design
        x_best = iterations[initial_count:]
        y_best = best_so_far[initial_count:]
        ax1.plot(
            x_best,
            y_best,
            linewidth=2,
            label=res_file,
        )

        if (tail_start > initial_count):
            ax2.plot(
                x_best[tail_start - initial_count:],
                y_best[tail_start - initial_count:],
                linewidth=2,
                label=res_file,
            )

    ax1.set_xlabel("Iteration", fontsize=11)
    ax1.set_ylabel("CVError", fontsize=11)
    ax2.set_xlabel("Iteration", fontsize=11)
    ax2.set_ylabel("CVError", fontsize=11)

    title1 = "Optimization Progress"
    ax1.set_title(title1, fontsize=12)

    title2 = "Tail of optimization Progress"
    ax2.set_title(title2, fontsize=12)

    ax1.legend(fontsize=10)
    ax2.legend(fontsize=10)
    # ax1.yscale("log")
    ax1.grid(True, alpha=0.3)
    ax2.grid(True, alpha=0.3)
    plt.savefig("Comparison_of_optimization_strategies.png")

    plt.tight_layout()
    plt.show()
