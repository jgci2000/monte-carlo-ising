import numpy as np
import matplotlib.pyplot as plt
import sys

jdos_exact = np.loadtxt("JDOS_exact_L4_SS.txt")

run_max = 100
f_exp_vals = np.arange(1, 9 + 1)
flatness = np.array([70, 80, 90, 95, 99])

error = np.zeros((len(flatness), len(f_exp_vals)))
abs_error = np.zeros((len(flatness), len(f_exp_vals)))
std_error = np.zeros((len(flatness), len(f_exp_vals)))
steps = np.zeros((len(flatness), len(f_exp_vals)))
std_steps = np.zeros((len(flatness), len(f_exp_vals)))
wall_time = np.zeros((len(flatness), len(f_exp_vals)))
std_wall_time = np.zeros((len(flatness), len(f_exp_vals)))

for i, p in enumerate(flatness):
    print(f"p = {i}/{len(flatness) - 1}", end='\r')
    for j, f_exp in enumerate(f_exp_vals):
        error_tmp = list()
        abs_error_tmp = list()
        steps_tmp = list()
        wall_time_tmp = list()
        
        for run in range(run_max):
            jdos = np.loadtxt("data/f_" + str(f_exp) + "/flatness_" + str(p) + "/" + str(run) + "_JDOS.txt")
            error_tmp.append(np.sum((jdos[jdos != 0] - jdos_exact[jdos != 0]) / jdos_exact[jdos != 0]))
            abs_error_tmp.append(np.sum(np.abs(jdos[jdos != 0] - jdos_exact[jdos != 0]) / jdos_exact[jdos != 0]))
            
            aux = np.loadtxt("data/f_" + str(f_exp) + "/flatness_" + str(p) + "/debug/" + str(run) + "_JDOS.txt", skiprows=1) 
            steps_tmp.append(np.sum(aux[:, 0]))
            wall_time_tmp.append(np.sum(aux[:, 1]))
        
        error[i, j] = np.mean(error_tmp)
        abs_error[i, j] = np.mean(abs_error_tmp)
        std_error[i, j] = np.std(error_tmp)
        steps[i, j] = np.mean(steps_tmp)
        std_steps[i, j] = np.std(steps_tmp)
        wall_time[i, j] = np.mean(wall_time_tmp)
        std_wall_time[i, j] = np.std(wall_time_tmp)

if sys.argv[1] == "save":
    np.savetxt("data/steps.csv", steps, delimiter=",")
    np.savetxt("data/abs_error.csv", abs_error, delimiter=",")
    np.savetxt("data/std_error.csv", std_error, delimiter=",")
    np.savetxt("data/wall_time.csv", wall_time, delimiter=",")
        
elif sys.argv[1] == "plot":
    plt.figure("Absolute Error WL vs f_final")
    plt.title(r"Absolute Error in WL vs $f_{final}$")
    plt.xlabel(r"$\log(f_{final})$")
    plt.ylabel(r"$|\epsilon|$")
    for i in range(len(flatness)):
        plt.plot(f_exp_vals, abs_error[i, :], 'o-', label=f"$p={flatness[i]}$")
    plt.legend()

    plt.figure("Error WL vs f_final")
    plt.title(r"Error in WL vs $f_{final}$")
    plt.xlabel(r"$\log(f_{final})$")
    plt.ylabel(r"$\epsilon$")
    for i in range(len(flatness)):
        plt.plot(f_exp_vals, error[i, :], 'o-', label=f"$p={flatness[i]}$")
    plt.legend()

    plt.figure("Absolute Error WL vs steps")
    plt.title(r"Absolute Error in WL vs steps")
    plt.xlabel(r"$\log(steps)$")
    plt.ylabel(r"$|\epsilon|$")
    for i in range(len(flatness)):
        plt.plot(np.log10(steps[i, :]), abs_error[i, :], 'o-', label=f"$p={flatness[i]}$")
    plt.legend()

    plt.figure("Wall Time WL vs steps")
    plt.title(r"Wall Time (s) in WL vs steps")
    plt.xlabel(r"$\log(steps)$")
    plt.ylabel(r"Wall Time (s)")
    for i in range(len(flatness)):
        plt.plot(np.log10(steps[i, :]), wall_time[i, :], 'o-', label=f"$p={flatness[i]}$")
    plt.legend()

    plt.show()


