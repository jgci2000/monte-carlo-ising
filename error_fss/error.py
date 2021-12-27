import numpy as np
import matplotlib.pyplot as plt

jdos_exact = np.loadtxt("JDOS_exact_L4_SS.txt")

run_max = 1000
rep_exp_vals = np.arange(2, 7 + 1)
skip_exp_vals = np.array([-2, -1, 0, 1, 2])

error = np.zeros((len(skip_exp_vals), len(rep_exp_vals)))
abs_error = np.zeros((len(skip_exp_vals), len(rep_exp_vals)))

for i, rep_exp in enumerate(rep_exp_vals):    
    print(f"rep = {i}/{len(rep_exp_vals) - 1}", end='\r')
    for j, skip_exp in enumerate(skip_exp_vals):
        error_tmp = list()
        abs_error_tmp = list()
        
        for run in range(run_max):
            jdos = np.loadtxt("data/REP_" + str(rep_exp) + "/skip_" + str(skip_exp) + "/" + str(run + 1) + "_JDOS.txt")
            error_tmp.append(np.sum((jdos[jdos != 0] - jdos_exact[jdos != 0]) / jdos_exact[jdos != 0]))
            abs_error_tmp.append(np.sum(np.abs(jdos[jdos != 0] - jdos_exact[jdos != 0]) / jdos_exact[jdos != 0]))

        error[j, i] = np.mean(error_tmp)
        abs_error[j, i] = np.mean(abs_error_tmp)

plt.figure("Absolute Error FSS vs f_final")
plt.title(r"Absolute Error in FSS vs $f_{final}$")
plt.xlabel(r"$\log(f_{final})$")
plt.ylabel(r"$|\epsilon|$")
plt.ylim([0, np.max(abs_error) * 1.1])
for i in range(len(skip_exp_vals)):
    plt.plot(rep_exp_vals, abs_error[i, :], 'o-', label=f"$skip={skip_exp_vals[i]}$")
plt.legend()

plt.figure("Error FSS vs f_final")
plt.title(r"Error in FSS vs $f_{final}$")
plt.xlabel(r"$\log(f_{final})$")
plt.ylabel(r"$\epsilon$")
for i in range(len(skip_exp_vals)):
    plt.plot(rep_exp_vals, error[i, :], 'o-', label=f"$skip={skip_exp_vals[i]}$")
plt.legend()

plt.show()
