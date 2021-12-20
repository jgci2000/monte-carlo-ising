import numpy as np
import matplotlib.pyplot as plt

jdos_exact = np.loadtxt("JDOS_exact_L4_SS.txt")

run_max = 1000
f_exp_vals = np.arange(3, 10 + 1)
flatness = np.array([70, 80, 90, 95, 99])

error = np.zeros((len(flatness), len(f_exp_vals)))
abs_error = np.zeros((len(flatness), len(f_exp_vals)))

for i, p in enumerate(flatness):
    print(f"p = {i}/{len(flatness) - 1}", end='\r')
    for j, f_exp in enumerate(f_exp_vals):
        error_tmp = list()
        abs_error_tmp = list()
        
        for run in range(run_max):
            jdos = np.loadtxt("data/" + str(p) + "/" + str(f_exp) + "/" + str(run) + "_JDOS.txt")
            error_tmp.append(np.sum((jdos[jdos != 0] - jdos_exact[jdos != 0]) / jdos_exact[jdos != 0]))
            abs_error_tmp.append(np.sum(np.abs(jdos[jdos != 0] - jdos_exact[jdos != 0]) / jdos_exact[jdos != 0]))

        error[i, j] = np.mean(error_tmp)
        abs_error[i, j] = np.mean(abs_error_tmp)
        
plt.figure("Absolute Error WL vs f_final")
plt.title(r"Absolute Error in WL vs $f_{final}$")
plt.xlabel(r"$\log(f_{final})$")
plt.ylabel(r"$|\epsilon|$")
plt.ylim([0, np.max(abs_error) * 1.1])
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

plt.show()
