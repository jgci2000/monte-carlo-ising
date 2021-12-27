import numpy as np
import matplotlib.pyplot as plt

jdos_exact = np.loadtxt("JDOS_exact_L4_SS.txt")

run_max = 1000
rep_exp_vals = np.arange(2, 7 + 1)
skip_exp_vals = np.array([-2, -1, 0, 1, 2])

error = np.zeros((len(skip_exp_vals), len(rep_exp_vals)))
abs_error = np.zeros((len(skip_exp_vals), len(rep_exp_vals)))
steps = np.zeros((len(skip_exp_vals), len(rep_exp_vals))) 
wall_time = np.zeros((len(skip_exp_vals), len(rep_exp_vals)))

for i, rep_exp in enumerate(rep_exp_vals):    
    print(f"rep = {i}/{len(rep_exp_vals) - 1}", end='\r')
    for j, skip_exp in enumerate(skip_exp_vals):
        error_tmp = list()
        abs_error_tmp = list()
        
        steps_tmp = list()
        wall_time_tmp = list()
        
        for run in range(run_max):
            jdos = np.loadtxt("data/REP_" + str(rep_exp) + "/skip_" + str(skip_exp) + "/" + str(run + 1) + "_JDOS.txt")
            error_tmp.append(np.sum((jdos[jdos != 0] - jdos_exact[jdos != 0]) / jdos_exact[jdos != 0]))
            abs_error_tmp.append(np.sum(np.abs(jdos[jdos != 0] - jdos_exact[jdos != 0]) / jdos_exact[jdos != 0]))
            
            aux = np.loadtxt("data/REP_" + str(rep_exp) + "/skip_" + str(skip_exp) + "/debug/" + str(run + 1) + "_JDOS.txt", skiprows=1) 
            steps_tmp.append(np.sum(aux[:, 0]))
            wall_time_tmp.append(np.sum(aux[:, 1]))
            

        error[j, i] = np.mean(error_tmp)
        abs_error[j, i] = np.mean(abs_error_tmp)
        
        steps[j, i] = np.mean(steps_tmp)
        wall_time[j, i] = np.mean(wall_time_tmp)

plt.figure("Absolute Error FSS vs REP")
plt.title(r"Absolute Error in FSS vs $REP$")
plt.xlabel(r"$\log(REP)$")
plt.ylabel(r"$|\epsilon|$")
plt.ylim([0, np.max(abs_error) * 1.1])
for i in range(len(skip_exp_vals)):
    plt.plot(rep_exp_vals, abs_error[i, :], 'o-', label=f"$skip={skip_exp_vals[i]}$")
plt.legend()

plt.figure("Error FSS vs REP")
plt.title(r"Error in FSS vs $REP$")
plt.xlabel(r"$\log(REP)$")
plt.ylabel(r"$\epsilon$")
for i in range(len(skip_exp_vals)):
    plt.plot(rep_exp_vals, error[i, :], 'o-', label=f"$skip={skip_exp_vals[i]}$")
plt.legend()

plt.figure("Wall Time FSS vs REP")
plt.title("Wall Time (s) in FSS vs $REP$")
plt.xlabel(r"$\log(REP)$")
plt.ylabel("Wall Time (s)")
for i in range(len(skip_exp_vals)):
    plt.plot(rep_exp_vals, wall_time[i, :], 'o-', label=f"$skip={skip_exp_vals[i]}$")
plt.legend()

plt.figure("Absolute Error FSS vs skip")
plt.title(r"Absolute Error in FSS vs $skip$")
plt.xlabel(r"$skip$")
plt.ylabel(r"$|\epsilon|$")
plt.ylim([0, np.max(abs_error) * 1.1])
for i in range(len(rep_exp_vals)):
    plt.plot(skip_exp_vals, abs_error[:, i], 'o-', label=f"$rep={rep_exp_vals[i]}$")
plt.legend()

plt.figure("Error FSS vs skip")
plt.title(r"Error in FSS vs $skip$")
plt.xlabel(r"$skip$")
plt.ylabel(r"$\epsilon$")
for i in range(len(rep_exp_vals)):
    plt.plot(skip_exp_vals, error[:, i], 'o-', label=f"$rep={rep_exp_vals[i]}$")
plt.legend()

plt.figure("Wall Time FSS vs skip")
plt.title("Wall Time (s) in FSS vs $skip$")
plt.xlabel(r"$skip$")
plt.ylabel("Wall Time (s)")
for i in range(len(rep_exp_vals)):
    plt.plot(skip_exp_vals, wall_time[:, i], 'o-', label=f"$rep={rep_exp_vals[i]}$")
plt.legend()

plt.show()
