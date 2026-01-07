# -*- coding: utf-8 -*-
"""
Created on Oct 26 2025
@author: jonas
"""
#################################################################
# Calculating several adapt-vqe runs in parallel with different #
# old data imputs - using threading or multiprocessing and      #
# pickle files as input and output storage                      #
#################################################################


# Import molecules
from adaptvqe.molecules import create_h6
from adaptvqe.molecules import create_h3
from adaptvqe.molecules import create_h4
from adaptvqe.molecules import create_lih
from adaptvqe.molecules import create_beh2

# Import pools
from adaptvqe.pools import QE
from adaptvqe.pools import FullPauliPool, TiledPauliPool

# Import adapt-vqe algorithm
from adaptvqe.algorithms.adapt_vqe import LinAlgAdapt

# Import additional packages
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import matplotlib.ticker as ticker
import cloudpickle as cp
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from typing import Callable, Dict, List, Any
import traceback
from typing import Callable, Dict, Any, List, Optional



##############################################################
# define calculation/test number
test_nr = 'RORP_H6_4A_100i_parallel_NR_002'

# define imput of old data calculations for input and the directory in which you place the input files
file_names_old_data = ["H6_r=4_QE_204i_gradient_sliced_until_23_1.pkl", "H6_r=4_QE_204i_gradient_sliced_until_23_2.pkl", "H6_r=4_QE_204i_gradient_sliced_until_23_3.pkl", "H6_r=4_QE_204i_gradient_sliced_until_23_4.pkl"] #"LiH_r=1.5_QE_9i_gradient.pkl"
file_directory = r'/home/jonas/Documents/03_dplts/paper_writing/calculations/codes/random_operator_random_pos_CPU_threaded_1.0/adaptvqe/algorithms/stored_adapt-vqe_runs/'
##############################################################

def load_data_r(filename):
    """Imports file data if the file exists."""
    print(file_directory)
    file_path = os.path.join(file_directory, filename)
    if os.path.exists(file_directory):
        with open(file_path, "rb") as f:
            data = pickle.load(f)
        print(f"Data loaded from '{filename}'.")
        print('############################')
        print(data)
        print(type(data))
        print('############################')
        return data
    else:
        print(f"File '{filename}' does not exist. Return: None.")
        return None  # If the file does not exist.


def run_parallel_measurements(
    measurement_fn: Callable[[], Dict[str, Any]],
    n: int,
    use_processes: bool = False,
    max_workers: int | None = None,
    file_names_old_data: Optional[List[str]] = None,
) -> Dict[str, List[Any]]:
    """
    Führt measurement_fn() n-mal parallel aus und sammelt die Rückgabewerte.
    measurement_fn: Funktion, die ein dict zurückgibt (wie in deiner Beschreibung).
    n: Anzahl Durchläufe.
    use_processes: False -> Threads (keine Pickle-Anforderung). True -> Prozesse (gute CPU-Parallelität, benötigt Pickle).
    max_workers: Anzahl paralleler Worker (None -> default).
    Rückgabe: dict mapping keys -> list of run-values (in Ausführungsreihenfolge).
    Zusätzlich wird bei Fehlern eine Liste 'errors' mit (index, exception_str) angelegt.
    """
    executor_cls = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
    results: List[Dict[str, Any]] = []
    errors: List[Dict[str, Any]] = []

    # Submit tasks
    with executor_cls(max_workers=max_workers) as ex:
        #futures = {ex.submit(measurement_fn): i for i in range(n)}
        futures = {ex.submit(measurement_fn, name): i for i, name in enumerate(file_names_old_data)}

        for fut in as_completed(futures):
            idx = futures[fut]
            try:
                res = fut.result()
                if not isinstance(res, dict):
                    errors.append({"index": idx, "error": f"Expected dict, got {type(res)}", "value": res})
                else:
                    results.append((idx, res))
            except Exception as e:
                tb = traceback.format_exc()
                errors.append({"index": idx, "error": str(e), "traceback": tb})

    # Sort results by original submission order (optional; as_completed is unordered)
    results.sort(key=lambda x: x[0])  # sort by index
    run_dicts = [r for (_, r) in results]

    # Merge into one dict of lists
    merged: Dict[str, List[Any]] = {}
    if run_dicts:
        # initialize keys from first run
        for k in run_dicts[0].keys():
            merged[k] = []
        # append each run's value (don't extend: each run's value stays as one element)
        for rd in run_dicts:
            for k, v in rd.items():
                merged[k].append(v)

    # Attach errors info
    merged["_run_errors"] = errors  # leer, falls keine Fehler

    return merged

def set_settings():
    global track_prep_g_var_in
    global threshold_in
    global r
    global molecule
    global pool
    global pre_app_setting_in
    global data_old_in
    global data_old_in
    track_prep_g_var_in = True  # True = show prepending gradients
    threshold_in = 5e-7
    r = 4 #distance in angström
    molecule = create_h6(r) #create_h2(r)
    pool = QE(molecule)
    pre_app_setting_in = 'a'

def measurement_random(name_old_data_in):
    
    if name_old_data_in != None:
        data_old_in = load_data_r(name_old_data_in)
        molecule = None
        molecule_for_pool = create_h6(r)
        pool = QE(molecule_for_pool)
    else:
        data_old_in = None
        molecule = create_h6(r)

    
    my_adapt=LinAlgAdapt(
        pool=pool,
        molecule=molecule,
        max_adapt_iter=100,
        recycle_hessian=False,
        tetris=False,
        verbose=True,
        threshold=threshold_in,
        track_prep_g=track_prep_g_var_in,
        pre_app_setting = pre_app_setting_in,
        previous_data = data_old_in,
        operator_pos_arr = None,
        storage_for_names = name_old_data_in
        #all_pos_gradient_arr = [] # This is the array that will be filled with the gradients of all positions
        )
    my_adapt.run()

    data = my_adapt.data
    # simuliert z.B. Messwert + Statuscode
    return {
        "data_evolution_errors" : data.evolution.errors,
        "data_evolution_energies" : data.evolution.energies, 
        "molecule" : molecule,
        "pool" : pool,
        "data_evolution_gradient_norms" : data.evolution.gradient_norms,
        "data_acc_cnot_counts" : data.acc_cnot_counts(pool),
        "data_acc_cnot_depths" : data.acc_cnot_depths(pool),
        "data_gradient_positions" : data.evolution.pos_gradient_arr,
        "data_indices" : data.evolution.indices,
        "data_coefficients" : data.evolution.coefficients,
    }



def save_data(prefix, append):
    # Speichere die beiden Listen als Pickle-Dateien, um die originale Struktur beizubehalten
    with open("/home/jonas/Documents/03_dplts/paper_writing/calculations/codes/random_operator_random_pos_CPU_threaded_1.0/RORP_H6_4A_100i_parallel_NR_002_insert.pkl", "wb") as f:
        cp.dump(append, f)

def load_data(prefix):
    with open("/home/jonas/Documents/03_dplts/paper_writing/calculations/codes/random_operator_random_pos_CPU_threaded_1.0/RORP_H6_4A_100i_parallel_NR_002_insert.pkl", "rb") as f:
        append = cp.load(f)
    return append

def main():
    if (os.path.exists("/home/jonas/Documents/03_dplts/paper_writing/calculations/codes/random_operator_random_pos_CPU_threaded_1.0/RORP_H6_4A_100i_parallel_NR_002_insert.pkl")):
        dictionary = load_data(test_nr)
        print('A')
    else:
        # Falls nicht vorhanden, führe run_adapt aus und speichere die Ergebnisse
        set_settings()
        data_output = run_parallel_measurements(measurement_random, n=len(file_names_old_data), use_processes=False, max_workers=15, file_names_old_data=file_names_old_data)
        dictionary = data_output
        save_data(test_nr, dictionary)
        print('B')
    print("finished")
    print(dictionary)


    #######################################################
    # Print iteration - error, gradient plot
    for index, wert in enumerate(dictionary["data_evolution_gradient_norms"]):
        gradient_norms = dictionary["data_evolution_gradient_norms"][index]
        errors = dictionary["data_evolution_errors"][index]

        # Ensure that both lists are the same length.
        iterations = list(range(1, len(gradient_norms) + 1))  # Starts at 1 instead of 0

        fig, ax1 = plt.subplots()

        # Left y-axis for gradient norms
        ax1.plot(iterations, gradient_norms, label="Gradient Norm", color="blue")
        ax1.set_yscale("log")
        ax1.set_xlabel("Iterations", fontsize=14)
        ax1.set_ylabel("Gradient Norm", color="blue", fontsize=16)
        ax1.tick_params(axis='y', labelcolor="blue", labelsize=16)
        ax1.tick_params(axis='x', labelsize=16)
        ax1.grid()

        # Right y-axis for errors
        ax2 = ax1.twinx()
        ax2.plot(iterations, errors, label="Errors", color="red", linestyle="dashed")
        ax2.set_yscale("log")
        ax2.set_ylabel("Errors", color="red", fontsize=16)
        ax2.tick_params(axis='y', labelcolor="red", labelsize=16)

        fig.tight_layout()  # Prevents overlaps
        plt.show()

    
    
    
    
    ###########################################################

    # Ensure that both lists are the same length.
    iterations = list(range(1, len(gradient_norms) + 1))  # Starts at 1 instead of 0

    fig, ax1 = plt.subplots()

    # Left y-axis for gradient norms
    for index, wert in enumerate(dictionary["data_evolution_gradient_norms"]):
        gradient_norms = dictionary["data_evolution_gradient_norms"][index]
        ax1.plot(iterations, gradient_norms, label="Gradient Norm", color="blue")
    ax1.set_yscale("log")
    ax1.set_xlabel("Iterations", fontsize=14)
    ax1.set_ylabel("Gradient Norm", color="blue", fontsize=16)
    ax1.tick_params(axis='y', labelcolor="blue", labelsize=16)
    ax1.tick_params(axis='x', labelsize=16)
    ax1.grid()

    # Right y-axis for errors
    errors = dictionary["data_evolution_errors"][index]
    for index, wert in enumerate(dictionary["data_evolution_gradient_norms"]):
        ax2 = ax1.twinx()
        ax2.plot(iterations, errors, label="Errors", color="red", linestyle="dashed")
    ax2.set_yscale("log")
    ax2.set_ylabel("Errors", color="red", fontsize=16)
    ax2.tick_params(axis='y', labelcolor="red", labelsize=16)

    fig.tight_layout()  # Prevents overlaps
    plt.show()

main()








































