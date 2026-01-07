# -*- coding: utf-8 -*-
"""
Created on Oct 26 2025
@author: jonas
"""

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

##############################################################
#define test number
test_nr = 'CALCULATION_H6_4A_OPT_OPT_10_OPERATORS_NR_001_100iter'

# not relevant for this application - allows to vary the positions for operator insertion in a tailored way - set 1 and None for this calculation
number_of_positions = 1 #How many postitons for the operator do you want to calculate?
position_arguments = None #-1 positions bevore the last element; if you want to add the operator to the end of the list you need to use the len(ansatz) number as an index for the list

# define imput of old data calculations for input and the directory in which you place the input files
file_name_old_data = "H6_r=4_QE_204i_gradient_sliced_until_23_1.pkl" #"LiH_r=1.5_QE_9i_gradient.pkl"
file_directory = r'/home/jonas/Documents/03_dplts/paper_writing/calculations/codes/optimized_position_app_choice_10operator_correct_v1.0/adaptvqe/algorithms/stored_adapt-vqe_runs/'
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

def analysis(number_of_positions, position_arguments):
    datei="/Documents/Julian_Calculations/Code_Gradient_Thoughs_01/ceo-adapt-vqe/Prepered_Codes_result" #xy.pickle ist die Datei, die die Daten enthält
    #if number_of_positions != len(position_arguments):
    #    print("ERROR: Number of operator postion arrays do not match number entered yourselfes.")
    return position_arguments

def calculation_xy(number_of_positions):
    list_data_evolution_errors=[]
    list_data_evolution_energies=[]
    list_data_evolution_gradient_norms=[]
    list_data_acc_cnot_counts=[]
    list_data_acc_cnot_depths=[]
    list_data_gradient_positions=[]
    list_data_indices = []
    list_data_coefficients = []

    print("CALC1")
    track_prep_g_var_in = True  # True = show prepending gradients
    threshold_in = 5e-7
    r = 4 #distance in angström
    molecule = create_h6(r) #create_h2(r)
    print(molecule)
    pool = QE(molecule)
    pre_app_setting_in = 'a'
    print(file_name_old_data != None)
    if file_name_old_data != None:
        data_old_in = load_data_r(file_name_old_data)
        molecule = None
    else:
        data_old_in = None
        molecule = molecule
    print("CALC2")
    for x in range(0,number_of_positions):
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
        storage_for_names = test_nr,
        )
        
        my_adapt.run()
        data = my_adapt.data
        list_data_evolution_errors.append(data.evolution.errors)
        list_data_evolution_energies.append(data.evolution.energies)    
        list_data_evolution_gradient_norms.append(data.evolution.gradient_norms)
        list_data_acc_cnot_counts.append(data.acc_cnot_counts(pool))
        list_data_acc_cnot_depths.append(data.acc_cnot_depths(pool))
        list_data_gradient_positions.append(data.evolution.pos_gradient_arr)
        list_data_indices.append(data.evolution.indices)
        list_data_coefficients.append(data.evolution.coefficients)
        print("CALC3")
        print(data.evolution.errors)
        print("CALC4")
    dictionary = {
        "data_evolution_errors" : list_data_evolution_errors,
        "data_evolution_energies" : list_data_evolution_energies, 
        "molecule" : molecule,
        "pool" : pool,
        "data_evolution_gradient_norms" : list_data_evolution_gradient_norms,
        "data_acc_cnot_counts" :  list_data_acc_cnot_counts,
        "data_acc_cnot_depths" : list_data_acc_cnot_depths,
        "data_gradient_positions" : list_data_gradient_positions,
        "data_indices" : list_data_indices,
        "data_coefficients" : list_data_coefficients,
        } 
    return dictionary

def plot_with_lines(input_values_y, names_y_functions, function_colors, y_scale_option, shaded_region, x_name, y_name, plot_name, x_def, x_type):
    # Create figure and axis picture size
    fig, ax = plt.subplots(figsize=(3, 2.75))

    if x_type == 'linspace':
        x = []
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        for numb_y in range(len(input_values_y)):
            x.append(np.linspace(x_def[0], x_def[1], x_def[2])) # X-axis scaling
    if x_type == 'values':
        x = x_def
    

    for index_input_val_y in range(len(input_values_y)):
        ax.plot(x[index_input_val_y], input_values_y[index_input_val_y], function_colors[index_input_val_y], linestyle='-', linewidth=0.5, marker='o', markersize=3, label=f"{names_y_functions[index_input_val_y]}")
    
    # different options: linear, log, symlog, asinh, logit, function, functionlog
    # Logarithmic scale for y-axis
    if y_scale_option != False:
        ax.set_yscale(y_scale_option)

    # Shaded region
    if shaded_region:
        ax.fill_between(x[0], 10**-4, 0.0015936193427667, color='blue', alpha=0.2)

    # Legende einfügen
    #ax.legend(fontsize=7, loc='upper right')  # fontsize und Position anpassen

    # Labels and title
    ax.set_xlabel(x_name, fontsize=7)
    ax.set_ylabel(y_name, fontsize=7)
    ax.set_title(plot_name, fontsize=8)

    plt.tight_layout()

    # Show plot
    plt.show()


def plot_with_lines_with_iterations_x(input_values_y, names_y_functions, function_colors, y_scale_option, shaded_region, x_name, y_name, plot_name, x_def, x_type):
    # Create figure and axis picture size
    fig, ax = plt.subplots(figsize=(3, 2.75))

    if x_type == 'linspace':
        x = []
        ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
        for numb_y in range(len(input_values_y)):
            x.append(np.linspace(x_def[0], x_def[1], x_def[2])) # X-axis scaling
    if x_type == 'values':
        x = x_def
    

    for index_input_val_y in range(len(input_values_y)):
        ax.plot(x[index_input_val_y], input_values_y[index_input_val_y], function_colors[index_input_val_y], linestyle='-', linewidth=0.5, marker='o', markersize=3, label=f"{names_y_functions[index_input_val_y]}")
    
    # different options: linear, log, symlog, asinh, logit, function, functionlog
    # Logarithmic scale for y-axis
    if y_scale_option != False:
        ax.set_yscale(y_scale_option)

    # Shaded region
    if shaded_region:
        ax.fill_between(x[0], 10**-4, 0.0015936193427667, color='blue', alpha=0.2)

    # Legende einfügen
    #ax.legend(fontsize=7, loc='upper right')  # fontsize und Position anpassen

    # Schritte auf der x-Achse festlegen (z. B. alle 5 Einheiten)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))

    # increase number size for all axes
    ax.tick_params(axis='both', which='major', labelsize=10)  # Ändert die Schriftgröße der Haupt-Ticks
    ax.tick_params(axis='both', which='minor', labelsize=10)  # Optional: Ändert die Schriftgröße der Neben-Ticks

    # Labels and title
    ax.set_xlabel(x_name, fontsize=10)
    ax.set_ylabel(y_name, fontsize=10)
    ax.set_title(plot_name, fontsize=11)

    plt.tight_layout()

    # Show plot
    plt.show()

def plot_multiple_lines(input_values_y, names_y_functions, function_colors, x_values, x_name, y_name, plot_name, y_scale_option=None):
    """
    Plots multiple lines on the same diagram.

    Parameters:
    - input_values_y: List of lists containing y-values for each line.
    - names_y_functions: List of names for each line (used in the legend).
    - function_colors: List of colors for each line.
    - x_values: List of x-values (shared or individual for each line).
    - x_name: Label for the x-axis.
    - y_name: Label for the y-axis.
    - plot_name: Title of the plot.
    - y_scale_option: Optional scale for the y-axis (e.g., 'log').
    """
    fig, ax = plt.subplots(figsize=(6, 4))

    # Plot each line
    for i, y_values in enumerate(input_values_y):
        ax.plot(x_values[i], y_values, color=function_colors[i], label=names_y_functions[i], marker='o', markersize=4)

    # Set y-axis scale if specified
    if y_scale_option:
        ax.set_yscale(y_scale_option)

    # Add labels, title, and legend
    ax.set_xlabel(x_name, fontsize=10)
    ax.set_ylabel(y_name, fontsize=10)
    ax.set_title(plot_name, fontsize=12)
    ax.legend(fontsize=8)

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()

def save_data(prefix, append):
    # Speichere die beiden Listen als Pickle-Dateien, um die originale Struktur beizubehalten
    with open("/home/jonas/Documents/03_dplts/paper_writing/calculations/codes/optimized_position_app_choice_10operator_correct_v1.0/CALCULATION_H6_4A_OPT_OPT_10_OPERATORS_NR_001_100iter_insert.pkl", "wb") as f:
        cp.dump(append, f)

def load_data(prefix):
    with open("/home/jonas/Documents/03_dplts/paper_writing/calculations/codes/optimized_position_app_choice_10operator_correct_v1.0/CALCULATION_H6_4A_OPT_OPT_10_OPERATORS_NR_001_100iter_insert.pkl", "rb") as f:
        append = cp.load(f)
    return append

def main():
    if (os.path.exists("/home/jonas/Documents/03_dplts/paper_writing/calculations/codes/optimized_position_app_choice_10operator_correct_v1.0/CALCULATION_H6_4A_OPT_OPT_10_OPERATORS_NR_001_100iter_insert.pkl")):
        dictionary = load_data(test_nr)
        print('A')
    else:
        # Falls nicht vorhanden, führe run_adapt aus und speichere die Ergebnisse
        positions = analysis(number_of_positions, position_arguments)
        data_output = calculation_xy(number_of_positions)
        dictionary = data_output
        save_data(test_nr, dictionary)
        print('B')
    print("finished")
    print(dictionary)

    #######################################################
    # print iteration - error, gradient plot
    for index, wert in enumerate(dictionary["data_evolution_gradient_norms"]):
        gradient_norms = dictionary["data_evolution_gradient_norms"][index]
        errors = dictionary["data_evolution_errors"][index]

        # Sicherstellen, dass beide Listen die gleiche Länge haben
        iterations = list(range(1, len(gradient_norms) + 1))  # Startet bei 1 statt 0

        fig, ax1 = plt.subplots()

        # Linke y-Achse für Gradient Norms
        ax1.plot(iterations, gradient_norms, label="Gradient Norm", color="blue")
        ax1.set_yscale("log")
        ax1.set_xlabel("Iterations", fontsize=14)
        ax1.set_ylabel("Gradient Norm", color="blue", fontsize=16)
        ax1.tick_params(axis='y', labelcolor="blue", labelsize=16)
        ax1.tick_params(axis='x', labelsize=16)
        ax1.grid()

        # Rechte y-Achse für Errors
        ax2 = ax1.twinx()
        ax2.plot(iterations, errors, label="Errors", color="red", linestyle="dashed")
        ax2.set_yscale("log")
        ax2.set_ylabel("Errors", color="red", fontsize=16)
        ax2.tick_params(axis='y', labelcolor="red", labelsize=16)

        fig.tight_layout()  # Verhindert Überlappungen
        plt.show()

main()