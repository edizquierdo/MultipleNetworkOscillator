Scripts to simulate recurrent subcircuits inside the VNC unit
main.cpp in script folder simulate independently the three subcircuits, 
Interunit connections are ablated in Worm.cpp file
In Worm.cpp there is also a modified timestep function that only update the nervios_system part of the model. So, no body variables are updated.

Each run of main.cpp receive three values from terminal ('./main v1 v2 v3') that goes as external inputs to the 3 neurons in each recurrent subcircuit.
The output is the total amount of change in derivative normalizad by simulation time (20s after 10s of adaptation) and number of neurons in the circuits (3).
The output is a 3 value vector, each value correspond to mean change in derivative per second per neuron in each subcircuit.

simulate_and_plot.py execule main for each selected indiviauls in a gris of 9x9x9 external input values.
Select the maximum value for each subcircuit (maximal oscillation) and generate a plot.

