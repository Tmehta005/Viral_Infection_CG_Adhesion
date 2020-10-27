# viral infection coarse grain adhesion model

## Background

This code provides a multi-scale modelling framework based on the biophysical mechanisms associated with the attachment of the virus with the host cell.  The unified approach helps to characterize the time scales associated with cellular mechanisms associated with the virus attack. This code uses coarse grained dynamic model for virus host cell interactions during cellular adhesion and entry. The square section of the membrane domain is represented as the beads connected with springs. Each of the bead is connected to its 6 neighbors, except for the ones at the end. Virus is represented as an analytical spherical surface.

## Code Description

spherical_100.py: It is the driver for the case of spherical virions, like Corona along with the Brownian dynamics calculations to describe the folding of the membrane around the virus. This file reads the inputs from an input file membrane_100.xlsx and imports the functions from graph_def.py. This main program performs the calculations at 5 different Temperatures (K) (280, 290, 300, 310, 320) and six different values of membrane tension (pN) (0, 50, 100, 150, 200, 250) and writes the output in the locationout_3D.xls 

## Output

Adhesion dynamics and effectiveness are calculated based on the large amounts of bead location data. The output file has two different results - (a) the time-dependent average proximity of the membrane from the virus for each of the Temperature-Tension combinations (b) The steady state values of the proximity index and time to attach.


