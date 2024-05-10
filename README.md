# AERONastran
@Author: Guillermo Peña Martínez
This code was created to elaborate an assigment for the advanced aeroelasticity course from the Msc in Aeronautical Engineering (UPM). The code is provided as is, and may not work perfectly well in some cases.

Features:
  - Creates a semi - wing mesh.
  - Plot the FEM nodes
  - Create the pannel distribution for DLM
  - Plot the pannels along with the mesh
  - Dump the nodes, pannels and parameters into a .bdf file to be read by Nastran.
  - Generates sol 144 entries.
  - Generates sol 145 entries for PKNL flutter and EIGRL node adquisition.
  - Generates sol 103 and punch file for stiffness and mass matrices.
  - Executes NASTRAN with the input file automatically.
