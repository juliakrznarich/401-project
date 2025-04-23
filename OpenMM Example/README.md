

# OpenMM Evaluation

This directory contains materials related to the evaluation of OpenMM, an open-source toolkit for molecular simulation.

**Description:**

OpenMM is a powerful and highly extensible library used in computational biochemistry, biophysics, and materials science to simulate the movement of atoms and molecules. It provides a flexible framework for setting up and running molecular dynamics (MD) simulations, energy minimization, and other related calculations. OpenMM is designed to leverage the power of modern hardware, including GPUs, to achieve high performance in these computationally intensive tasks.

**Use in Science and Engineering:**

In scientific research and engineering, OpenMM is crucial for:

* **Understanding Biomolecular Systems:** Researchers use it to study the dynamics of proteins (like myoglobin, as in this project), nucleic acids, lipids, and carbohydrates at the atomic level. This helps in understanding protein folding, enzyme mechanisms, drug binding, and other essential biological processes.
* **Drug Discovery and Development:** MD simulations with OpenMM can be used to screen potential drug candidates, understand their binding interactions with target proteins, and optimize their properties.
* **Materials Science:** OpenMM can simulate the behavior of various materials, aiding in the design and development of new materials with desired properties.
* **Method Development:** Its open-source nature allows researchers to extend its capabilities by implementing new force fields, algorithms, and analysis tools.

**Classification:**

OpenMM is best classified as a **programming tool** or a **library**. While it can be used through Python scripts and other programming interfaces to build applications, it is not a standalone end-user application in itself. It provides the building blocks and functionalities that scientists and engineers use to create and run their simulations. It acts as a foundation upon which more specialized simulation workflows and analysis pipelines can be built.

**Installation**

1. Deactivate any currently active environment (optional, but good practice)
conda deactivate

2. Navigate to your home directory or project directory
cd ~

3. Create a new Conda environment 
conda create --name openmm_env python=3.10  # Or your preferred Python version

4. Activate your new environment
conda activate openmm_env

5. Install OpenMM from the conda-forge channel within your active environment
conda install -c conda-forge openmm

# Example: Protein in Water

**Download the Protein Structure**

wget https://files.rcsb.org/download/1AKI.pdb

**create a .py with the following contents**
**Load the PDB file into OpenMM**

from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

pdb = PDBFile("1AKI.pdb")

**Define the forcefield**

forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

**Cleanup**

modeller = Modeller(pdb.topology, pdb.positions)
modeller.deleteWater()
residues=modeller.addHydrogens(forcefield)

**Solvate**

modeller.addSolvent(forcefield, padding=1.0*nanometer)

**Setup system and Integrator**

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

**Local energy minimization**

print("Minimizing energy")
simulation.minimizeEnergy()

**Setup reporting**

simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(StateDataReporter("md_log.txt", 100, step=True,
        potentialEnergy=True, temperature=True, volume=True))

**NVT equillibration**

print("Running NVT")
simulation.step(10000)

**NPT production MD**

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
simulation.context.reinitialize(preserveState=True)


print("Running NPT")
simulation.step(10000)

**Analysis**

import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("md_log.txt", delimiter=',')

step = data[:,0]
potential_energy = data[:,1]
temperature = data[:,2]
volume = data[:,3]

plt.plot(step, potential_energy)
plt.xlabel("Step")
plt.ylabel("Potential energy (kJ/mol)")
plt.show()
plt.plot(step, temperature)
plt.xlabel("Step")
plt.ylabel("Temperature (K)")
plt.show()
plt.plot(step, volume)
plt.xlabel("Step")
plt.ylabel("Volume (nm^3)")
plt.show()

 **References**
 https://openmm.github.io/openmm-cookbook/latest/notebooks/tutorials/protein_in_water.html