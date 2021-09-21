"""
This is a sample simulation that does not represent any particular biological system. It is just a showcase 
of how create a Simulation object, add forces, and initialize the reporter. 

In this simulation, a simple polymer chain of 10,000 monomers is 
"""


import os, sys
import polychrom
from polychrom import simulation, starting_conformations, forces, forcekits
import openmm
import os
from polychrom.hdf5_format import HDF5Reporter
import numpy as np

#N=10000
N=100

#np.random.seed(10)
#reporter = HDF5Reporter(folder="trajectory", max_data_length=5, overwrite=True) # check: hdf5_format.py

sim = simulation.Simulation(
    platform="CPU", 
    #integrator="verlet",
    integrator="verlet",
    #integrator="variableLangevin", #se metto questo devo mettere errot_tol e collision rate al posto di timestep
    #error_tol=0.003,
    GPU="0",
    timestep=0.01,
    #timestep=0.01, #used for verlet
    collision_rate=0.03, 
    N=N,
    save_decimals=5,
    PBCbox=False,
    #reporters=[reporter],
    #temperatura di default a 300K
)

#polymer = starting_conformations.grow_cubic(N, 2, method="linear") # This function grows a ring or linear polymer 
                                                        # on a cubic lattice 
                                                        # in the cubic box of size boxSize. 

polymer = starting_conformations.create_random_walk(1.1, N)
sim.set_data(polymer, center=True)  # loads a polymer, puts a center of mass at zero

#sim.add_force(forces.spherical_confinement(sim, density=0.85, k=1))

sim.add_force( #anche add_forceandtorque
    forcekits.polymer_chains(
        sim,
        chains=[(0, None, False)],
        # By default the library assumes you have one polymer chain
        # If you want to make it a ring, or more than one chain, use self.setChains
        # self.setChains([(0,50,True),(50,None,False)]) will set a 50-monomer ring and a chain from monomer 50 to the end
        bond_force_func=forces.harmonic_bonds,
        bond_force_kwargs={
            "bondLength": 1.0,
            "bondWiggleDistance": 0.1,  # Bond distance will fluctuate +- 0.1 on average
        },
        angle_force_func=None,
        angle_force_kwargs={
            "k": 0.0,
            # K is more or less arbitrary, k=4 corresponds to presistence length of 4,
            # k=1.5 is recommended to make polymer realistically flexible; k=8 is very stiff
        },
        #nonbonded_force_func=forces.polynomial_repulsive, #inserire grosberg
        #nonbonded_force_kwargs={
        #    "trunc": 10.0,  # this will let chains cross sometimes
            #'trunc':10.0, # this will resolve chain crossings and will not let chain cross anymore
        #},
        grosberg_force_func=forces.grosberg_repulsive_force,
        grosberg_force_kwargs={
            "trunc": None, "radiusMult": 1.0,
        },
        except_bonds=True,
    )
)
sim.local_energy_minimization()

for _ in range(1000):  # Do 10 blocks. meglio 1000 blocchi da uno step ciascuno
    sim.do_block(100)  # Of 100 timesteps each. Data is saved automatically. 
sim.print_stats()  # In the end, print very simple statistics

#reporter.dump_data()  # always need to run in the end to dump the block cache to the disk
                      # check: hdf5_format.py

#to convert to ascii: h5dump -o file_name.asci -y -w 400 file_name.h5
