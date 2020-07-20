import numpy as np
from osp.core.namespaces import QE
from osp.core.utils import pretty_print
from osp.wrappers.quantumespresso.qe_session import qeSession

# Creates simulation
sim = QE.Simulation()
k = QE.K_POINTS(vector = (7, 7, 7), unit = "")

# Creates a cell, the element Silicon, a pseudopotential, two atoms and cell parameters
SiCell = QE.Cell()
Si = QE.Element(name = "Si")
SiPseudo = QE.PSEUDOPOTENTIAL(name = "Si.pbe-n-kjpaw_psl.1.0.0.UPF")
Si1 = QE.Atom()
SiParams = QE.CellParams()
celldm1 = QE.Celldm1(value = 5.43070, unit = "au")

# Adds pseudopotential and atoms to the element
# Describes element's mass
# Adds atoms and cell parameters to the cell
# Positions the atoms
Si.add(SiPseudo, Si1)
Si.add(QE.Mass(value = 28.085, unit = "amu"))
SiCell.add(Si1, SiParams)
Si1.add(QE.Position(vector = (0, 0, 0), unit = ""))
SiCell.add(celldm1)

# Specifies the values of the cell parameters
SiParams.add(QE.CellParameterX(vector = (0.5, 0.5, 0), unit = ""),
             QE.CellParameterY(vector = (0.5, 0, 0.5), unit = ""),
             QE.CellParameterZ(vector = (0, 0.5, 0.5), unit = ""))

# Adds cell and element to simulation
sim.add(SiCell)
sim.add(Si)
sim.add(k)
sim.add(QE.Pressure(value = 100, unit = "kbar"))
sim.add(QE.StressTensor(tensor2 = np.zeros((3, 3)), unit = "kbar"))
root = ""
pretty_print(sim)
with qeSession(root) as session:
    # Adds session to wrapper
    quantum_espresso_wrapper = QE.QEWrapper(session = session)
    # Adds simulation to wrapper
    sim = quantum_espresso_wrapper.add(sim)
    # pretty_print(sim)
    # Creates a qeUtil object and creates an input file based off of the simulation
    print("Running calculation...")
    
    # Runs the simulation
    quantum_espresso_wrapper.session._run(prefix = "si", command_type = "pw.x", calculation_type = "scf")
    quantum_espresso_wrapper.session._run(prefix = "si", command_type = "pw.x", calculation_type = "bands")
    quantum_espresso_wrapper.session._run(prefix = "si", command_type = "bands.x", calculation_type = "")
    quantum_espresso_wrapper.session._run(prefix = "si", command_type = "pw.x", calculation_type = "relax", IONS = {'ion_dynamics': "'bfgs'"})

    # print("Results: ")
    # Pretty prints the simulation
    pretty_print(sim)