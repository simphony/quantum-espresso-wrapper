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
SiCell.add(QE.Volume(value = 22, unit = "au^3"))
sim.add(QE.TotalEnergy(value = -434, unit = "Ry"))
q = QE.QPoint(vector = (0, 0, 0), unit = "", calculate = True)
sim.add(q)
q.add(QE.Mode(number = 3))
q.add(QE.Mode(number = 2))
q.add(QE.Mode(number = 1))

sim2 = QE.Simulation()
fd = QE.Cell()
sim2.add(fd)

fd.add(QE.Volume(value = 33, unit = "au^3"))
sim2.add(QE.TotalEnergy(value = -432, unit = "Ry"))
 
with qeSession(root) as session:
    # Adds session to wrapper
    quantum_espresso_wrapper = QE.QEWrapper(session = session)
    # Adds simulation to wrapper
    sim = quantum_espresso_wrapper.add(sim)
    # pretty_print(sim)
    # Creates a qeUtil object and creates an input file based off of the simulation
    print("Running calculation...")
    # Runs the simulation
    # pretty_print(quantum_espresso_wrapper)
    # quantum_espresso_wrapper.session._run(simulation = sim, prefix = "si", command_type = "pw.x", calculation_type = "scf")
    # quantum_espresso_wrapper.session._run(simulation = sim, prefix = "si", command_type = "pw.x", calculation_type = "bands")
    # quantum_espresso_wrapper.session._run(simulation = sim, prefix = "si", command_type = "bands.x", calculation_type = "")
    # quantum_espresso_wrapper.session._run(simulation = sim, prefix = "si", command_type = "pw.x", calculation_type = "relax", IONS = {'ion_dynamics': "'bfgs'"})
    # quantum_espresso_wrapper.session._run(simulation = sim, prefix = "si", command_type = "pw.x", calculation_type = "scf", SYSTEM = {'occupations': "'tetrahedra'"})
    # quantum_espresso_wrapper.session._run(simulation = sim, prefix = "si", command_type = "dos.x", calculation_type = "")
    # quantum_espresso_wrapper.session._run(simulation = sim, prefix = "si", command_type = "pp.x", calculation_type = "9")
    # quantum_espresso_wrapper.session._run(simulation = [sim, sim2], prefix = 'si', command_type = "ev.x", calculation_type = '1')
    quantum_espresso_wrapper.session._run(simulation = sim, prefix = "si", command_type = "ph.x", calculation_type = "")
    
    pretty_print(sim)
    # pretty_print(sim2)
    # print("Results: ")
    # Pretty prints the simulation

