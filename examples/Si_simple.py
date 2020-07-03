from osp.core import Parser
p = Parser()
p.parse("ontology.qe.yml")
from osp.core import QE
from osp.core.utils import pretty_print
from osp.wrappers.quantumespresso.qe_session import qeSession
from osp.wrappers.quantumespresso.qe_utils import qeUtils

# Creates simulation
sim = QE.Simulation()
k = QE.K_points(vector = (7, 7, 7), unit = "")

# Creates a cell, the element Silicon, a pseudopotential, two atoms and cell parameters
SiCell = QE.Cell()
Si = QE.Element(name = "Si")
SiPseudo = QE.Pseudopotential(name = "Si.pbe-n-kjpaw_psl.1.0.0.UPF")
Si1 = QE.ATOM()
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

with qeSession() as session:
    # Adds session to wrapper
    quantum_espresso_wrapper = QE.Qe_wrapper(session = session)
    # Adds simulation to wrapper
    sim = quantum_espresso_wrapper.add(sim)
    # pretty_print(sim)
    # Creates a qeUtil object and creates an input file based off of the simulation
    print("Runnng calculation...")
    util = qeUtils(qeSession)
    util._create_input('testinput.in', sim, CONTROL = {"outdir": "'.'"})
    
    # Runs the simulation
    quantum_espresso_wrapper.session.run({"-input": "testinput.in"}, {">": "basicSiCrystal.out"})
    
    #Adds output to simulation
    sim.add(QE.TotalEnergy(value = util._read_output('basicSiCrystal.out'), unit = "Ry"))

    # print("Results: ")
    # Pretty prints the simulation
    pretty_print(sim)
#%%
