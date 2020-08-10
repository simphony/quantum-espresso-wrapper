import numpy as np 

from osp.core.namespaces import QE
from osp.core.utils import pretty_print
from osp.wrappers.quantumespresso.qe_session import qeSession

sim = QE.Simulation()
k = QE.K_POINTS(vector = (7, 7, 7), unit = "")

SiCell = QE.Cell()
Si = QE.Element(name = "Si")
SiPseudo = QE.PSEUDOPOTENTIAL(name = "Si.pbe-n-kjpaw_psl.1.0.0.UPF")
Si1 = QE.Atom()
SiParams = QE.CellParams()
celldm1 = QE.Celldm1(value = 5.43070, unit = "au")

Si.add(SiPseudo, Si1)
Si.add(QE.Mass(value = 28.085, unit = "amu"))
SiCell.add(Si1, SiParams)
Si1.add(QE.Position(vector = (0, 0, 0), unit = ""))
SiCell.add(celldm1)

SiParams.add(QE.CellParameterX(vector = (0.5, 0.5, 0), unit = ""),
             QE.CellParameterY(vector = (0.5, 0, 0.5), unit = ""),
             QE.CellParameterZ(vector = (0, 0.5, 0.5), unit = ""))

sim.add(SiCell)
sim.add(Si)
sim.add(k)

sim.add(QE.Pressure(value = 100, unit = "kbar"))
sim.add(QE.StressTensor(tensor2 = np.zeros((3, 3)), unit = "kbar"))

pretty_print(sim)

session = qeSession(root = "../../../")
quantum_espresso_wrapper = QE.QEWrapper(session = session)
quantum_espresso_wrapper.add(sim)
print("Running calculation...")

quantum_espresso_wrapper.session._run(simulation = sim, prefix = "si", command_type = "pw.x", calculation_type = "scf", root = "/mnt/c/iwm/quantum-espresso-wrapper/", CONTROL = {'pseudo_dir': "'.'"})
