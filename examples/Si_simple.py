from osp.core import Parser
p = Parser()
p.parse("ontology.qe.yml")
from osp.core import QE
from osp.core.utils import pretty_print
from osp.wrappers.quantumespresso.qe_session import qeSession
from osp.wrappers.quantumespresso.qe_utils import qeUtils


sim = QE.Simulation()

SiCell = QE.Cell()
Si = QE.Element(name = "Si")
SiPseudo = QE.Pseudopotential(name = "Si.pbe-n-kjpaw_psl.1.0.0.UPF")
Si1 = QE.ATOM()
Si2 = QE.ATOM()
SiParams = QE.CellParams()

Si.add(SiPseudo, Si1, Si2)
Si.add(QE.Mass(value = 28.085, unit = "amu"))
SiCell.add(Si1, Si2, SiParams)
Si1.add(QE.Position(vector = (0, 0, 0), unit = ""))
Si2.add(QE.Position(vector = (0.25, 0.25, 0.25), unit = ""))

SiParams.add(QE.CellParameterX(vector = (0.5, 0.5, 0), unit = ""),
             QE.CellParameterY(vector = (0.5, 0, 0.5), unit = ""),
             QE.CellParameterZ(vector = (0, 0.5, 0.5), unit = ""))

sim.add(SiCell)
sim.add(Si)

with qeSession() as session:
    quantum_espresso_wrapper = QE.Qe_wrapper(session = session)
    sim = quantum_espresso_wrapper.add(sim)

    pretty_print(sim)
    print("Runnng calculation...")
    util = qeUtils(qeSession)
    util._create_input('testinput.in', sim)
    quantum_espresso_wrapper.session.run({"-input": "testinput.in"}, {">": "basicSiCrystal.out"})
    sim.add(QE.TotalEnergy(value = util._read_output('basicSiCrystal.out'), unit = "Ry"))
#     quantum_espresso_wrapper.session.run()
    print("Results: ")
    pretty_print(sim)
