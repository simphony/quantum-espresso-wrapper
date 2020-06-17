from osp.core import Parser
p = Parser()
p.parse("ontology.qe.yml")
from osp.core import QE
from osp.core.utils import pretty_print
from osp.wrappers.quantumespresso.qe_session import qeSession


sim = QE.Simulation()

# Create the cell
SiCell = QE.Cell()

Si = QE.Element()
SiP = QE.Pseudopotential()
# Create the atoms
Si1 = QE.Atom()
Si1.add(QE.Position(vector=(0,0,0), unit = ""),
        QE.Mass(value=28.08500, unit = "amu"))

Si2 = QE.Atom()
Si2.add(QE.Position(vector=(0.25,0.25,0.25), unit = ""),
        QE.Mass(value=28.08500, unit = "amu"))

SiCell.add(Si1, Si2)

sim.add(SiCell)

with qeSession() as session:
    quantum_espresso_wrapper = QE.Qe_wrapper(session = session)
    sim = quantum_espresso_wrapper.add(sim)

    pretty_print(sim)
    print("Runnng calculation...")
    quantum_espresso_wrapper.session.run()
    print("Results: ")
