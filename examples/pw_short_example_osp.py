import numpy as np
from osp.core.namespaces import QE
from osp.wrappers.quantumespresso.qe_session import qeSession
from osp.core.utils import pretty_print

sim  = QE.Simulation()
root = ""
session = qeSession(root)
quantum_espresso_wrapper = QE.QEWrapper(session = session)
quantum_espresso_wrapper.add(sim)

cell = QE.Cell()
alat = 4.
cellParams = cell.add(QE.CellParams(tensor2 = 
[[1., 0., 0.,],
[0., 1., 0.,],
[0., 0., 1.,]], unit = "alat"))
cell.add(QE.Celldm1(value = alat, unit = ""))

O = QE.Element(name = "O")
Ba = QE.Element(name = "Ba")
Ti = QE.Element(name = "Ti")

O.add(QE.Mass(value = 15.999, unit = "amu"))
Ba.add(QE.Mass(value = 137.327, unit = "amu"))
Ti.add(QE.Mass(value = 47.867, unit = "amu"))

O.add(QE.PSEUDOPOTENTIAL(path = "O.pbe-n-kjpaw_psl.1.0.0.UPF"))
Ba.add(QE.PSEUDOPOTENTIAL(path = "Ba.pbe-spn-kjpaw_psl.1.0.0.UPF"))
Ti.add(QE.PSEUDOPOTENTIAL(path = "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF"))

O1 = O.add(QE.Atom())
O2 = O.add(QE.Atom())
O3 = O.add(QE.Atom())
Ba1 = Ba.add(QE.Atom())
Ti1 = Ti.add(QE.Atom())

O1.add(QE.Position(vector = [0.5, 0.5, 0.], unit = ""))
O2.add(QE.Position(vector = [0.5, 0., 0.5], unit = ""))
O3.add(QE.Position(vector = [0., 0.5, 0.5], unit = ""))
Ba1.add(QE.Position(vector = [0., 0., 0.], unit = ""))
Ti1.add(QE.Position(vector = [0.5, 0.5, 0.5], unit = ""))

cell.add(O1, O2, O3, Ba1, Ti1)

kpoints = QE.K_POINTS(vector6 = (4, 4, 4, 0, 0, 0), unit = "automatic")

sim.add(cell, O, Ba, Ti, kpoints)

paramdict = {
        'CONTROL': {
            'calculation': 'scf',
            'restart_mode': 'from_scratch',
            'wf_collect': '.true.',
        },
        'SYSTEM': {
            'ecutwfc': 30.,
            'ecutrho': 240.,
        },
        'ELECTRONS': {
            'conv_thr': 1.e-6,
        }
    }

pretty_print(sim)
session._run(simulation = sim, prefix = "BaTiO3", command_type="pw.x", calculation_type="scf", root = root, **paramdict)
pretty_print(sim)
