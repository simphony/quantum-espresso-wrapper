from osp.core.session import WrapperSession
from osp.core import QE
from osp.wrappers.quantumespresso.qe_engine import SimulationEngine
from osp.wrappers.quantumespresso.qe_utils import qeUtils


class qeSession(WrapperSession):
    
    def __init__(self, engine = None, **kwargs):
        #File names
        self._infile = "inputfile.in"

        engine = engine or SimulationEngine(self)
        self._qe_utils = qeUtils(self)
        super().__init__(engine, **kwargs)


    def __str__(self):
        return "Quantum Espresso Wrapper Session"

    def _apply_added(self, root_obj, buffer):
        return super()._apply_added(root_obj, buffer)

    def _apply_updated(self, root_obj, buffer):
        return super()._apply_updated(root_obj, buffer)

    def _apply_deleted(self, root_obj, buffer):
        return super()._apply_deleted(root_obj, buffer)

    def _load_from_backend(self, uids, expired=None):
        for uid in uids:
            if uid in self._registry:
                yield self._registry[uid]
            else:
                yield None
    
    def run(self, input_files, output_files):
        root = self._registry.get(self.root)
        simulation = root.get(oclass = QE.Simulation)

        self._engine.run(input_files = input_files, output_files = output_files)

    
