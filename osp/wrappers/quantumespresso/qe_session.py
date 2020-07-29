from osp.core.session import SimWrapperSession
from osp.core.namespaces import QE
from osp.wrappers.quantumespresso.qe_engine import SimulationEngine
from osp.wrappers.quantumespresso.qe_utils import qeUtils
from osp.core.utils import simple_search

class qeSession(SimWrapperSession):
    
    def __init__(self, root, engine = None, **kwargs):

        # Engine and file utils
        engine = engine or SimulationEngine(self)
        self._qe_utils = qeUtils(self, root)
        super().__init__(engine, **kwargs)

    def __str__(self):
        return "Quantum Espresso Wrapper Session"

    def _run(self, simulation, prefix, command_type = "pw.x", calculation_type = "scf", **kwargs):
        self._prefix = prefix
        self._command_type = command_type
        self._calculation_type = calculation_type

        # Sets input and output files
        self._input_file = f"{self._prefix}.{self._command_type[:-2]}{self._calculation_type}.in"
        self._output_file = f"{self._prefix}.{self._command_type[:-2]}{self._calculation_type}.out"

        # Creates input, runs, and updates the cuds structure
        self._qe_utils._create_input(simulation, **kwargs)
        self._engine.run()
        self._qe_utils._update_cuds(simulation)

    def _load_from_backend(self, uids, expired=None):
        for uid in uids:
            try:
                yield self._registry.get(uid)
            except KeyError:
                yield None

    def _apply_added(self, root_obj, buffer):
        return super()._apply_added(root_obj, buffer)

    def _apply_deleted(self, root_obj, buffer):
        return super()._apply_deleted(root_obj, buffer)

    def _apply_updated(self, root_obj, buffer):
        return super()._apply_updated(root_obj, buffer)
        

    
