from osp.core.session import SimWrapperSession
from osp.core import QE
from osp.wrappers.quantumespresso.qe_engine import SimulationEngine
from osp.wrappers.quantumespresso.qe_utils import qeUtils
from osp.core.utils import simple_search


class qeSession(SimWrapperSession):
    
    def __init__(self, root, engine = None, **kwargs):
        #File names
        self._infile = "inputfile.in"

        engine = engine or SimulationEngine(self)
        self._qe_utils = qeUtils(self, root)
        super().__init__(engine, **kwargs)


    def __str__(self):
        return "Quantum Espresso Wrapper Session"

    def _run(self, input_files, output_files):
        """Run the energy calculation

        Args:
            input_files (dict): should be in format {"-input": "inputfilename"}
            output_files (dict): should be in format {">": "outputfilename"}
        """
        root = self._registry.get(self.root)
        simulation = root.get(oclass = QE.Simulation)

        self._engine.run(input_files = input_files, output_files = output_files)

    def _load_from_backend(self, uids, expired=None):
        for uid in uids:
            try:
                yield self._registry.get(uid)
            except KeyError:
                yield None

    def _update_atom_from_backend(self, atom):
        force = atom.get(oclass = qe.Force)
        position = atom.get(oclass = qe.Position)
        if not position:
            self._update_velocity_from_backend(None, atom)
        if not force:
            self._update_force_from_backend(None, atom)

    def _update_force_from_backend(self, force, cuds_atom = None):
        if cuds_atom is None:
            cuds_atom = self._parent_atom(force)
    
    def _update_position_from_backend(self, position, cuds_atom = None):
        if cuds_atom is None:
            cuds_atom = self._parent_atom(position)

    def _parent_atom(self, cuds_object):
        cuds_atom = cuds_object.get(rel = CUBA.Relationship, oclass = QE.Atom)[0]
        return cuds_atom



    def _apply_added(self, root_obj, buffer):
        if self._ran:
            for cuds_object in buffer.values:
                self._add_by_type(cuds_object)
        else:
            self._add_initial_cuds(buffer)

    def _apply_updated(self, root_obj, buffer):
        for cuds_object in buffer.values():
            self._update_by_type(cuds_object)

    def _apply_deleted(self, root_obj, buffer):
        for cuds_object in ():
            self._remove_by_type(cuds_object)

    
