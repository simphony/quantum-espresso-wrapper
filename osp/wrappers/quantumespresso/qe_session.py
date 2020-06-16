from osp.core import Parser
p = Parser()
p.parse("ontology.qe.yml")

from osp.core.session import SimWrapperSession
from qe_engine import qeEngine
from osp.core import QE


class QESession(SimWrapperSession):
    

    def __init__(self, engine, **kwargs):
        super().__init__(engine, **kwargs)

    def __str__(self):
        return "Quantum Espresso Wrapper Session"

    def _apply_added(self, root_obj, buffer):
        return super()._apply_added(root_obj, buffer)
        for obj in buffer.values():
            if obj.is_a():

                self._engine.add_atom()

    def _apply_updated(self, root_obj, buffer):
        for obj in buffer.values():
            if obj.is_a(QE.CELL):
                self._engine.

    def _apply_deleted(self, root_obj, buffer):
        return super()._apply_deleted(root_obj, buffer)

    def _load_from_backend(self, uids, expired=None):
        return super()._load_from_backend(uids, expired=expired)
    
    def _run(self, root_cuds_object):
        return super()._run(root_cuds_object)


    
