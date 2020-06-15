from osp.core import Parser
p = Parser()
p.parse("ontology.qe.yml")
from osp.core import cuba
from osp.core import QE_ONTOLOGY as qe

from osp.core.session import SimWrapperSession

class QESession(SimWrapperSession):
    

    def __init__(self):
        