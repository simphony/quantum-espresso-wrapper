from osp.core import Parser
p = Parser()
p.parse("ontology.qe.yml")
from osp.core import QE
from osp.core.utils import pretty_print
from osp.wrappers.quantumespresso import qe_session
