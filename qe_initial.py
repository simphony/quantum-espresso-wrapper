from osp.core import Parser
p = Parser()
p.parse("ontology.qe.yml")
from osp.core import cuba
from osp.core import QE_ONTOLOGY as qe

from qe_syntactic_layer import run_pw

run_pw('basicsuperminusone.in')
