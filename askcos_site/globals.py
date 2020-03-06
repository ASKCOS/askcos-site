"""
Data and model instances for global use in ``askcos_site``.
"""

from pymongo import MongoClient
# Setting logging low
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

import makeit.global_config as gc
from makeit.prioritization.precursors.scscore import SCScorePrecursorPrioritizer
from makeit.retrosynthetic.transformer import RetroTransformer
from makeit.utilities.buyable.pricer import Pricer


################################################################################
# Database client
db_client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect=gc.MONGO['connect'])

################################################################################
# Database collections
reaction_db = db_client[gc.REACTIONS['database']][gc.REACTIONS['collection']]

chemical_db = db_client[gc.CHEMICALS['database']][gc.CHEMICALS['collection']]

buyables_db = db_client[gc.BUYABLES['database']][gc.BUYABLES['collection']]

solvent_db = db_client[gc.SOLVENTS['database']][gc.SOLVENTS['collection']]

retro_templates = db_client[gc.RETRO_TEMPLATES['database']][gc.RETRO_TEMPLATES['collection']]

forward_templates = db_client[gc.FORWARD_TEMPLATES['database']][gc.FORWARD_TEMPLATES['collection']]

################################################################################
# Retro Transformer
retro_transformer = RetroTransformer(template_prioritizer=None, precursor_prioritizer=None, fast_filter=None)
retro_transformer.load()
RETRO_CHIRAL_FOOTNOTE = 'Using {} chiral retrosynthesis templates from {}/{}'.format(
    gc.RELEVANCE_TEMPLATE_PRIORITIZATION['reaxys']['output_size'],
    gc.RETRO_TEMPLATES['database'],
    gc.RETRO_TEMPLATES['collection']
)

################################################################################
# Pricer
pricer = Pricer()
pricer.load()
print('Loaded global pricer.')

################################################################################
# SCScorer
scscorer = SCScorePrecursorPrioritizer()
scscorer.load_model(model_tag='1024bool')
print('Loaded global SCScorer.')
