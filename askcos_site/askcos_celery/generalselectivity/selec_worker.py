"""
selec selectivity predictor
"""

from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger

from ..atom_mapper.atom_mapping_worker import get_atom_mapping
from ..descriptors.descriptors_worker import get_descriptors
from askcos.synthetic.selectivity.general_selectivity import GeneralSelectivityPredictor

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'selec_worker'


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A GENERAL SELECTIVITY PREDICTOR WORKER ###')
    print('Initialized')


@shared_task
def get_selec(reac,  mapped=False, mode='qm-gnn', outcomes_included=False):
    print('site selectivity got a request {}'.format(reac))

    def mapper(rxnsmiles):
        result = get_atom_mapping.delay(rxnsmiles)
        return result.get(10)

    def descriptors(smis):
        result = get_descriptors.delay(smis)
        return result.get(5)

    selec_pred = GeneralSelectivityPredictor(atom_mapper=mapper, descriptor_predictor=descriptors)
    res = selec_pred.predict(reac, mapped=mapped, mode=mode, outcomes_included=outcomes_included)
    return res
