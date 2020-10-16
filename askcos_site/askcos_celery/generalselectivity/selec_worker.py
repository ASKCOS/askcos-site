"""
selec selectivity predictor
"""

from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger

from ..atom_mapper.atom_mapping_worker import get_atom_mapping
from ..descriptors.descriptors_worker import get_descriptors
from askcos.synthetic.selectivity.general_selectivity import QmGnnGeneralSelectivityPredictor

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'selec_worker'
wln_atom_mapper = None
transformer_mapper = None
descriptor_predictor = None
selec_pred = None

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A GENERAL SELECTIVITY PREDICTOR WORKER ###')

    from askcos_site.askcos_celery.torchserve_api import TorchserveAPI
    from askcos.synthetic.atom_mapper.wln_mapper import WLN_AtomMapper
    global wln_atom_mapper
    global transformer_mapper
    global descriptor_predictor
    global selec_pred

    wln_atom_mapper = WLN_AtomMapper()
    transformer_mapper = TorchserveAPI(hostname='ts-rxnmapper', model_name='rxnmapper')
    descriptor_predictor = TorchserveAPI(hostname='ts-descriptors', model_name='descriptors')
    mapper_func = wln_atom_mapper.evaluate

    def descriptors_predictor_wrapper(smiles):
        return descriptor_predictor.predict(smiles.split('.'))
    descriptors = descriptors_predictor_wrapper

    selec_pred = QmGnnGeneralSelectivityPredictor(atom_mapper=mapper_func, descriptor_predictor=descriptors)

    print('Initialized')


@shared_task
def get_selec(reac,  mapped=False, mode='qm-gnn', all_outcomes=False, verbose=True, mapper='Transformer', no_map_reagents=False):
    print('site selectivity got a request {}'.format(reac))

    global wln_atom_mapper
    global transformer_mapper
    global descriptor_predictor
    global selec_pred

    if mapper == 'WLN atom mapper':
        mapper_func = wln_atom_mapper.evaluate
    else:
        def transformer_mapper_wrapper(rxnsmiles):
            return transformer_mapper.predict([rxnsmiles])[0]['mapped_rxn']
        mapper_func = transformer_mapper_wrapper

    def descriptors_predictor_wrapper(smiles):
        return descriptor_predictor.predict(smiles.split('.'))
    descriptors = descriptors_predictor_wrapper


    #try:
    res = selec_pred.predict(reac, mapped=mapped, mode=mode, all_outcomes=all_outcomes, verbose=verbose, no_map_reagents=no_map_reagents)
    #except Exception as e:
    #    return str(e)
    #else:
    return res
