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
selec_pred_gnn = None
selec_pred_qm_gnn = None
selec_pred_qm_gnn_no_reagent = None

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
    from ..atom_mapper.atom_mapping_worker import get_atom_mapping
    from ..descriptors.descriptors_worker import get_descriptors
    from askcos.synthetic.selectivity.general_selectivity import GnnGeneralSelectivityPredictor, \
        QmGnnGeneralSelectivityPredictor, QmGnnGeneralSelectivityPredictorNoReagent

    global wln_atom_mapper
    global transformer_mapper
    global descriptor_predictor
    global selec_pred_gnn
    global selec_pred_qm_gnn
    global selec_pred_qm_gnn_no_reagent

    wln_atom_mapper = WLN_AtomMapper()
    transformer_mapper = TorchserveAPI(hostname='ts-rxnmapper', model_name='rxnmapper')
    descriptor_predictor = TorchserveAPI(hostname='ts-descriptors', model_name='descriptors')

    selec_pred_gnn = GnnGeneralSelectivityPredictor(atom_mapper='mapper', descriptor_predictor='descriptor')
    selec_pred_qm_gnn = QmGnnGeneralSelectivityPredictor(atom_mapper='mapper', descriptor_predictor='descriptor')
    selec_pred_qm_gnn_no_reagent = QmGnnGeneralSelectivityPredictorNoReagent(atom_mapper='mapper', descriptor_predictor='descriptor')

    print('Initialized')


@shared_task
def get_selec(reac,  mapped=False, mode='qm_GNN', all_outcomes=False, verbose=True, mapper='Transformer', no_map_reagents=False):
    print('site selectivity got a request {}'.format(reac))

    global wln_atom_mapper
    global transformer_mapper
    global descriptor_predictor
    global selec_pred_gnn
    global selec_pred_qm_gnn
    global selec_pred_qm_gnn_no_reagent

    if mapper == 'WLN atom mapper':
        mapper_func = wln_atom_mapper.evaluate
    else:
        def transformer_mapper_wrapper(rxnsmiles):
            return transformer_mapper.predict([rxnsmiles])[0]['mapped_rxn']
        mapper_func = transformer_mapper_wrapper

    def descriptors_predictor_wrapper(smiles):
        return descriptor_predictor.predict(smiles.split('.'))
    descriptors = descriptors_predictor_wrapper

    try:
        if mode == 'GNN':
            res = selec_pred_gnn.predict(reac, mapper_func, descriptors, mapped=mapped, all_outcomes=all_outcomes,
                                         verbose=verbose, no_map_reagents=no_map_reagents)
        elif mode == 'qm_GNN':
            _, reagent, _ = reac.split('>')

            if reagent:
                res = selec_pred_qm_gnn.predict(reac, mapper_func, descriptors, mapped=mapped, all_outcomes=all_outcomes,
                                                   verbose=verbose, no_map_reagents=no_map_reagents)
            else:
                res = selec_pred_qm_gnn_no_reagent.predict(reac, mapper_func, descriptors, mapped=mapped, all_outcomes=all_outcomes,
                                        verbose=verbose, no_map_reagents=no_map_reagents)
    except Exception as e:
        return str(e)
    else:
        return res
