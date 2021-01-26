from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger

from askcos.synthetic.impurity.impurity_predictor import ImpurityPredictor
from ..atom_mapper.atom_mapping_worker import get_atom_mapping
from ..treebuilder.tfx_fast_filter import TFXFastFilter

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'impurity_worker'
wln_predictor = None


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A IMPURITY PREDICTOR WORKER ###')
    from askcos.synthetic.evaluation.template_free import TemplateFreeNeuralNetScorer
    global wln_predictor

    try:
        wln_predictor = TemplateFreeNeuralNetScorer()
    except Exception as e:
        print(e)
        raise
    print('Initialized')


@shared_task(bind=True)
def get_impurities(self, reactants, reagents='', products='', solvents='',
                   predictor_selection='WLN forward predictor',
                   inspector_selection='WLN forward inspector',
                   mapper_selection='WLN atom mapper',
                   top_k=3, threshold=0.2, check_mapping=True):

    if predictor_selection == 'WLN forward predictor':
        predictor = wln_predictor.evaluate
    else:
        raise NotImplementedError(f'Unsupported predictor: ${predictor_selection}')

    if inspector_selection == 'WLN forward inspector':
        inspector = None
    elif inspector_selection == 'Reaxys inspector':
        fast_filter_hostname = 'fast-filter'
        fast_filter = TFXFastFilter(fast_filter_hostname, 'fast_filter')
        inspector = fast_filter.predict
    else:
        raise NotImplementedError(f'Unsupported predictor: ${predictor_selection}')

    def mapper(*args):
        kwargs = {'mapper': mapper_selection}
        return get_atom_mapping.apply_async(args, kwargs, priority=2).get(30)

    impurity_predictor = ImpurityPredictor(predictor, inspector, mapper,
                                           topn_outcome=top_k, insp_threshold=threshold,
                                           celery_task=self, check_mapping=check_mapping)
    # make prediction
    return impurity_predictor.predict(reactants, reagents=reagents, products=products, solvents=solvents)
