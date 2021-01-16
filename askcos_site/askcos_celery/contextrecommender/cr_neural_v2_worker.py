"""
The role of a context_worker is to take in an attempted reaction and return
a set of conditions to (try to) run the reaction in. Each worker will
load a pre-trained neural network model. For each request, this worker
must query the database to get details about the instance.
"""

from celery import shared_task
from celery.signals import celeryd_init

import askcos.global_config as gc
from askcos.synthetic.context.v2.reaction_context_predictor import ReactionContextRecommenderWLN, ReactionContextRecommenderFP

from . import lru_cache

CORRESPONDING_QUEUE = 'cr_neural_v2_worker'
MODEL_CACHE_CAPACITY = 2

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A NEURAL NETWORK CONTEXT V2 RECOMMENDER WORKER ###')

    global model_cache
    model_cache = lru_cache.LRUCache(capacity=MODEL_CACHE_CAPACITY)

    # Setting logging low
    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    
    print('Context recommendation v2 model will be loaded on-demand.')
    print('MODEL_CACHE_CAPACITY={}'.format(MODEL_CACHE_CAPACITY))

    print('### NEURAL NETWORK CONTEXT V2 RECOMMENDER STARTED UP ###')


@shared_task
def get_n_conditions(*args, **kwargs):
    """Retrieve a context recommendation given the reaction to attempt.

    reaction (str): Where each is a list of SMILES.
    n = Number of contexts to return. It is also the beam width.
    """
    model_name = kwargs.pop('model_name', 'graph')
    if model_name in gc.CONTEXT_V2['default-models']:
        model_name = gc.CONTEXT_V2['default-models'][model_name]

    print('Context recommender v2 worker got a request: {}, {}'.format(args, kwargs))
    try:
        model = model_cache.get(model_name)
    except KeyError:
        print('Context recommender v2 worker is loading new model: {}'.format(model_name))
        model_config = gc.CONTEXT_V2['models'].get(model_name, None)
        
        if model_config is None:
            print('Context recommender v2 worker cannot load model: model_config does not exist')
            return None
        
        if model_name.startswith('graph'):
            model = ReactionContextRecommenderWLN(**model_config)
        elif model_name.startswith('fp'):
            model = ReactionContextRecommenderFP(**model_config)
        else:
            print('Context recommender v2 worker cannot determine model type')
            return None
        
        model_cache.set(model_name, model)
    
    res = model.predict(*args, **kwargs)

    print('Task completed, returning results.')
    return res
