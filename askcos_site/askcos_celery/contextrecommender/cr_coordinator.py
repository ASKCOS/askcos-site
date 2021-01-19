from celery import shared_task
from celery.result import allow_join_result
from celery.signals import celeryd_init

import askcos.global_config as gc
from askcos_site.askcos_celery.contextrecommender.cr_network_worker import get_n_conditions as network_get_n_conditions
from askcos_site.askcos_celery.contextrecommender.cr_nn_worker import get_n_conditions as neighbor_get_n_conditions
from askcos_site.askcos_celery.contextrecommender.cr_network_v2_worker import get_n_conditions as network_v2_get_n_conditions

CORRESPONDING_QUEUE = 'cr_coordinator'


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A CONTEXT RECOMMENDER COORDINATOR ###')
    print('### CONTEXT RECOMMENDER COORDINATOR STARTED UP ###')


@shared_task
def get_context_recommendations(*args, **kwargs):
    """Retrieve a context recommendation given the reaction to attempt.

    rxn = [reacants, products], Where each is a list of SMILES.
    n = Number of contexts to return.
    """

    context_recommender = kwargs.pop('context_recommender', gc.nearest_neighbor)

    print('Context context_recommender worker got a request: {}, {}'.format(args, kwargs))
    with allow_join_result():
        if context_recommender == gc.nearest_neighbor:
            res = neighbor_get_n_conditions.delay(*args, **kwargs)
            return res.get()
        elif context_recommender == gc.neural_network:
            res = network_get_n_conditions.delay(*args, **kwargs)
            return res.get()
        elif context_recommender == gc.context_neural_network_v2:
            res = network_v2_get_n_conditions.delay(*args, **kwargs)
            return res.get()
        else:
            raise NotImplementedError

@shared_task
def get_recommender_types():
    return [gc.nearest_neighbor,gc.neural_network,gc.context_neural_network_v2]
