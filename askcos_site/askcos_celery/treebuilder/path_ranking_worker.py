"""
Celery worker for performing pathway ranking and clustering.
"""

from celery import shared_task
from celery.signals import celeryd_init

from askcos.retrosynthetic.pathway_ranker import PathwayRanker
from ..torchserve_api import TorchserveAPI

CORRESPONDING_QUEUE = 'path_ranking_worker'


class TSPathwayRanker(TorchserveAPI, PathwayRanker):
    """
    Version of pathway ranker using torchserve to perform inference.

    Uses ``TorchserveAPI.predict`` to override ``PathwayRanker.predict``.
    """
    pass


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    """Initializes coordinator for pathway ranking.

    Args:
        options (dict, optional): Used to check if the queue is correct.
            (default: {{}})
        **kwargs: Unused.
    """
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP PATH RANKING WORKER ###')
    print('### PATH RANKING WORKER READY ###')


@shared_task
def get_scores(trees, **kwargs):
    """
    Wrapper for ``PathwayRanker.scorer`` function.

    Returns:
        tree_status ((int, int, dict)): Result of tree_status().
        trees (list of dict): List of dictionaries, where each dictionary
            defines a synthetic route.
    """
    ranker = TSPathwayRanker(hostname='ts-pathway-ranker', model_name='pathway-ranker')
    output = ranker.scorer(trees, **kwargs)
    return output
