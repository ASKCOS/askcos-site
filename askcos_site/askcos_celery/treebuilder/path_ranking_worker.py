"""
Celery worker for performing pathway ranking and clustering.
"""

import scipy.sparse
from celery import shared_task
from celery.signals import celeryd_init

from askcos.retrosynthetic.pathway_ranker import PathwayRanker
from askcos.retrosynthetic.pathway_ranker.utils import convert_askcos_trees, tree_to_input, merge_into_batch
from ..torchserve_api import TorchserveAPI

CORRESPONDING_QUEUE = 'path_ranking_worker'


class TSPathwayRanker(TorchserveAPI, PathwayRanker):
    """
    Version of pathway ranker using torchserve to perform inference.
    """

    def preprocess(self, trees):
        """
        Convert trees into fingerprints for model input. Fingerprint arrays
        are converted to sparse arrays in Compressed Sparse Column format
        to reduce request size when querying the torchserve API.

        One-step trees are removed from the input because they cannot be ranked.

        The original indices of the remaining trees are also returned to help
        map the results to the original list of trees.

        Args:
            trees (list): list of askcos trees to process

        Returns:
            original_indices (list): list of the original indices of the selected trees
            batch (dict): dictionary of generated fingerprints and metadata
        """
        def compress(mat):
            csc_mat = scipy.sparse.csc_matrix(mat)
            return csc_mat.data.tolist(), csc_mat.indices.tolist(), csc_mat.indptr.tolist()

        output = convert_askcos_trees(trees)

        # Separate out one-step trees because they can't be ranked
        original_indices, remaining_trees = zip(*((i, tree) for i, tree in enumerate(output) if tree['depth'] > 1))

        # Convert trees to set of input fingerprints and metadata
        batch = merge_into_batch([tree_to_input(tree) for tree in remaining_trees])

        # Convert fingerprint arrays into tensors
        batch['pfp'] = compress(batch['pfp'])
        batch['rxnfp'] = compress(batch['rxnfp'])
        batch['node_order'] = batch['node_order'].tolist()
        batch['adjacency_list'] = batch['adjacency_list'].tolist()
        batch['edge_order'] = batch['edge_order'].tolist()

        return list(original_indices), batch

    def postprocess(self, data):
        return data


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
