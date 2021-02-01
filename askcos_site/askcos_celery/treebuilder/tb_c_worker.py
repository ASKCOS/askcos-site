"""A worker to build a retrosynthetic tree.

The role of a treebuilder worker is to take a target compound
and apply all retrosynthetic templates to it. The top results
are returned based on the defined mincount and max_branching. The
heuristic chemical scoring function, defined in the transformer
class, is used for prioritization. Each worker pre-loads a
transformer and grabs templates from the database.
"""

import numpy as np
from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger

from askcos.prioritization.precursors.relevanceheuristic import RelevanceHeuristicPrecursorPrioritizer
from askcos_site.globals import scscorer, retro_templates
from .tfx_relevance_template_prioritizer import TFXRelevanceTemplatePrioritizer
from .tfx_fast_filter import TFXFastFilter

relevance_heuristic_prioritizer = RelevanceHeuristicPrecursorPrioritizer()
relevance_heuristic_prioritizer.load_model()

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
CORRESPONDING_QUEUE = 'tb_c_worker'
CORRESPONDING_RESERVABLE_QUEUE = 'tb_c_worker_reservable'
retroTransformer = None

db_comparison_map = {'>': '$gt', '>=': '$gte', '<': '$lt', '<=': '$lte', '==': '$eq'}


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    """Configures worker and instantiates RetroTransformer.

    Args:
        options (dict, optional): Used ensure correct queue. (default: {{}})
        **kwargs: Unused.
    """
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER WORKER ###')

    from askcos_site.askcos_celery.treebuilder.retro_transformer_celery import RetroTransformerCelery

    # Instantiate and load retro transformer
    global retroTransformer
    retroTransformer = RetroTransformerCelery(template_prioritizer=None, fast_filter=None, scscorer=scscorer.get_max_score_from_joined_smiles)
    retroTransformer.load(load_templates=False)
    print('### TREE BUILDER WORKER STARTED UP ###')


@shared_task
def get_top_precursors(
        smiles, precursor_prioritizer=None,
        template_set='reaxys', template_prioritizer_version=None, max_num_templates=1000,
        max_cum_prob=1, fast_filter_threshold=0.75,
        cluster=True, cluster_method='kmeans', cluster_feature='original',
        cluster_fp_type='morgan', cluster_fp_length=512, cluster_fp_radius=1,
        postprocess=False, selec_check=False, attribute_filter=[]
    ):
    """Get the precursors for a chemical defined by its SMILES.

    Args:
        smiles (str): SMILES of node to expand.
        precursor_prioritizer (optional, callable): Use to override
            prioritizer created during initialization. This can be
            any callable function that reorders a list of precursor
            dictionary objects. (default: {None})
        template_set (str): Name of template set to be used during prediction. This 
            determines which documents to use from mongodb. (default: {reaxys})
        template_prioritizer_version (str): Version number for tensorflow serving API. If
            None, the latest model is used. (default: {None})
        max_num_templates (int, optional): Maximum number of templates to consider.
            (default: {1000})
        max_cum_prob (float, optional): Maximum cumulative probability of
            selected relevant templates. (default: {1})
        fast_filter_threshold (float, optional): Threshold to use for fast filter.
            (default: {0.75})
        cluster (bool, optional): Whether to cluster results. (default: {True}). This is passed along to RetroResult.return_top()
        cluster_method (str, optional): Clustering method to use ['kmeans', 'hdbscan']. (default: {'kmeans'})
        cluster_feature (str, optional): Features to use for clustering ['original', 'outcomes', 'all']. 'Original' means features that disappear from original target. 'Outcomes' means new features that appear in predicted precursor outcomes. 'All' means the logical 'or' of both. (default: {'original'})
        cluster_fp_type (str, optional): Type of fingerprint to use. Curretnly only 'morgan' is supported. (default: {'morgan'})
        cluster_fp_length (int, optional): Fixed-length folding to use for fingerprint generation. (default: {512})
        cluster_fp_radius (int, optional): Radius to use for fingerprint generation. (default: {1})
        postprocess (bool): Flag for performing post processing.
        selec_check (bool, optional): apply selectivity checking for precursors to find other outcomes. (default: False)
        attribute_filter (list[dict], optional): template atrtibute filter to apply before template application

    Returns:
        2-tuple of (str, list of dict): SMILES string of input and top
            precursors found.
    """

    template_relevance_hostname = 'template-relevance-{}'.format(template_set)
    template_prioritizer = TFXRelevanceTemplatePrioritizer(
        hostname=template_relevance_hostname, model_name='template_relevance', version=template_prioritizer_version
    )

    fast_filter_hostname = 'fast-filter'
    fast_filter = TFXFastFilter(fast_filter_hostname, 'fast_filter').predict

    cluster_settings = {
        'cluster_method': cluster_method,
        'feature': cluster_feature,
        'fp_type': cluster_fp_type,
        'fp_length': cluster_fp_length,
        'fp_radius': cluster_fp_radius,
    }

    if precursor_prioritizer == 'SCScore':
        precursor_prioritizer = scscorer.reorder_precursors
    else:
        precursor_prioritizer = relevance_heuristic_prioritizer.reorder_precursors

    global retroTransformer
    result = retroTransformer.get_outcomes(
        smiles, template_set=template_set,
        max_num_templates=max_num_templates, max_cum_prob=max_cum_prob, 
        fast_filter_threshold=fast_filter_threshold, template_prioritizer=template_prioritizer,
        precursor_prioritizer=precursor_prioritizer, fast_filter=fast_filter,
        cluster_precursors=cluster, cluster_settings=cluster_settings, selec_check=selec_check,
        attribute_filter=attribute_filter
    )
    
    if postprocess:
        for r in result:
            r['templates'] = r.pop('tforms')
        return result
    else:
        return smiles, result

@shared_task
def template_relevance(smiles, max_num_templates, max_cum_prob,
                       template_set='reaxys', template_prioritizer_version=None,
                       return_templates=False, attribute_filter=None):
    """
    Celery task for template_relevance prediction.

    If return_templates = True, return as list of dictionaries containing
    template information.
    """
    global retroTransformer

    hostname = 'template-relevance-{}'.format(template_set)
    template_prioritizer = TFXRelevanceTemplatePrioritizer(
        hostname=hostname,
        model_name='template_relevance',
        version=template_prioritizer_version,
    )

    scores, indices = template_prioritizer.predict(smiles, max_num_templates=None, max_cum_prob=None)

    # Filter results by template attributes if needed
    if attribute_filter:
        query = {
            'index': {'$in': indices.tolist()},
            'template_set': template_set
        }
        for item in attribute_filter:
            query['attributes.{}'.format(item['name'])] = {
                db_comparison_map[item['logic']]: item['value']
            }

        result_proj = ['_id', 'index', 'reaction_smarts'] if return_templates else ['index']
        cursor = retro_templates.find(query, result_proj)

        template_map = {x['index']: x for x in cursor}
        bool_mask = np.array([i in template_map for i in indices.tolist()])
        scores = scores[bool_mask]
        indices = indices[bool_mask]
    else:
        template_map = None

    # Then trim results based on max_num_templates and max_cum_prob
    if max_num_templates is not None:
        scores = scores[:max_num_templates]
        indices = indices[:max_num_templates]

    if max_cum_prob is not None:
        bool_mask = np.cumsum(scores) <= max_cum_prob
        scores = scores[bool_mask]
        indices = indices[bool_mask]

    if not isinstance(scores, list):
        scores = scores.tolist()
    if not isinstance(indices, list):
        indices = indices.tolist()

    if return_templates:
        if template_map is None:
            cursor = retro_templates.find({
                    'index': {'$in': indices},
                    'template_set': template_set
                }, ['_id', 'index', 'reaction_smarts'])

            template_map = {x['index']: x for x in cursor}

        templates = [template_map[i] for i in indices]

        # Add score to template document
        for i, (score, template) in enumerate(zip(scores, templates)):
            template['rank'] = i + 1
            template['score'] = score

        return templates
    else:
        return scores, indices

@shared_task
def apply_one_template_by_idx(*args, **kwargs):
    """Wrapper function for ``RetroTransformer.apply_one_template_by_idx``.

    Returns:
        list of 5-tuples of (int, str, int, list, float): Result of
            applying given template to the molecule.
    """
    global retroTransformer

    postprocess = kwargs.pop('postprocess', False)

    template_set = kwargs.get('template_set', 'reaxys')
    template_prioritizer_version = kwargs.pop('template_prioritizer_version', None)

    hostname = 'template-relevance-{}'.format(template_set)
    template_prioritizer = TFXRelevanceTemplatePrioritizer(
        hostname=hostname, model_name='template_relevance', version=template_prioritizer_version
    )

    fast_filter_hostname = 'fast-filter'
    fast_filter = TFXFastFilter(fast_filter_hostname, 'fast_filter').predict

    kwargs.update({
        'template_prioritizer': template_prioritizer,
        'fast_filter': fast_filter
    })

    result = retroTransformer.apply_one_template_by_idx(*args, **kwargs)

    if postprocess:
        # Unpack result into list of dicts
        # First item is the tree builder path ID which is not needed
        result = [{
            'smiles': r[1],
            'template_idx': r[2],
            'precursors': r[3],
            'ffscore': r[4],
        } for r in result]

    return result

@shared_task
def fast_filter_check(*args, **kwargs):
    """Wrapper for fast filter check.

    These workers will already have it initialized. Best way to allow
    independent queries.

    Returns:
        list: Reaction outcomes.
    """
    print('got request for fast filter')
    fast_filter_hostname = 'fast-filter'
    fast_filter = TFXFastFilter(fast_filter_hostname, 'fast_filter')
    return fast_filter.predict(*args, **kwargs)


@shared_task(bind=True)
def reserve_worker_pool(self):
    """Reserves pool of workers.

    Called by a tb_coordinator to reserve this pool of workers to do a tree
    expansion. This is accomplished by changing what queue(s) this pool
    listens to.

    Returns:
        str: Name of the new queue the worker pool listens to.
    """
    hostname = self.request.hostname
    private_queue = CORRESPONDING_QUEUE + '_' + hostname
    print('Tried to reserve this worker!')
    print('I am {}'.format(hostname))
    print('Telling myself to ignore the {} and {} queues'.format(
        CORRESPONDING_QUEUE, CORRESPONDING_RESERVABLE_QUEUE))
    from askcos_site.celery import app
    app.control.cancel_consumer(CORRESPONDING_QUEUE, destination=[hostname])
    app.control.cancel_consumer(
        CORRESPONDING_RESERVABLE_QUEUE, destination=[hostname])

    # *** purge the queue in case old jobs remain
    import celery.bin.amqp
    amqp = celery.bin.amqp.amqp(app=app)
    amqp.run('queue.purge', private_queue)
    print('Telling myself to only listen to the new {} queue'.format(private_queue))
    app.control.add_consumer(private_queue, destination=[hostname])
    return private_queue


@shared_task(bind=True)
def unreserve_worker_pool(self):
    """Releases this worker pool so it can listen to the original queues.

    Returns:
        True
    """
    hostname = self.request.hostname
    private_queue = CORRESPONDING_QUEUE + '_' + hostname
    print('Tried to unreserve this worker!')
    print('I am {}'.format(hostname))
    print('Telling myself to ignore the {} queue'.format(private_queue))
    from askcos_site.celery import app
    app.control.cancel_consumer(private_queue, destination=[hostname])
    print('Telling myself to only listen to the {} and {} queues'.format(
        CORRESPONDING_QUEUE, CORRESPONDING_RESERVABLE_QUEUE))
    app.control.add_consumer(CORRESPONDING_QUEUE, destination=[hostname])
    app.control.add_consumer(
        CORRESPONDING_RESERVABLE_QUEUE, destination=[hostname])
    return True
