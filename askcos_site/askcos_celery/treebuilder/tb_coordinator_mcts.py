"""
The role of a treebuilder coordinator is to take a target compound
and build up the retrosynthetic tree by sending individual chemicals
to workers, which each apply the full set of templates. The coordinator
will keep track of the dictionary and ensure unique IDs in addition
to keeping track of the chemical prices using the Pricer module.

The coordinator, finally, returns a set of buyable trees obtained
from an IDDFS.
"""

from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger

from askcos_site.globals import db_client
from askcos_site.main.models import SavedResults

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'tb_coordinator_mcts'

results_collection = db_client['results']['results']


def update_result_state(id_, state):
    result = SavedResults.objects.get(result_id=id_)
    result.result_state = state
    result.save()
    return


def save_results(result, settings, task_id):
    doc = {
        '_id': task_id,
        'result': result,
        'settings': settings
    }
    results_collection.insert_one(doc)


@celeryd_init.connect
def configure_coordinator(options={}, **kwargs):
    """Initializes coordinator for MCTS tree building.

    Args:
        options (dict, optional): Used to check if the queue is correct.
            (default: {{}})
        **kwargs: Unused.
    """
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TREE BUILDER MCTS COORDINATOR ###')

    from .tree_builder_celery import MCTSCelery

    global treeBuilder

    treeBuilder = MCTSCelery(celery=True, nproc=8)  # 8 active pathways
    print('Finished initializing treebuilder MCTS coordinator')


@shared_task(trail=False)
def get_buyable_paths(*args, **kwargs):
    """Wrapper for ``MCTSTreeBuilder.get_buyable_paths`` function.

    Returns:
        tree_status ((int, int, dict)): Result of tree_status().
        trees (list of dict): List of dictionaries, where each dictionary
            defines a synthetic route.
    """
    run_async = kwargs.pop('run_async', False)
    paths_only = kwargs.pop('paths_only', False)

    settings = {'smiles': args[0], 'version': 1}  # Refers to tree builder version
    settings.update(kwargs)

    kwargs['json_format'] = 'nodelink'

    template_prioritizer_version = kwargs.pop('template_prioritizer_version', None)
    if template_prioritizer_version:
        treeBuilder.template_prioritizer_version = template_prioritizer_version


    print('Treebuilder MCTS coordinator was asked to expand {}'.format(args[0]))
    _id = get_buyable_paths.request.id
    try:
        paths, status, graph = treeBuilder.get_buyable_paths(*args, **kwargs)
        result_doc = {
            'status': status,
            'paths': paths,
            'graph': graph,
            'version': 2,  # Refers to graph version
        }
    except:
        if run_async:
            update_result_state(_id, 'failed')
        raise
    if run_async:
        update_result_state(_id, 'completed')
        save_results(result_doc, settings, _id)
    print('Task completed, returning results.')

    if paths_only:
        return paths
    else:
        return status, paths
