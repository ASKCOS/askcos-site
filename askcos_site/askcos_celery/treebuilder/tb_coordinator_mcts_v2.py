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

from askcos_site.globals import db_client, pricer, scscorer, chemical_db
from askcos_site.main.models import SavedResults
from askcos.retrosynthetic.mcts.v2.tree_builder import MCTS
from askcos.utilities.historian.chemicals import ChemHistorian
from .tfx_relevance_template_prioritizer import TFXRelevanceTemplatePrioritizer
from .tfx_fast_filter import TFXFastFilter

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

CORRESPONDING_QUEUE = 'tb_coordinator_mcts_v2'

results_collection = db_client['results']['results']

retroTransformer = None
chemhistorian = None


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

    from askcos_site.askcos_celery.treebuilder.retro_transformer_celery import RetroTransformerCelery

    global retro_transformer
    retro_transformer = RetroTransformerCelery(template_prioritizer=None, fast_filter=None, scscorer=None)
    retro_transformer.load()

    global chemhistorian
    chemhistorian = ChemHistorian(CHEMICALS_DB=chemical_db, use_db=True)

    print('Finished initializing treebuilder MCTS coordinator')


@shared_task(trail=False)
def get_buyable_paths(smiles, **kwargs):
    """Wrapper for ``MCTSTreeBuilder.get_buyable_paths`` function.

    Returns:
        tree_status ((int, int, dict)): Result of tree_status().
        trees (list of dict): List of dictionaries, where each dictionary
            defines a synthetic route.
    """
    run_async = kwargs.pop('run_async', False)
    paths_only = kwargs.pop('paths_only', False)

    print('Treebuilder MCTS coordinator was asked to expand {}'.format(smiles))
    _id = get_buyable_paths.request.id

    template_prioritizer = kwargs.pop('template_prioritizer', 'reaxys')

    hostname = 'template-relevance-{}'.format(template_prioritizer)
    template_prioritizer = TFXRelevanceTemplatePrioritizer(
        hostname=hostname, model_name='template_relevance'
    )

    fast_filter_hostname = 'fast-filter'
    fast_filter = TFXFastFilter(fast_filter_hostname, 'fast_filter').predict

    global retro_transformer
    global chemhistorian

    kwargs.update({
        'pricer': pricer,
        'scscorer': scscorer,
        'chemhistorian': chemhistorian,
        'retro_transformer': retro_transformer,
        'template_prioritizer': template_prioritizer,
        'fast_filter': fast_filter,
    })

    tree_builder = MCTS(**kwargs)

    try:
        paths = tree_builder.get_buyable_paths(smiles, **kwargs)
        status = (len(tree_builder.chemicals), len(tree_builder.reactions))
        graph = tree_builder.dump_tree()
        result_doc = {
            'status': status,
            'paths': paths,
            'graph': graph,
            'version': 2,
        }
    except:
        if run_async:
            update_result_state(_id, 'failed')
        raise
    if run_async:
        update_result_state(_id, 'completed')
        settings = {'smiles': smiles}
        settings.update(kwargs)
        save_results(result_doc, settings, _id)
    print('Task completed, returning results.')

    if paths_only:
        return paths
    else:
        return status, paths
