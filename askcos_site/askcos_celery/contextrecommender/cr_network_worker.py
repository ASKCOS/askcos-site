"""
The role of a context_worker is to take in an attempted reaction and return
a set of conditions to (try to) run the reaction in. Each worker will
load a pre-trained neural network model. For each request, this worker
must query the database to get details about the instance.
"""

from celery import shared_task
from celery.signals import celeryd_init

CORRESPONDING_QUEUE = 'cr_network_worker'


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A NEURAL NETWORK CONTEXT RECOMMENDER WORKER ###')

    global recommender

    # Setting logging low
    from rdkit import RDLogger
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    from askcos.synthetic.context.neuralnetwork import NeuralNetContextRecommender
    try:
        recommender = NeuralNetContextRecommender()
        recommender.load()
    except Exception as e:
        print(e)
    print('Loaded context recommendation model')

    print('### NEURAL NETWORK CONTEXT RECOMMENDER STARTED UP ###')


@shared_task
def get_n_conditions(*args, **kwargs):
    """Retrieve a context recommendation given the reaction to attempt.

    rxn = [reacants, products], Where each is a list of SMILES.
    n = Number of contexts to return.
    """
    postprocess = kwargs.pop('postprocess', False)

    print('Context recommender worker got a request: {}, {}'.format(args, kwargs))
    res = recommender.get_n_conditions(*args, **kwargs)

    if postprocess:
        if kwargs.get('return_scores'):
            contexts, scores = res
        else:
            contexts, scores = res, None

        res = []
        for context in contexts:
            res.append({
                'temperature': context[0],
                'solvent': context[1],
                'reagent': context[2],
                'catalyst': context[3]
            })

        if scores is not None:
            for c, s in zip(res, scores):
                c['score'] = s

    print('Task completed, returning results.')
    return res
