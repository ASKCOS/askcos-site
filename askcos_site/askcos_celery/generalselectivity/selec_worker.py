"""
selec selectivity predictor
"""

from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

selec_pred = None
CORRESPONDING_QUEUE = 'selec_worker'


@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A GENERAL SELECTIVITY PREDICTOR WORKER ###')
    global selec_pred
    # Import as needed
    from askcos.synthetic.selectivity.general_selectivity import GeneralSelectivityPredictor
    try:
        selec_pred = GeneralSelectivityPredictor()
    except Exception as e:
        raise (e)
    print('Initialized')
    configure_reac = '[Br:1][Br:2].[NH2:3][c:4]1[n:5][cH:6][n:7][c:8]2[nH:9][cH:10][n:11][c:12]12>O>' \
                     '[Br:2][c:10]1[nH:9][c:8]2[n:7][cH:6][n:5][c:4]([NH2:3])[c:12]2[n:11]1.' \
                     '[Br:2][c:6]1[n:5][c:4]([NH2:3])[c:12]2[c:8]([n:7]1)[nH:9][cH:10][n:11]2'
    print(selec_pred.predict(configure_reac))
    print('Finished configuring sites worker')


@shared_task
def get_selec(reac):
    global selec_pred
    print('site selectivity got a request {}'.format(reac))
    res = selec_pred.predict(reac)
    return res
