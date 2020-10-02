'''
descritptors predictor
'''
from __future__ import absolute_import, unicode_literals, print_function
from celery import shared_task
from celery.signals import celeryd_init
from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

descriptor_pred = None
CORRESPONDING_QUEUE = 'descriptors_worker'

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A DESCRIPTOR PREDICTOR WORKER ###')
    global desc_pred
    # Import as needed
    from ..torchserve_api import TorchserveAPI
    try:
        desc_pred = TorchserveAPI(hostname='descriptors', model_name='descriptors')
    except Exception as e:
        raise(e)
    print('Initialized')
    print('Finished configuring descriptor worker')


@shared_task
def get_descriptors(smi):
    global desc_pred
    print('descriptor predictor got a request {}'.format(smi))
    desc_pred = desc_pred.predict(smi.split('.'))
    return desc_pred