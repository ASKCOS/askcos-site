from __future__ import absolute_import, unicode_literals, print_function
import os
from celery import Celery
from django.conf import settings

# Set the Django settings module
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'askcos_site.settings')

# Define readable names for celery workers for status reporting
READABLE_NAMES = {
    'cr_network_worker': 'Context Recommender Worker',
    'tb_c_worker': 'One-Step/Tree Builder Retrosynthesis Worker',
    'tb_coordinator_mcts': 'Tree Builder Coordinator',
    'sites_worker': 'Site Selectivity Worker',
    'impurity_worker': 'Impurity Worker',
    'atom_mapping_worker': 'Atom Mapping Worker',
    'tffp_worker': 'Template-free Forward Predictor',
    'selec_worker': 'General Selectivity Worker',
}

# Note: cannot use guest for authenticating with broker unless on localhost
REDIS_HOST = os.environ.get('REDIS_HOST', 'localhost')
REDIS_PORT = os.environ.get('REDIS_PORT', '6379')
RABBIT_HOST = os.environ.get('RABBIT_HOST', 'localhost')
RABBIT_PORT = os.environ.get('RABBITMQ_NODE_PORT', '5672')
app = Celery('askcos_site', broker='amqp://{}:{}'.format(RABBIT_HOST, RABBIT_PORT),
    backend='redis://{}:{}'.format(REDIS_HOST, REDIS_PORT),
    include=[
        'askcos_site.askcos_celery.treebuilder.tb_c_worker',
        'askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts',
        'askcos_site.askcos_celery.contextrecommender.cr_network_worker',
        'askcos_site.askcos_celery.treeevaluator.template_free_forward_predictor_worker',
        'askcos_site.askcos_celery.siteselectivity.sites_worker',
        'askcos_site.askcos_celery.impurity.impurity_worker',
        'askcos_site.askcos_celery.impurity.impurity_predictor_worker',
        'askcos_site.askcos_celery.atom_mapper.atom_mapping_worker',
        'askcos_site.askcos_celery.generalselectivity.selec_worker',
    ]
)

# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')
app.conf.task_queue_max_priority = 20 # necessary for new tb_worker queues to be priority

if __name__ == '__main__':
    app.start()
