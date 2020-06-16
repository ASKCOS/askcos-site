import os
from celery import Celery

# Set the Django settings module
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'askcos_site.settings')

# Define readable names for celery workers for status reporting
READABLE_NAMES = {
    'cr_network_worker': 'Context Recommender Worker',
    'tb_c_worker': 'One-Step/Tree Builder Retrosynthesis Worker',
    'tb_coordinator_mcts': 'Tree Builder Coordinator',
    'tb_coordinator_mcts_v2': 'Tree Builder v2 Coordinator',
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
app = Celery(
    main='askcos_site',
    broker='amqp://{}:{}'.format(RABBIT_HOST, RABBIT_PORT),
    backend='redis://{}:{}'.format(REDIS_HOST, REDIS_PORT),
)

# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
app.config_from_object('askcos_site.askcos_celery.celeryconfig')


if __name__ == '__main__':
    app.start()
