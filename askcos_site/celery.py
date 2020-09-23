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
    'path_ranking_worker': 'Pathway Ranking Worker',
    'descriptors_worker': 'Descriptor predictor worker'
}

# Note: cannot use guest for authenticating with broker unless on localhost
redis_host = os.environ.get('REDIS_HOST', 'localhost')
redis_port = os.environ.get('REDIS_PORT', '6379')
redis_password = os.environ.get('REDIS_PASSWORD', '')
if redis_password:
    redis_password = ':{0}@'.format(redis_password)
redis_url = 'redis://{password}{host}:{port}'.format(
    password=redis_password,
    host=redis_host,
    port=redis_port,
)

rabbit_host = os.environ.get('RABBIT_HOST', 'localhost')
rabbit_port = os.environ.get('RABBITMQ_NODE_PORT', '5672')
rabbit_user = os.environ.get('RABBIT_USER', 'guest')
rabbit_password = os.environ.get('RABBIT_PASSWORD', 'guest')
rabbit_url = 'amqp://{user}:{password}@{host}:{port}'.format(
    user=rabbit_user,
    password=rabbit_password,
    host=rabbit_host,
    port=rabbit_port,
)

app = Celery(
    main='askcos_site',
    broker=rabbit_url,
    backend=redis_url,
)

# Using a string here means the worker don't have to serialize
# the configuration object to child processes.
app.config_from_object('askcos_site.askcos_celery.celeryconfig')


if __name__ == '__main__':
    app.start()
