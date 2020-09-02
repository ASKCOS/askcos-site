# Default serializers
task_serializer = 'json'
result_serializer = 'json'

# Allowed content types - other message types are discarded
accept_content = ['json']

# Timezone for message dates and times (set to match django settings)
timezone = 'UTC'

# Convert dates and times to UTC
enable_utc = True

# Time that results remain in queue
result_expires = 1800  # 30 min
result_persistent = False

# Interval between sending worker heartbeats
broker_heartbeat = 0

# Maximum allowed priority (larger values are capped)
# Having more priority levels will consume more CPU resources
# Max priority of 2 enables 3 levels - 0, 1, 2
task_queue_max_priority = 2

# Default priority for tasks
task_default_priority = 1

# Number of messages to prefetch from queue
# Prefetching fewer tasks gives queue more opportunity to reorder by priority
worker_prefetch_multiplier = 1

# Modules to be imported by each worker
imports = [
    'askcos_site.askcos_celery.treebuilder.tb_c_worker',
    'askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts',
    'askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts_v2',
    'askcos_site.askcos_celery.treebuilder.path_ranking_worker',
    'askcos_site.askcos_celery.contextrecommender.cr_network_worker',
    'askcos_site.askcos_celery.treeevaluator.template_free_forward_predictor_worker',
    'askcos_site.askcos_celery.siteselectivity.sites_worker',
    'askcos_site.askcos_celery.impurity.impurity_worker',
    'askcos_site.askcos_celery.impurity.impurity_predictor_worker',
    'askcos_site.askcos_celery.atom_mapper.atom_mapping_worker',
    'askcos_site.askcos_celery.generalselectivity.selec_worker',
]

# Task routes (to make sure workers are task-specific)
# Key is pattern matched against task name to determine queue
task_routes = {
    'askcos_site.askcos_celery.treebuilder.tb_c_worker.*': {'queue': 'tb_c_worker'},
    'askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts.*': {'queue': 'tb_coordinator_mcts'},
    'askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts_v2.*': {'queue': 'tb_coordinator_mcts_v2'},
    'askcos_site.askcos_celery.treebuilder.path_ranking_worker.*': {'queue': 'path_ranking_worker'},
    'askcos_site.askcos_celery.contextrecommender.cr_network_worker.*': {'queue': 'cr_network_worker'},
    'askcos_site.askcos_celery.treeevaluator.template_free_forward_predictor_worker.*': {'queue': 'tffp_worker'},
    'askcos_site.askcos_celery.siteselectivity.sites_worker.*': {'queue': 'sites_worker'},
    'askcos_site.askcos_celery.impurity.impurity_worker.*': {'queue': 'impurity_worker'},
    'askcos_site.askcos_celery.impurity.impurity_predictor_worker.*': {'queue': 'atom_mapping_worker'},
    'askcos_site.askcos_celery.atom_mapper.atom_mapping_worker.*': {'queue': 'atom_mapping_worker'},
    'askcos_site.askcos_celery.generalselectivity.selec_worker.*': {'queue': 'selec_worker'},
}
