from celery.result import AsyncResult
from rest_framework import serializers
from rest_framework.decorators import api_view
from rest_framework.response import Response

from askcos_site.celery import app

READABLE_NAMES = {
    'cr_network_worker': 'Context Recommender Worker',
    'tb_c_worker': 'One-Step/Tree Builder Retrosynthesis Worker',
    'tb_c_worker_preload': 'One-Step/Tree Builder Retrosynthesis Worker (Pre-loaded)',
    'tb_coordinator_mcts': 'Tree Builder Coordinator',
    'sites_worker': 'Site Selectivity Worker',
    'impurity_worker': 'Impurity worker',
    'atom_mapping_worker': 'Atom mapping worker',
    'tffp_worker': 'Template-free Forward Predictor'
}


class TaskSerializer(serializers.Serializer):
    """Serializer for celery task lookup."""
    task_id = serializers.CharField()


@api_view(['GET'])
def celery_status(request):
    """API endpoint for retrieving celery worker status."""

    resp = {}
    status = {}
    stats = app.control.inspect().stats()
    active = app.control.inspect().active()

    if not stats or not active:
        return Response(resp)

    for worker in stats:
        name, server = worker.split('@')
        if not status.get(name):
            status[name] = {'available': 0, 'busy': 0}
        status[name]['busy'] += len(active[worker])
        status[name]['available'] += stats[worker]['pool']['max-concurrency'] - status[name]['busy']

    status_list = []
    for key in status:
        status_list.append({
            'name': READABLE_NAMES.get(key),
            'queue': key,
            'busy': status[key]['busy'],
            'available': status[key]['available']
        })

    for key, val in READABLE_NAMES.items():
        if key not in status:
            status_list.append({
                'name': READABLE_NAMES.get(key),
                'queue': key,
                'busy': 0,
                'available': 0
            })

    resp['queues'] = sorted(status_list, key=lambda x: x['name'])

    return Response(resp)


@api_view(['GET'])
def task_status(request):
    """API endpoint for retrieving celery task status."""
    serializer = TaskSerializer(data=request.query_params)
    serializer.is_valid(raise_exception=True)
    task_id = serializer.validated_data['task_id']

    resp = {}

    result = AsyncResult(task_id)
    if not result:
        resp['error'] = 'Cannot find task with task_id: {}'.format(task_id)
        return Response(resp, status=400)

    state = result.state

    try:
        info = result.info
        resp['percent'] = info.get('percent')
        resp['message'] = info.get('message')
    except AttributeError:
        # info is weird, unsure how to handle it
        pass

    if state == 'running' or state == 'PENDING':
        resp['complete'] = False
    elif state == 'failed' or state == 'FAILURE':
        resp['complete'] = False
        resp['failed'] = True
    else:
        resp['state'] = state
        resp['complete'] = True
        resp['percent'] = 1
        resp['message'] = 'Task complete!'
        try:
            outcomes = result.get(10)  # should not take very long to get results of a completed task
        except Exception as e:
            resp['error'] = str(e)
            result.revoke()
            return Response(resp, status=400)

        resp['results'] = outcomes

    return Response(resp)
