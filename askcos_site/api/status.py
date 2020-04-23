from django.http import JsonResponse
from askcos_site.celery import app, READABLE_NAMES
from celery.result import AsyncResult
from celery.exceptions import TimeoutError


def celery_status(request):
    resp = {}
    status = {}
    stats = app.control.inspect().stats()
    active = app.control.inspect().active()
    if not stats or not active:
        return JsonResponse(resp)
    worker_names = stats.keys()
    for worker in worker_names:
        name, server = worker.split('@')
        if not status.get(name):
            status[name] = {'available': 0, 'busy': 0}
        status[name]['busy'] += len(active[worker])
        status[name]['available'] += stats[worker]['pool']['max-concurrency'] - status[name]['busy']
    status_list = []
    for key in status:
        status_list.append({
            'name': READABLE_NAMES.get(key, key),
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
    return JsonResponse(resp)


def task_status(request):
    resp = {}
    task_id = request.GET.get('task_id')
    if not task_id:
        resp['error'] = 'Must provide "task_id" to API call'
        return JsonResponse(resp, status=400)
    
    result = AsyncResult(task_id)
    if not result:
        resp['error'] = 'Cannot find task with task_id: {}'.format(task_id)
        return JsonResponse(resp, status=400)

    state = result.state
    
    try:
        info = result.info
        resp['percent'] = info.get('percent')
        resp['message'] = info.get('message')
    except AttributeError:
        # info is weird, unsure how to handle it
        pass


    if state == 'running' or state=='PENDING':
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
            outcomes = result.get(10) # should not take very long to get results of a completed task
        except Exception as e:
            resp['error'] = str(e)
            result.revoke()
            return JsonResponse(resp, status=400)

        resp['results'] = outcomes
    
    return JsonResponse(resp)