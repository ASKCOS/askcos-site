from celery.exceptions import TimeoutError
from celery.result import AsyncResult
from rest_framework.generics import GenericAPIView
from rest_framework.response import Response
from rest_framework.viewsets import ViewSet

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


class CeleryTaskAPIView(GenericAPIView):
    """
    Base API view for a celery task.

    Method: POST

    Returns:

    - `task_id`: ID of celery task if async request
    - `output`: output of celery task if not async request
    """

    TIMEOUT = 30

    def post(self, request, *args, **kwargs):
        """
        Handle POST requests for a generic celery task endpoint.
        """
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        result = self.execute(data)

        resp = {'request': data}

        if data['async']:
            resp['task_id'] = result.id
            return Response(resp)

        try:
            output = result.get(self.TIMEOUT)
        except TimeoutError:
            resp['error'] = 'API request timed out (limit {}s).'.format(self.TIMEOUT)
            result.revoke()
            return Response(resp, status=408)
        except Exception as e:
            resp['error'] = str(e)
            result.revoke()
            return Response(resp, status=400)

        resp['output'] = self.process(data, output)

        return Response(resp)

    def execute(self, data):
        """
        Execute the celery task and return a celery result object.
        """
        raise NotImplementedError('Should be implemented by child class.')

    def process(self, data, output):
        """
        Post-process output from a celery task.
        """
        return output


class CeleryStatusAPIView(GenericAPIView):
    """
    API endpoint for retrieving celery worker status.

    Method: GET

    Returns:

    - `queues`: list of worker information for each celery queue
    """

    def get(self, request, *args, **kwargs):
        """
        Handle GET requests for scscore prediction.
        """
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


class CeleryTaskViewSet(ViewSet):
    """
    API endpoint for retrieving status and result of a celery task.

    For a particular task, specified as URI parameter (/api/v2/celery/task/<task id>/):

    Method: GET

    Returns:

    - `complete`: boolean indicating whether job is complete
    - `failed`: boolean indicating if job failed
    - `percent`: completion percent of job
    - `message`: message regarding job status
    - `state`: state of the job
    - `error`: any error message if encountered
    - `results`: result of celery task if complete
    """

    def retrieve(self, request, pk):
        """Get the status and result of a single celery task by task_id."""
        resp = {}

        result = AsyncResult(pk)
        if not result:
            resp['error'] = 'Cannot find task with task_id: {}'.format(pk)
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


celery_status = CeleryStatusAPIView.as_view()
