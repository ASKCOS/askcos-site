from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.viewsets import ViewSet

from askcos_site.globals import db_client
from askcos_site.main.models import SavedResults

results_collection = db_client['results']['results']


class ResultsViewSet(ViewSet):
    """
    A ViewSet for accessing a user's job results.
    """

    def list(self, request):
        """Get all results belonging to the authenticated user."""
        resp = {'results': []}

        results = SavedResults.objects.filter(user=request.user)

        for result in results:
            resp['results'].append({
                'id': result.result_id,
                'state': result.result_state,
                'description': result.description,
                'created': result.created.strftime('%B %d, %Y %H:%M:%S %p'),
                'type': result.result_type,
            })

        return Response(resp)

    def retrieve(self, request, pk):
        """Get a particular result instance."""
        resp = {'id': pk, 'result': None, 'error': None}

        try:
            result = SavedResults.objects.get(user=request.user, result_id=pk)
        except SavedResults.DoesNotExist:
            resp['error'] = 'Result not found!'
            return Response(resp, status=404)
        else:
            if result.result_state == 'completed':
                result_doc = results_collection.find_one({'_id': pk})
                resp['result'] = result_doc
            else:
                resp['error'] = 'Job not yet complete.'

        return Response(resp)

    @action(detail=True, methods=['GET'])
    def check(self, request, pk):
        """Get status of a particular result instance."""
        resp = {'id': pk, 'state': None, 'error': None}

        try:
            result = SavedResults.objects.get(user=request.user, result_id=pk)
        except SavedResults.DoesNotExist:
            resp['error'] = 'Result not found!'
            return Response(resp, status=404)
        else:
            resp['state'] = result.result_state

        return Response(resp)

    def destroy(self, request, pk):
        """Delete a particular result instance."""
        resp = {'success': False, 'error': None}

        try:
            result = SavedResults.objects.get(user=request.user, result_id=pk)
        except SavedResults.DoesNotExist:
            resp['error'] = 'Result not found!'
            return Response(resp, status=404)
        else:
            try:
                result.delete()
                results_collection.delete_one({'_id': pk})
            except:
                resp['error'] = 'Could not delete result.'

        resp['success'] = True

        return Response(resp)
