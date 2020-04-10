from rest_framework.generics import GenericAPIView
from rest_framework.response import Response
from rest_framework.reverse import reverse


class RootAPIView(GenericAPIView):
    """
    API endpoint for listing all API endpoints.

    Method: GET

    Returns: name and url of all API endpoints
    """
    name = 'API Root'

    def get(self, request, *args, **kwargs):
        """
        Handle GET requests for API endpoint list.
        """
        resp = {
            '/celery': reverse('v2:celery_api', request=request),
            '/celery/task': reverse('v2:celery_task_api-list', request=request),
            '/atom-mapper': reverse('v2:atom_mapper_api', request=request),
            '/context': reverse('v2:context_api', request=request),
            '/fast-filter': reverse('v2:fast_filter_api', request=request),
            '/forward': reverse('v2:forward_api', request=request),
            '/impurity': reverse('v2:impurity_api', request=request),
            '/retro': reverse('v2:retro_api', request=request),
            '/selectivity': reverse('v2:selectivity_api', request=request),
            '/tree-builder': reverse('v2:tree_builder_api', request=request),
            '/cluster': reverse('v2:cluster_api', request=request),
            '/reactions': reverse('v2:reactions_api', request=request),
            '/scscore': reverse('v2:scscore_api', request=request),
            '/rdkit/smiles/canonicalize': reverse('v2:smiles_api-canonicalize', request=request),
            '/rdkit/smiles/validate': reverse('v2:smiles_api-validate', request=request),
            '/rdkit/smiles/from-molfile': reverse('v2:smiles_api-from-molfile', request=request),
            '/rdkit/smiles/to-molfile': reverse('v2:smiles_api-to-molfile', request=request),
            '/buyables': reverse('v2:buyables_api-list', request=request),
            '/buyables/upload': reverse('v2:buyables_api-upload', request=request),
            '/template': reverse('v2:template_api-list', request=request),
            '/results': reverse('v2:results_api-list', request=request),
            '/blacklist/chemicals': reverse('v2:blacklist_chemicals_api-list', request=request),
            '/blacklist/reactions': reverse('v2:blacklist_reactions_api-list', request=request),
            '/token-auth': reverse('v2:token_auth_api', request=request),
            '/token-refresh': reverse('v2:token_refresh_api', request=request),
        }

        return Response(resp)


apiroot = RootAPIView.as_view()
