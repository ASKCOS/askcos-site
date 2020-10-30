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
    namespace = ''
    urlpatterns = []

    def get(self, request, *args, **kwargs):
        """
        Handle GET requests for API endpoint list.
        """
        endpoints = {}

        for url in self.urlpatterns:
            name = self.namespace + ':' + url.name
            to_skip = ['-detail', '-export', '-check', '-activate', '-deactivate', '-ipp', '-tree']
            if any(name.endswith(suffix) for suffix in to_skip):
                # Don't include any detail views in the endpoint list
                continue

            pattern = str(url.pattern).strip('^').strip('$').strip('/')
            if not pattern:
                # Don't include root endpoint in the list
                continue

            endpoints[pattern] = reverse(name, request=request)

        return Response(endpoints)
