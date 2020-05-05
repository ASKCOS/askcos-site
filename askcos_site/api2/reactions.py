from rest_framework import serializers
from rest_framework.generics import GenericAPIView
from rest_framework.response import Response

from askcos_site.globals import reaction_db


class ReactionsSerializer(serializers.Serializer):
    """Serializer for reaction lookup parameters."""
    ids = serializers.ListField(child=serializers.CharField())
    template_set = serializers.CharField(required=False)


class ReactionsAPIView(GenericAPIView):
    """
    API endpoint for reaction lookup task.

    Method: POST

    Parameters:

    - `ids` (list): list of reaction ids to retrieve
    - `template_set` (str, optional): template set to search within

    Returns:

    - `reactions`: list of reactions
    """

    serializer_class = ReactionsSerializer

    def post(self, request, *args, **kwargs):
        """
        Handle POST requests for reactions lookup.
        """
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        query = {'reaction_id': {'$in': data['ids']}}
        if 'template_set' in data:
            query['template_set'] = data['template_set']

        reactions_by_ids = list(reaction_db.find(query))
        resp = {'reactions': reactions_by_ids}

        return Response(resp)


reactions = ReactionsAPIView.as_view()
