from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import serializers

from askcos_site.globals import reaction_db


class ReactionsSerializer(serializers.Serializer):
    """Serializer for reaction lookup parameters."""
    ids = serializers.ListField(child=serializers.CharField())
    template_set = serializers.CharField(required=False)


@api_view(['POST'])
def reactions(request):
    """API endpoint for reaction lookup task."""
    serializer = ReactionsSerializer(data=request.data)
    serializer.is_valid(raise_exception=True)
    data = serializer.validated_data

    query = {'reaction_id': {'$in': data['ids']}}
    if 'template_set' in data:
        query['template_set'] = data['template_set']

    reactions_by_ids = list(reaction_db.find(query))
    resp = {'reactions': reactions_by_ids}

    return Response(resp)
