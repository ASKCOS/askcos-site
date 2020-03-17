from rdkit import Chem
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import serializers

from askcos_site.globals import scscorer


class SCScorerSerializer(serializers.Serializer):
    """Serializer for SCScorer parameters."""
    smiles = serializers.CharField()

    def validate_smiles(self, value):
        """Verify that the requested smiles is valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse smiles with rdkit.')
        return value


@api_view(['POST'])
def scscore(request):
    """API endpoint for scscore prediction task."""
    serializer = SCScorerSerializer(data=request.data)
    serializer.is_valid(raise_exception=True)
    data = serializer.validated_data

    resp = {
        'request': data,
        'score': scscorer.get_score_from_smiles(data['smiles'], noprice=True),
    }

    return Response(resp)
