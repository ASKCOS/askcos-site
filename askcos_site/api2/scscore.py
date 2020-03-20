from rdkit import Chem
from rest_framework import serializers
from rest_framework.generics import GenericAPIView
from rest_framework.response import Response

from askcos_site.globals import scscorer


class SCScorerSerializer(serializers.Serializer):
    """Serializer for SCScorer parameters."""
    smiles = serializers.CharField()

    def validate_smiles(self, value):
        """Verify that the requested smiles is valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse smiles with rdkit.')
        return value


class SCScorerAPIView(GenericAPIView):
    """
    API endpoint for scscore prediction task.
    """

    serializer_class = SCScorerSerializer

    def post(self, request, *args, **kwargs):
        """
        Handle POST requests for scscore prediction.
        """
        serializer = SCScorerSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        resp = {
            'request': data,
            'score': scscorer.get_score_from_smiles(data['smiles'], noprice=True),
        }

        return Response(resp)


scscore = SCScorerAPIView.as_view()
