from celery.exceptions import TimeoutError
from rdkit import Chem
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import serializers

from askcos_site.askcos_celery.siteselectivity.sites_worker import get_sites


class SelectivitySerializer(serializers.Serializer):
    """Serializer for site selectivity task parameters."""
    smiles = serializers.CharField()

    def validate_smiles(self, value):
        """Verify that the requested smiles is valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse smiles with rdkit.')
        return value


@api_view(['POST'])
def selectivity(request):
    """API endpoint for site selectivity prediction task."""
    serializer = SelectivitySerializer(data=request.data)
    serializer.is_valid(raise_exception=True)
    data = serializer.validated_data

    resp = {'request': data}

    res = get_sites.delay(data['smiles'])

    try:
        result = res.get(30)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit 30s)'
        res.revoke()
        return Response(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return Response(resp, status=400)

    resp['result'] = result

    return Response(resp)
