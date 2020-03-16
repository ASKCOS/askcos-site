from rdkit import Chem
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import serializers

from askcos_site.askcos_celery.treebuilder.tb_c_worker import fast_filter_check

TIMEOUT = 30


class FastFilterSerializer(serializers.Serializer):
    """Serializer for fast-filter task parameters."""
    reactants = serializers.CharField()
    products = serializers.CharField()

    def validate_reactants(self, value):
        """Verify that the requested reactants are valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
        return value

    def validate_products(self, value):
        """Verify that the requested products are valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse products smiles with rdkit.')
        return value


@api_view(['POST'])
def fast_filter(request):
    """API endpoint for fast-filter prediction task."""
    serializer = FastFilterSerializer(data=request.data)
    serializer.is_valid(raise_exception=True)
    data = serializer.validated_data

    reactants = data['reactants']
    products = data['products']

    resp = {'request': data, 'error': None}

    res = fast_filter_check.delay(reactants, products)

    try:
        outcome = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(TIMEOUT)
        res.revoke()
        return Response(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return Response(resp, status=400)
    
    resp['score'] = outcome

    return Response(resp)
