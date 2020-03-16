from celery.exceptions import TimeoutError
from rdkit import Chem
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import serializers

from askcos_site.askcos_celery.contextrecommender.cr_network_worker import get_n_conditions as network_get_n_conditions

TIMEOUT = 30


class ContextRecommenderSerializer(serializers.Serializer):
    """Serializer for context recommendation task parameters."""
    reactants = serializers.CharField()
    products = serializers.CharField()
    with_smiles = serializers.BooleanField(default=True)
    singleSlvt = serializers.BooleanField(default=True)
    return_scores = serializers.BooleanField(default=False)
    num_results = serializers.IntegerField(default=10)

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
def neural_network(request):
    """API endpoint for context recommendation prediction using neural network model."""
    serializer = ContextRecommenderSerializer(data=request.data)
    serializer.is_valid(raise_exception=True)
    data = serializer.validated_data

    reactants = data['reactants']
    products = data['products']
    with_smiles = data['with_smiles']
    singleSlvt = data['singleSlvt']
    return_scores = data['return_scores']
    n = data['num_results']
    rxn = reactants + '>>' + products

    resp = {'request': data}

    res = network_get_n_conditions.delay(rxn, n, singleSlvt, with_smiles, return_scores)

    try:
        if return_scores:
            contexts, scores = res.get(TIMEOUT)
        else:
            contexts = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(TIMEOUT)
        res.revoke()
        return Response(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return Response(resp, status=400)

    json_contexts = []
    for context in contexts:
        json_contexts.append({
            'temperature': context[0],
            'solvent': context[1],
            'reagent': context[2],
            'catalyst': context[3]
        })

    if return_scores:
        for c, s in zip(json_contexts, scores):
            c['score'] = s

    resp['contexts'] = json_contexts

    return Response(resp)
