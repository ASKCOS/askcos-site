from celery.exceptions import TimeoutError
from rdkit import Chem
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import serializers

from askcos_site.askcos_celery.treeevaluator.template_free_forward_predictor_worker import get_outcomes

TIMEOUT = 30


class ForwardPredictorSerializer(serializers.Serializer):
    """Serializer for fast-filter task parameters."""
    reactants = serializers.CharField()
    reagents = serializers.CharField(default='')
    solvent = serializers.CharField(default='')
    num_results = serializers.IntegerField(default=100)

    def validate_reactants(self, value):
        """Verify that the requested reactants are valid. Returns canonicalized SMILES."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
        return Chem.MolToSmiles(mol, isomericSmiles=True)

    def validate_reagents(self, value):
        """Verify that the requested reagents are valid. Returns canonicalized SMILES."""
        if value:
            mol = Chem.MolFromSmiles(value)
            if not mol:
                raise serializers.ValidationError('Cannot parse reagents smiles with rdkit.')
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        return value

    def validate_solvent(self, value):
        """Verify that the requested solvent are valid. Returns canonicalized SMILES."""
        if value:
            mol = Chem.MolFromSmiles(value)
            if not mol:
                raise serializers.ValidationError('Cannot parse solvent smiles with rdkit.')
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        return value


@api_view(['POST'])
def template_free(request):
    """API endpoint for template-free forward prediction task."""
    serializer = ForwardPredictorSerializer(data=request.data)
    serializer.is_valid(raise_exception=True)
    data = serializer.validated_data

    reactants = data['reactants']
    reagents = data['reagents']
    solvent = data['solvent']
    num_results = data['num_results']

    resp = {'request': data, 'error': None}

    combined_smiles = reactants
    if reagents:
        combined_smiles += '.{}'.format(reagents)
    if solvent:
        combined_smiles += '.{}'.format(solvent)

    res = get_outcomes.delay(combined_smiles, top_n=num_results)

    try:
        outcomes = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(TIMEOUT)
        res.revoke()
        return Response(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return Response(resp, status=400)

    resp['outcomes'] = outcomes

    return Response(resp)
