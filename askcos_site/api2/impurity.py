from rdkit import Chem
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import serializers

from askcos_site.askcos_celery.impurity.impurity_worker import get_impurities

TIMEOUT = 30


class ImpurityPredictorSerializer(serializers.Serializer):
    """Serializer for impurity prediction task parameters."""
    reactants = serializers.CharField()
    reagents = serializers.CharField(default='')
    products = serializers.CharField(default='')
    solvent = serializers.CharField(default='')
    top_k = serializers.IntegerField(default=3)
    threshold = serializers.FloatField(default=0.75)
    predictor = serializers.CharField(default='WLN forward predictor')
    inspector = serializers.CharField(default='Reaxys predictor')
    mapper = serializers.CharField(default='WLN atom mapper')
    check_mapping = serializers.BooleanField(default=True)

    def check_smiles(self, value):
        """Verify that the requested SMILES string is valid. Returns canonicalized SMILES."""
        if value:
            mol = Chem.MolFromSmiles(value)
            if not mol:
                raise serializers.ValidationError('Cannot parse smiles with rdkit.')
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        return value

    def validate_reactants(self, value):
        """Verify that the requested reactants are valid. Returns canonicalized SMILES."""
        return self.check_smiles(value)

    def validate_reagents(self, value):
        """Verify that the requested reagents are valid. Returns canonicalized SMILES."""
        return self.check_smiles(value)

    def validate_products(self, value):
        """Verify that the requested products are valid. Returns canonicalized SMILES."""
        return self.check_smiles(value)

    def validate_solvent(self, value):
        """Verify that the requested solvent is valid. Returns canonicalized SMILES."""
        return self.check_smiles(value)


@api_view(['POST'])
def impurity_predict(request):
    """API endpoint for impurity prediction task."""
    serializer = ImpurityPredictorSerializer(data=request.data)
    serializer.is_valid(raise_exception=True)
    data = serializer.validated_data

    resp = {'request': data}

    # TODO: need work
    res = get_impurities.delay(
        data['reactants'],
        reagents=data['reagents'],
        products=data['products'],
        solvents=data['solvent'],
        top_k=data['top_k'],
        threshold=data['threshold'],
        predictor_selection=data['predictor'],
        inspector_selection=data['inspector'],
        mapper_selection=data['mapper'],
        check_mapping=data['check_mapping']
    )

    resp['task_id'] = res.task_id

    return Response(resp)
