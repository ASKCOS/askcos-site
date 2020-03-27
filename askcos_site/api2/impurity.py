from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.impurity.impurity_worker import get_impurities
from .celery import CeleryTaskAPIView


class ImpurityPredictorSerializer(serializers.Serializer):
    """Serializer for impurity prediction task parameters."""
    reactants = serializers.CharField()
    reagents = serializers.CharField(default='')
    products = serializers.CharField(default='')
    solvent = serializers.CharField(default='')
    top_k = serializers.IntegerField(default=3)
    threshold = serializers.FloatField(default=0.75)
    predictor = serializers.CharField(default='WLN forward predictor')
    inspector = serializers.CharField(default='Reaxys inspector')
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


class ImpurityAPIView(CeleryTaskAPIView):
    """
    API endpoint for impurity prediction task.

    Method: POST

    Parameters:

    - `reactants` (str): SMILES string of reactants
    - `reagents` (str, optional): SMILES string of reagents
    - `products` (str, optional): SMILES string of products
    - `solvent` (str, optional): SMILES string of solvent
    - `top_k` (int, optional): max number of results to return
    - `threshold` (float, optional): probability threshold
    - `predictor` (str, optional): forward predictor to use
    - `inspector` (str, optional): reaction scorer to use
    - `mapper` (str, optional): reaction atom mapper to use
    - `check_mapping` (bool, optional): whether to check atom mapping

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = ImpurityPredictorSerializer

    def execute(self, data):
        """
        Execute an impurity task and return the celery result object.
        """
        result = get_impurities.delay(
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

        return result


impurity_predict = ImpurityAPIView.as_view()
