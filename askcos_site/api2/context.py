from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.contextrecommender.cr_network_worker import get_n_conditions as network_get_n_conditions
from .celery import CeleryTaskAPIView


class ContextRecommenderSerializer(serializers.Serializer):
    """Serializer for context recommendation task parameters."""
    reactants = serializers.CharField()
    products = serializers.CharField()
    with_smiles = serializers.BooleanField(default=True)
    single_solvent = serializers.BooleanField(default=True)
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


class ContextRecommenderAPIView(CeleryTaskAPIView):
    """
    API endpoint for context recommendation prediction using neural network model.

    Method: POST

    Parameters:

    - `reactants` (str): SMILES string of reactants
    - `products` (str): SMILES string of products
    - `with_smiles` (bool, optional): whether to use SMILES for prediction
    - `single_solvent` (bool, optional): whether to use single solvent for prediction
    - `return_scores` (bool, optional): whether to also return scores
    - `num_results` (int, optional): max number of results to return

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = ContextRecommenderSerializer

    def execute(self, request, data):
        """
        Execute context recommendation task.
        """
        rxn = data['reactants'] + '>>' + data['products']

        result = network_get_n_conditions.delay(
            rxn,
            n=data['num_results'],
            singleSlvt=data['single_solvent'],
            with_smiles=data['with_smiles'],
            return_scores=data['return_scores'],
            postprocess=True,
        )

        return result


neural_network = ContextRecommenderAPIView.as_view()
