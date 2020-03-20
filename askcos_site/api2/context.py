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
    async = serializers.BooleanField(default=True)

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
    - `async` (bool, optional): whether to directly return celery task id instead of waiting for result

    Returns:

    - `output`: list of reaction conditions
    """

    serializer_class = ContextRecommenderSerializer

    def execute(self, data):
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
        )

        return result

    def process(self, data, output):
        """
        Post-process output from context recommendation task.
        """
        if data['return_scores']:
            contexts, scores = output
        else:
            contexts = output

        json_contexts = []
        for context in contexts:
            json_contexts.append({
                'temperature': context[0],
                'solvent': context[1],
                'reagent': context[2],
                'catalyst': context[3]
            })

        if data['return_scores']:
            for c, s in zip(json_contexts, scores):
                c['score'] = s

        return json_contexts


neural_network = ContextRecommenderAPIView.as_view()
