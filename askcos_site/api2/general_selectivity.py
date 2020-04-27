from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.generalselectivity.selec_worker import get_selec
from .celery import CeleryTaskAPIView


class GeneralSelectivitySerializer(serializers.Serializer):
    """Serializer for selectivity task parameters."""
    rxnsmiles = serializers.CharField()

    def validate_rxnsmiles(self, value):
        """Verify that the requested reaction smiles is valid."""
        try:
            reactants, agents, products = value.split('>')
        except ValueError:
            raise serializers.ValidationError('Cannot parse reaction smiles.')
        if not Chem.MolFromSmiles(reactants):
            raise serializers.ValidationError('Cannot parse reactants using rdkit.')
        if not Chem.MolFromSmiles(products):
            raise serializers.ValidationError('Cannot parse products using rdkit.')
        return value


class SelectivityAPIView(CeleryTaskAPIView):
    """
    API endpoint for general selectivity prediction task.

    Method: POST

    Parameters:

    - `reaction_smiles` (str): reaction smiles with map atom number

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = GeneralSelectivitySerializer

    def execute(self, request, data):
        """
        Execute site selectivity task and return celery result object.
        """
        result = get_selec.delay(data['rxnsmiles'])
        return result


selectivity = SelectivityAPIView.as_view()
