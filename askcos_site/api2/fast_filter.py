from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.treebuilder.tb_c_worker import fast_filter_check
from .celery import CeleryTaskAPIView


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


class FastFilterAPIView(CeleryTaskAPIView):
    """
    API endpoint for fast-filter prediction task.

    Method: POST

    Parameters:

    - `reactants` (str): SMILES string of reactants
    - `products` (str): SMILES string of products

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = FastFilterSerializer

    def execute(self, request, data):
        """
        Execute fast filter task and return celery result object.
        """
        result = fast_filter_check.delay(data['reactants'], data['products'])
        return result


fast_filter = FastFilterAPIView.as_view()
