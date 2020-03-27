from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.atom_mapper.atom_mapping_worker import get_atom_mapping
from .celery import CeleryTaskAPIView


class AtomMapperSerializer(serializers.Serializer):
    """Serializer for atom mapping task parameters."""
    rxnsmiles = serializers.CharField()
    mapper = serializers.CharField(default='WLN atom mapper')

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


class AtomMapperAPIView(CeleryTaskAPIView):
    """
    API endpoint for atom mapping task.

    Method: POST

    Parameters:

    - `rxnsmiles` (str): reaction SMILES string
    - `mapper` (str, optional): atom mapping backend to use (currently only 'WLN atom mapper')

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = AtomMapperSerializer

    def execute(self, data):
        """
        Execute fast filter task and return celery result object.
        """
        result = get_atom_mapping.delay(data['rxnsmiles'], mapper=data['mapper'])
        return result


atom_mapper = AtomMapperAPIView.as_view()
