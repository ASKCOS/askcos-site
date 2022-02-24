from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.atom_mapper.atom_mapping_worker import get_atom_mapping
from .celery import CeleryTaskAPIView


class AtomMapperSerializer(serializers.Serializer):
    """Serializer for atom mapping task parameters."""
    rxnsmiles = serializers.CharField()
    mapper = serializers.CharField(default='WLN atom mapper')
    priority = serializers.IntegerField(default=1)

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
    - `priority` (int, optional): set priority for celery task (0 = low, 1 = normal (default), 2 = high)

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = AtomMapperSerializer

    def execute(self, request, data):
        """
        Execute fast filter task and return celery result object.
        """
        args = (data['rxnsmiles'],)
        kwargs = {
            'mapper': data['mapper'],
        }
        result = get_atom_mapping.apply_async(args, kwargs, priority=data['priority'])
        return result


atom_mapper = AtomMapperAPIView.as_view()
