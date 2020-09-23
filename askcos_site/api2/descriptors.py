from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.descriptors.ts_descriptors_worker import get_descriptors
from .celery import CeleryTaskAPIView

ALLOWED_ATOMS = ['C', 'H', 'O', 'N', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B']


class DescriptorsSerializer(serializers.Serializer):
    """Serializer for descriptor-predictor task parameters."""
    smiles = serializers.CharField()

    def validate_smiles(self, value):
        """Verify that the requested reactants are valid."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')

        for a in mol.GetAtoms():
            if a.GetSymbol() not in ALLOWED_ATOMS:
                raise serializers.ValidationError('invalid element {} found in smiles {}'.format(a.GetSymbol(), value))
            if a.GetFormalCharge() != 0 or a.GetNumRadicalElectrons() != 0:
                raise serializers.ValidationError('Charge or radical found in smiles {}'.format(value))

        return value


class DescriptorAPIView(CeleryTaskAPIView):
    """
    API endpoint for descriptor-predictor prediction task.

    Method: POST

    Parameters:

    - `smiles` (str): SMILES string

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = DescriptorsSerializer

    def execute(self, request, data):
        """
        Execute fast filter task and return celery result object.
        """
        result = get_descriptors.delay(data['smiles'])
        return result


descriptors = DescriptorAPIView.as_view()
