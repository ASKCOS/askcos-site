from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.generalselectivity.selec_worker import get_selec
from .celery import CeleryTaskAPIView


class GeneralSelectivitySerializer(serializers.Serializer):
    """Serializer for selectivity task parameters."""
    reactants = serializers.CharField()
    reagents = serializers.CharField(default='')
    solvent = serializers.CharField(default='')
    product = serializers.CharField()

    def validate_reactants(self, value):
        """Verify that the requested reactants smiles is valid."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
        return Chem.MolToSmiles(mol, isomericSmiles=True)

    def validate_reagents(self, value):
        """Verify that the requested reagents smiles is valid."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
        return Chem.MolToSmiles(mol, isomericSmiles=True)

    def validate_solvent(self, value):
        """Verify that the requested solvent smiles is valid."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
        return Chem.MolToSmiles(mol, isomericSmiles=True)

    def validate_product(self, value):
        """Verify that the requested product smiles is valid."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
        return Chem.MolToSmiles(mol, isomericSmiles=True)


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
        reactants = data['reactants']
        reagents = data['reagents']
        solvent = data['solvent']
        product = data['product']

        combined_smiles = reactants
        if reagents and solvent:
            combined_smiles += '>{}.{}'.format(reagents, solvent)
        elif reagents:
            combined_smiles += '>{}'.format(reagents)
        elif solvent:
            combined_smiles += '>{}'.format(solvent)
        combined_smiles += '>{}'.format(product)

        result = get_selec.delay(combined_smiles)
        return result


selectivity = SelectivityAPIView.as_view()
