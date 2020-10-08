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
    mapped = serializers.BooleanField(default=False)
    all_outcomes = serializers.BooleanField(default=False)
    verbose = serializers.BooleanField(default=True)
    mapper = serializers.CharField(default='Transformer')
    no_map_reagents = serializers.BooleanField(default=True)

    def validate_reactants(self, value):
        """Verify that the requested reactants smiles is valid."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
        return Chem.MolToSmiles(mol, isomericSmiles=True)

    def validate_reagents(self, value):
        """Verify that the requested reagents smiles is valid."""
        if value:
            mol = Chem.MolFromSmiles(value)
            if not mol:
                raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        return value

    def validate_solvent(self, value):
        """Verify that the requested solvent smiles is valid."""
        if value:
            mol = Chem.MolFromSmiles(value)
            if not mol:
                raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        return value

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
        mapped = data['mapped']
        all_outcomes = data['all_outcomes']
        verbose = data['verbose']
        mapper = data['mapper']
        no_map_reagents = data['no_map_reagents']

        combined_smiles = reactants + '>'
        if reagents and solvent:
            combined_smiles += '{}.{}'.format(reagents, solvent)
        elif reagents:
            combined_smiles += '{}'.format(reagents)
        elif solvent:
            combined_smiles += '{}'.format(solvent)
        combined_smiles += '>{}'.format(product)

        result = get_selec.delay(combined_smiles, mapped=mapped, all_outcomes=all_outcomes, verbose=verbose, mapper=mapper, no_map_reagents=no_map_reagents)
        return result


selectivity = SelectivityAPIView.as_view()
