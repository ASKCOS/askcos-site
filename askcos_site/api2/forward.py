from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.treeevaluator.template_free_forward_predictor_worker import get_outcomes
from .celery import CeleryTaskAPIView


class ForwardPredictorSerializer(serializers.Serializer):
    """Serializer for fast-filter task parameters."""
    reactants = serializers.CharField()
    reagents = serializers.CharField(default='')
    solvent = serializers.CharField(default='')
    num_results = serializers.IntegerField(default=100)
    atommap = serializers.BooleanField(default=False)

    def validate_reactants(self, value):
        """Verify that the requested reactants are valid. Returns canonicalized SMILES."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
        return Chem.MolToSmiles(mol, isomericSmiles=True)

    def validate_reagents(self, value):
        """Verify that the requested reagents are valid. Returns canonicalized SMILES."""
        if value:
            mol = Chem.MolFromSmiles(value)
            if not mol:
                raise serializers.ValidationError('Cannot parse reagents smiles with rdkit.')
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        return value

    def validate_solvent(self, value):
        """Verify that the requested solvent are valid. Returns canonicalized SMILES."""
        if value:
            mol = Chem.MolFromSmiles(value)
            if not mol:
                raise serializers.ValidationError('Cannot parse solvent smiles with rdkit.')
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        return value


class ForwardPredictorAPIView(CeleryTaskAPIView):
    """
    API endpoint for template-free forward prediction task.

    Method: POST

    Parameters:

    - `reactants` (str): SMILES string of reactants
    - `reagents` (str, optional): SMILES string of reagents (default='')
    - `solvent` (str, optional): SMILES string of solvent (default='')
    - `num_results` (int, optional): max number of results to return (default=100)
    - `atommap` (bool, optional): Flag to keep atom mapping from the prediction (default=False)

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = ForwardPredictorSerializer

    def execute(self, request, data):
        """
        Execute forward prediction task and return celery result object.
        """
        reactants = data['reactants']
        reagents = data['reagents']
        solvent = data['solvent']
        num_results = data['num_results']
        atommap = data['atommap']

        combined_smiles = reactants
        if reagents:
            combined_smiles += '.{}'.format(reagents)
        if solvent:
            combined_smiles += '.{}'.format(solvent)

        result = get_outcomes.delay(combined_smiles, top_n=num_results, atommap=atommap)

        return result


template_free = ForwardPredictorAPIView.as_view()
