from rdkit import Chem
from rest_framework import serializers

import askcos.global_config as gc
from askcos_site.askcos_celery.contextrecommender.cr_network_v2_worker import get_n_conditions
from .celery import CeleryTaskAPIView


class StringListField(serializers.ListField):
    child = serializers.CharField()
    allow_empty = True


def validate_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        return False
    if mol is None:
        return False
    return True


def validate_smiles_list(smiles_list):
    for smiles in smiles_list:
        if not validate_smiles(smiles):
            return False
    return True


def has_atom_map(smiles):
    '''Only perform basic check, i.e. mapno exists
    '''
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        return False
    if mol is None:
        return False
    
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() != 0:
            return True
    return False


class ContextRecommenderSerializer(serializers.Serializer):
    """Serializer for context recommendation task parameters."""
    reactants = serializers.CharField()
    products = serializers.CharField()
    reagents = StringListField(default=None)
    num_results = serializers.IntegerField(default=10)
    model = serializers.CharField(default='graph')

    def validate(self, data):
        """Verify that the requested reaction is atom-mapped for graph models."""
        if data['model'].startswith('graph'):
            if not has_atom_map(data['reactants']):
                raise serializers.ValidationError('Reactants smiles does not contains atom mapping.')
            if not has_atom_map(data['products']):
                raise serializers.ValidationError('Products smiles does not contains atom mapping.')
        return data

    def validate_reactants(self, value):
        """Verify that the requested reactants are valid."""
        if not validate_smiles(value):
            raise serializers.ValidationError('Cannot parse reactants smiles with rdkit.')
        return value

    def validate_products(self, value):
        """Verify that the requested products are valid."""
        if not validate_smiles(value):
            raise serializers.ValidationError('Cannot parse products smiles with rdkit.')
        return value

    def validate_reagents(self, value):
        """Verify that the input reagents are valid."""
        if value is None:
            return value
        
        if not validate_smiles_list(value):
            raise serializers.ValidationError('Cannot parse reagents smiles with rdkit.')
        # TODO: validate reagents are in the reagents list
        return value
    
    def validate_model(self, value):
        valid_models = list(gc.CONTEXT_V2['default-models'].keys()) + list(gc.CONTEXT_V2['models'].keys())
        if value not in valid_models:
            raise serializers.ValidationError('Invalid model type. Expects fp or graph.')
        if value in gc.CONTEXT_V2['default-models']:
            value = gc.CONTEXT_V2['default-models'][value]
        return value


class ContextRecommenderAPIView(CeleryTaskAPIView):
    """
    API endpoint for context recommendation prediction using neural network model.

    Method: POST

    Parameters:

    - `reactants` (str): SMILES string of reactants
    - `products` (str): SMILES string of products
    - `reagents` (list of str, optional): predefined reagents
    - `num_results` (int, optional): max number of results to return. It is also the beam width in beam search.
    - `model` (str, optional): 'fp' or 'graph'. Default is 'graph'.

    Returns:

    - `task_id`: celery task ID

    Test:
    - curl -k https://localhost/api/v2/context-v2/ -X POST -d 'reactants=[N:1]#[C:2][CH:3]([C:4](=O)[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1)[c:11]1[cH:12][cH:13][cH:14][cH:15][cH:16]1.[NH2:17][NH2:18]&products=[NH2:1][c:2]1[nH:18][n:17][c:4]([c:3]1-[c:11]1[cH:16][cH:15][cH:14][cH:13][cH:12]1)-[c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1'

    - curl -k https://localhost/api/v2/celery/task/{task_id}/
    """

    serializer_class = ContextRecommenderSerializer

    def execute(self, request, data):
        """
        Execute context recommendation task.
        """
        reaction = data['reactants'] + '>>' + data['products']
        reagents = data['reagents']
        model_name = data['model']

        result = get_n_conditions.delay(
            smiles=reaction,
            reagents=reagents,
            model_name=model_name,
            beam_size=data['num_results'],
        )

        return result


context = ContextRecommenderAPIView.as_view()
