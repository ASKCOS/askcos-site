from celery.exceptions import TimeoutError
from rdkit import Chem
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import serializers

from askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts import get_buyable_paths as get_buyable_paths_mcts


class TreeBuilderSerializer(serializers.Serializer):
    """Serializer for tree builder task parameters."""
    smiles = serializers.CharField()
    async = serializers.BooleanField(default=True)
    max_depth = serializers.IntegerField(default=4)
    max_branching = serializers.IntegerField(default=25)
    expansion_time = serializers.IntegerField(default=60)
    max_ppg = serializers.IntegerField(default=10)
    template_count = serializers.IntegerField(default=100)
    max_cum_prob = serializers.FloatField(default=0.995)

    chemical_property_logic = serializers.CharField(default='none')
    max_chemprop_c = serializers.IntegerField(required=False)
    max_chemprop_n = serializers.IntegerField(required=False)
    max_chemprop_o = serializers.IntegerField(required=False)
    max_chemprop_h = serializers.IntegerField(required=False)

    chemical_popularity_logic = serializers.CharField(default='none')
    min_chempop_reactants = serializers.IntegerField(required=False)
    min_chempop_products = serializers.IntegerField(required=False)

    filter_threshold = serializers.FloatField(default=0.75)
    template_prioritizer = serializers.CharField(default='reaxys')
    template_set = serializers.CharField(default='reaxys')
    hashed_historian = serializers.BooleanField(required=False)
    return_first = serializers.BooleanField(default=True)

    blacklisted_reactions = serializers.ListField(child=serializers.CharField(), required=False)
    forbidden_molecules = serializers.ListField(child=serializers.CharField(), required=False)

    def validate_smiles(self, value):
        """Verify that the requested smiles is valid. Returns canonicalized SMILES."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse smiles with rdkit.')
        return Chem.MolToSmiles(mol)

    def validate_chemical_property_logic(self, value):
        """Verify that the the specified chemical_property_logic is valid."""
        if value not in ['none', 'and', 'or']:
            raise serializers.ValidationError("Logic should be one of ['none', 'and', 'or'].")
        return value

    def validate_chemical_popularity_logic(self, value):
        """Verify that the the specified chemical_popularity_logic is valid."""
        if value not in ['none', 'and', 'or']:
            raise serializers.ValidationError("Logic should be one of ['none', 'and', 'or'].")
        return value


@api_view(['POST'])
def tree_builder(request):
    """API endpoint for tree builder prediction task."""
    serializer = TreeBuilderSerializer(data=request.data)
    serializer.is_valid(raise_exception=True)
    data = serializer.validated_data

    chemical_property_logic = data['chemical_property_logic']
    if chemical_property_logic != 'none':
        param_dict = {
            'C': 'max_chemprop_c',
            'N': 'max_chemprop_n',
            'O': 'max_chemprop_o',
            'H': 'max_chemprop_h',
        }
        max_natom_dict = {k: data[v] for k, v in param_dict if v in data}
        max_natom_dict['logic'] = chemical_property_logic
    else:
        max_natom_dict = None

    chemical_popularity_logic = data['chemical_popularity_logic']
    if chemical_popularity_logic != 'none':
        min_chemical_history_dict = {
            'logic': chemical_popularity_logic,
            'as_reactant': data.get('min_chempop_reactants', 5),
            'as_product': data.get('min_chempop_products', 5),
        }
    else:
        min_chemical_history_dict = None

    res = get_buyable_paths_mcts.delay(
        data['smiles'],
        max_depth=data['max_depth'],
        max_branching=data['max_branching'],
        expansion_time=data['expansion_time'],
        max_trees=500,
        max_ppg=data['max_ppg'],
        known_bad_reactions=data.get('blacklisted_reactions'),
        forbidden_molecules=data.get('forbidden_molecules'),
        template_count=data['template_count'],
        max_cum_template_prob=data['max_cum_prob'],
        max_natom_dict=max_natom_dict,
        min_chemical_history_dict=min_chemical_history_dict,
        apply_fast_filter=data['filter_threshold'] > 0,
        filter_threshold=data['filter_threshold'],
        template_prioritizer=data['template_prioritizer'],
        template_set=data['template_set'],
        hashed=data.get('hashed_historian', data['template_set'] == 'reaxys'),
        return_first=data['return_first'],
    )

    resp = {'request': data}

    if data['async']:
        resp['id'] = res.id
        resp['state'] = res.state
        return Response(resp)

    try:
        (tree_status, trees) = res.get(data['expansion_time'] * 3)
    except TimeoutError:
        resp['error'] = 'API request timed out (after {})'.format(data['expansion_time'] * 3)
        res.revoke()
        return Response(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return Response(resp, status=400)

    resp['trees'] = trees

    return Response(resp)
