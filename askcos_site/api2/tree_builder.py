from django.utils import timezone
from rdkit import Chem
from rest_framework import serializers
from rest_framework.exceptions import NotAuthenticated

from askcos_site.main.models import BlacklistedReactions, BlacklistedChemicals, SavedResults
from askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts import get_buyable_paths as get_buyable_paths_v1
from .celery import CeleryTaskAPIView


class TreeBuilderSerializer(serializers.Serializer):
    """Serializer for tree builder task parameters."""
    smiles = serializers.CharField()
    version = serializers.IntegerField(default=1)
    max_depth = serializers.IntegerField(default=4)
    max_branching = serializers.IntegerField(default=25)
    expansion_time = serializers.IntegerField(default=60)
    template_count = serializers.IntegerField(default=100)
    max_cum_prob = serializers.FloatField(default=0.995)

    buyable_logic = serializers.CharField(default='none')
    max_ppg_logic = serializers.CharField(default='none')
    max_ppg = serializers.IntegerField(required=False)

    max_scscore_logic = serializers.CharField(default='none')
    max_scscore = serializers.FloatField(required=False)

    chemical_property_logic = serializers.CharField(default='none')
    max_chemprop_c = serializers.IntegerField(required=False)
    max_chemprop_n = serializers.IntegerField(required=False)
    max_chemprop_o = serializers.IntegerField(required=False)
    max_chemprop_h = serializers.IntegerField(required=False)

    chemical_popularity_logic = serializers.CharField(default='none')
    min_chempop_reactants = serializers.IntegerField(required=False)
    min_chempop_products = serializers.IntegerField(required=False)

    filter_threshold = serializers.FloatField(default=0.75)
    template_set = serializers.CharField(default='reaxys')
    template_prioritizer_version = serializers.IntegerField(default=0)
    buyables_source = serializers.ListField(child=serializers.CharField(allow_blank=True), required=False, allow_empty=True)
    return_first = serializers.BooleanField(default=True)
    max_trees = serializers.IntegerField(default=500)

    store_results = serializers.BooleanField(default=False)
    description = serializers.CharField(default='')

    banned_reactions = serializers.ListField(child=serializers.CharField(), required=False)
    banned_chemicals = serializers.ListField(child=serializers.CharField(), required=False)

    priority = serializers.IntegerField(default=1)

    def validate_smiles(self, value):
        """Verify that the requested smiles is valid. Returns canonicalized SMILES."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse smiles with rdkit.')
        return Chem.MolToSmiles(mol)

    def validate_buyable_logic(self, value):
        """Verify that the the specified buyable_logic is valid."""
        if value not in ['none', 'and', 'or']:
            raise serializers.ValidationError("Logic should be one of ['none', 'and', 'or'].")
        return value

    def validate_max_ppg_logic(self, value):
        """Verify that the the specified max_ppg_logic is valid."""
        if value not in ['none', 'and', 'or']:
            raise serializers.ValidationError("Logic should be one of ['none', 'and', 'or'].")
        return value

    def validate_max_scscore_logic(self, value):
        """Verify that the the specified max_scscore_logic is valid."""
        if value not in ['none', 'and', 'or']:
            raise serializers.ValidationError("Logic should be one of ['none', 'and', 'or'].")
        return value

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

    def validate_banned_chemicals(self, value):
        """
        Verify that the provided SMILES is valid. Returns canonicalized SMILES.
        """
        new_value = []
        for v in value:
            mol = Chem.MolFromSmiles(v)
            if not mol:
                raise serializers.ValidationError('Cannot parse smiles with rdkit.')
            new_value.append(Chem.MolToSmiles(mol, isomericSmiles=True))
        return new_value

    def validate_banned_reactions(self, value):
        """
        Verify that the provided SMILES is valid. Returns canonicalized SMILES.
        """
        new_value = []
        for v in value:
            try:
                reactants, agents, products = v.split('>')
            except ValueError:
                raise serializers.ValidationError('Cannot parse reaction smiles.')
            try:
                reactants = standardize(reactants, isomericSmiles=True)
            except ValueError:
                raise serializers.ValidationError('Cannot parse reaction reactants.')
            try:
                agents = standardize(agents, isomericSmiles=True)
            except ValueError:
                raise serializers.ValidationError('Cannot parse reaction agents.')
            try:
                products = standardize(products, isomericSmiles=True)
            except ValueError:
                raise serializers.ValidationError('Cannot parse reaction products.')
            new_value.append(reactants + '>' + agents + '>' + products)
        return new_value


def standardize(smiles, isomericSmiles=True):
    """
    Split input SMILES into individual molecules, canonicalizes each, then
    sorts and re-combines canonicalized SMILES into single SMILES.
    """
    parts = smiles.split('.')
    canonicalized_parts = []
    for part in parts:
        mol = Chem.MolFromSmiles(part)
        if not mol:
            raise ValueError()
        canonicalized_parts.append(Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles))
    canonicalized_parts.sort()
    return '.'.join(canonicalized_parts)


class TreeBuilderAPIView(CeleryTaskAPIView):
    """
    API endpoint for tree builder prediction task.

    Method: POST

    Parameters:

    - `smiles` (str): SMILES string of target
    - `version` (int): tree builder version to use
    - `max_depth` (int, optional): maximum depth of returned pathways
    - `max_branching` (int, optional): maximum branching during pathway exploration
    - `expansion_time` (int, optional): time limit for tree expansion
    - `template_count` (int, optional): number of templates to consider
    - `max_cum_prob` (float, optional): maximum cumulative probability of templates
    - `buyable_logic` (str, optional): logic type for buyable termination (none/and/or)
    - `max_ppg_logic` (str, optional): logic type for price based termination (none/and/or)
    - `max_ppg` (int, optional): maximum price for price based termination
    - `max_scscore_logic` (str, optional): logic type for synthetic complexity termination (none/and/or)
    - `max_scscore` (int, optional): maximum scscore for synthetic complexity termination
    - `chemical_property_logic` (str, optional): logic type for chemical property termination (none/and/or)
    - `max_chemprop_c` (int, optional): maximum carbon count for termination
    - `max_chemprop_n` (int, optional): maximum nitrogen count for termination
    - `max_chemprop_o` (int, optional): maximum oxygen count for termination
    - `max_chemprop_h` (int, optional): maximum hydrogen count for termination
    - `chemical_popularity_logic` (str, optional): logic type for chemical popularity termination (none/and/or)
    - `min_chempop_reactants` (int, optional): minimum reactant precedents for termination
    - `min_chempop_products` (int, optional): minimum product precedents for termination
    - `filter_threshold` (float, optional): fast filter threshold
    - `template_set` (str, optional): template set to use
    - `template_prioritizer_version` (int, optional): version number of template relevance model to use
    - `buyables_source` (str, optional): source(s) to consider when looking up buyables (accepts comma delimited list)
    - `return_first` (bool, optional): whether to return upon finding the first pathway
    - `max_trees` (int, optional): maximum number of pathways to return
    - `store_results` (bool, optional): whether to permanently save this result
    - `description` (str, optional): description to associate with stored result
    - `banned_reactions` (list, optional): list of reactions to not consider
    - `banned_chemicals` (list, optional): list of molecules to not consider
    - `priority` (int, optional): set priority for celery task (0 = low, 1 = normal (default), 2 = high)

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = TreeBuilderSerializer

    def execute(self, request, data):
        """
        Execute tree builder task and return celery result object.
        """
        if data['store_results'] and not request.user.is_authenticated:
            raise NotAuthenticated('You must be authenticated to store tree builder results.')

        termination_logic = {'and': [], 'or': []}

        buyable_logic = data['buyable_logic']
        if buyable_logic != 'none':
            termination_logic[buyable_logic].append('buyable')

        max_ppg_logic = data['max_ppg_logic']
        if max_ppg_logic != 'none':
            max_ppg = data.get('max_ppg')
            termination_logic[max_ppg_logic].append('max_ppg')
        else:
            max_ppg = None

        max_scscore_logic = data['max_scscore_logic']
        if max_scscore_logic != 'none':
            max_scscore = data.get('max_scscore')
            termination_logic[max_scscore_logic].append('max_scscore')
        else:
            max_scscore = None

        chemical_property_logic = data['chemical_property_logic']
        if chemical_property_logic != 'none':
            param_dict = {
                'C': 'max_chemprop_c',
                'N': 'max_chemprop_n',
                'O': 'max_chemprop_o',
                'H': 'max_chemprop_h',
            }
            max_elements = {k: data[v] for k, v in param_dict.items() if v in data}
            termination_logic[chemical_property_logic].append('max_elements')
        else:
            max_elements = None

        chemical_popularity_logic = data['chemical_popularity_logic']
        if chemical_popularity_logic != 'none':
            min_history = {
                'as_reactant': data.get('min_chempop_reactants', 5),
                'as_product': data.get('min_chempop_products', 5),
            }
            termination_logic[chemical_popularity_logic].append('min_history')
        else:
            min_history = None

        # Clean buyables source
        buyables_source = data.get('buyables_source')
        if buyables_source is not None and 'none' in buyables_source:
            # Include both null and empty string source in query
            buyables_source.remove('none')
            buyables_source.extend([None, ''])

        # Retrieve user specific banlists
        banned_reactions = data.get('banned_reactions', [])
        banned_chemicals = data.get('banned_chemicals', [])
        if request.user.is_authenticated:
            banned_reactions += list(set(
                [x.smiles for x in BlacklistedReactions.objects.filter(user=request.user, active=True)]))
            banned_chemicals += list(set(
                [x.smiles for x in BlacklistedChemicals.objects.filter(user=request.user, active=True)]))

        args = (data['smiles'],)
        kwargs = {
            'max_depth': data['max_depth'],
            'max_branching': data['max_branching'],
            'expansion_time': data['expansion_time'],
            'template_count': data['template_count'],
            'max_cum_template_prob': data['max_cum_prob'],
            'max_ppg': max_ppg,
            'max_scscore': max_scscore,
            'max_elements': max_elements,
            'min_history': min_history,
            'termination_logic': termination_logic,
            'filter_threshold': data['filter_threshold'],
            'template_set': data['template_set'],
            'template_prioritizer_version': data['template_prioritizer_version'],
            'buyables_source': buyables_source,
            'known_bad_reactions': banned_reactions,
            'forbidden_molecules': banned_chemicals,
            'return_first': data['return_first'],
            'max_trees': data['max_trees'],
            'run_async': data['store_results'],
            'paths_only': True,
        }

        if data['version'] == 1:
            result = get_buyable_paths_v1.apply_async(args, kwargs, priority=data['priority'])

        if data['store_results']:
            now = timezone.now()
            saved_result = SavedResults.objects.create(
                user=request.user,
                created=now,
                dt=now.strftime('%B %d, %Y %H:%M:%S %p'),
                result_id=result.id,
                result_state='pending',
                result_type='tree_builder',
                description=data['description']
            )

        return result


tree_builder = TreeBuilderAPIView.as_view()
