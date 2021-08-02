from django.utils import timezone
from rdkit import Chem
from rest_framework import serializers
from rest_framework.exceptions import NotAuthenticated

from askcos_site.main.models import BlacklistedReactions, BlacklistedChemicals, SavedResults
from askcos_site.askcos_celery.treebuilder.tb_coordinator_mcts import get_buyable_paths as get_buyable_paths_mcts
from askcos_site.main.utils import is_banned
from .celery import CeleryTaskAPIView


class TreeBuilderSerializer(serializers.Serializer):
    """Serializer for tree builder task parameters."""
    smiles = serializers.CharField()
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
    template_set = serializers.CharField(default='reaxys')
    template_prioritizer_version = serializers.IntegerField(default=0)
    return_first = serializers.BooleanField(default=True)

    store_results = serializers.BooleanField(default=False)
    description = serializers.CharField(default='')

    banned_reactions = serializers.ListField(child=serializers.CharField(), required=False)
    banned_chemicals = serializers.ListField(child=serializers.CharField(), required=False)

    def validate_smiles(self, value):
        """Verify that the requested smiles is valid. Returns canonicalized SMILES."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse smiles with rdkit.')
        if is_banned(self.context['request'], value):
            raise serializers.ValidationError('ASKCOS does not provide results for compounds on restricted lists such as the CWC and DEA schedules.')
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
    - `max_depth` (int, optional): maximum depth of returned trees
    - `max_branching` (int, optional): maximum branching in returned trees
    - `expansion_time` (int, optional): time limit for tree expansion
    - `max_ppg` (int, optional): maximum price for buyable termination
    - `template_count` (int, optional): number of templates to consider
    - `max_cum_prob` (float, optional): maximum cumulative probability of templates
    - `chemical_property_logic` (str, optional): logic type for chemical property termination
    - `max_chemprop_c` (int, optional): maximum carbon count for termination
    - `max_chemprop_n` (int, optional): maximum nitrogen count for termination
    - `max_chemprop_o` (int, optional): maximum oxygen count for termination
    - `max_chemprop_h` (int, optional): maximum hydrogen count for termination
    - `chemical_popularity_logic` (str, optional): logic type for chemical popularity termination
    - `min_chempop_reactants` (int, optional): minimum reactant precedents for termination
    - `min_chempop_products` (int, optional): minimum product precedents for termination
    - `filter_threshold` (float, optional): fast filter threshold
    - `template_set` (str, optional): template set to use
    - `template_prioritizer_version` (int, optional): version number of template relevance model to use
    - `return_first` (bool, optional): whether to return upon finding the first pathway
    - `store_results` (bool, optional): whether to permanently save this result
    - `description` (str, optional): description to associate with stored result
    - `banned_reactions` (list, optional): list of reactions to not consider
    - `banned_chemicals` (list, optional): list of molecules to not consider

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

        chemical_property_logic = data['chemical_property_logic']
        if chemical_property_logic != 'none':
            param_dict = {
                'C': 'max_chemprop_c',
                'N': 'max_chemprop_n',
                'O': 'max_chemprop_o',
                'H': 'max_chemprop_h',
            }
            max_natom_dict = {k: data[v] for k, v in param_dict.items() if v in data}
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

        # Retrieve user specific banlists
        banned_reactions = data.get('banned_reactions', [])
        banned_chemicals = data.get('banned_chemicals', [])
        if request.user.is_authenticated:
            banned_reactions += list(set(
                [x.smiles for x in BlacklistedReactions.objects.filter(user=request.user, active=True)]))
            banned_chemicals += list(set(
                [x.smiles for x in BlacklistedChemicals.objects.filter(user=request.user, active=True)]))

        result = get_buyable_paths_mcts.delay(
            data['smiles'],
            max_depth=data['max_depth'],
            max_branching=data['max_branching'],
            expansion_time=data['expansion_time'],
            max_trees=500,
            max_ppg=data['max_ppg'],
            known_bad_reactions=banned_reactions,
            forbidden_molecules=banned_chemicals,
            template_count=data['template_count'],
            max_cum_template_prob=data['max_cum_prob'],
            max_natom_dict=max_natom_dict,
            min_chemical_history_dict=min_chemical_history_dict,
            apply_fast_filter=data['filter_threshold'] > 0,
            filter_threshold=data['filter_threshold'],
            template_prioritizer_version=data['template_prioritizer_version'],
            template_set=data['template_set'],
            return_first=data['return_first'],
            paths_only=True,
            run_async=data['store_results']
        )

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
