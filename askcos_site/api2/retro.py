import requests
from rdkit import Chem
from rest_framework import serializers
from rest_framework.generics import GenericAPIView
from rest_framework.response import Response

from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors, template_relevance, apply_one_template_by_idx
from .celery import CeleryTaskAPIView


class AttributeFilterSerializer(serializers.Serializer):
    """Serializer for individual attribute filter object"""
    name = serializers.CharField()
    logic = serializers.CharField()
    value = serializers.FloatField()

    def validate_logic(self, value):
        if value not in ['>', '>=', '<', '<=', '==']:
            raise serializers.ValidationError('Attribute filter logic "{}" not supported.'.format(value))
        return value


class RetroSerializer(serializers.Serializer):
    """Serializer for retrosynthesis task parameters."""
    target = serializers.CharField()
    num_templates = serializers.IntegerField(default=100)
    max_cum_prob = serializers.FloatField(min_value=0.0, max_value=1.0, default=0.995)
    filter_threshold = serializers.FloatField(default=0.75)
    template_set = serializers.CharField(default='reaxys')
    template_prioritizer_version = serializers.IntegerField(default=0)
    precursor_prioritizer = serializers.CharField(default='RelevanceHeuristic')

    cluster = serializers.BooleanField(default=True)
    cluster_method = serializers.CharField(default='kmeans')
    cluster_feature = serializers.CharField(default='original')
    cluster_fp_type = serializers.CharField(default='morgan')
    cluster_fp_length = serializers.IntegerField(default=512)
    cluster_fp_radius = serializers.IntegerField(default=1)

    selec_check = serializers.BooleanField(default=True)

    attribute_filter = AttributeFilterSerializer(default=[], many=True)

    priority = serializers.IntegerField(default=1)

    def validate_target(self, value):
        """Verify that the requested target is valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse target smiles with rdkit.')
        return value


class TFXRetroModelsSerializer(serializers.Serializer):
    """Serializer for available retro models parameters."""
    template_set = serializers.CharField()


class TemplateRelevanceSerializer(serializers.Serializer):
    """Serializer for template relevance prediction parameters"""
    smiles = serializers.CharField()
    num_templates = serializers.IntegerField(default=100)
    max_cum_prob = serializers.FloatField(min_value=0.0, max_value=1.0, default=0.995)
    template_set = serializers.CharField(default='reaxys')
    template_prioritizer_version = serializers.IntegerField(default=0)

    def validate_smiles(self, value):
        """Verify that the requested target is valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse target smiles with rdkit.')
        return value


class ApplyOneTemplateByIdxSerializer(serializers.Serializer):
    """Serializer for applying one template by index"""
    smiles = serializers.CharField()
    template_idx = serializers.IntegerField()
    template_set = serializers.CharField(default='reaxys')

    def validate_smiles(self, value):
        """Verify that the requested target is valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse target smiles with rdkit.')
        return value


class RetroAPIView(CeleryTaskAPIView):
    """
    API endpoint for single-step retrosynthesis task.

    Method: POST

    Parameters:

    - `target` (str): SMILES string of target
    - `num_templates` (int, optional): number of templates to consider
    - `max_cum_prob` (float, optional): maximum cumulative probability of templates
    - `filter_threshold` (float, optional): fast filter threshold
    - `template_set` (str, optional): reaction template set to use
    - `template_prioritizer_version` (int, optional): version number of template relevance model to use
    - `precursor_prioritizer` (str, optional): name of precursor prioritizer to use (Relevanceheuristic or SCScore)
    - `cluster` (bool, optional): whether or not to cluster results
    - `cluster_method` (str, optional): method for clustering results
    - `cluster_feature` (str, optional): which feature to use for clustering
    - `cluster_fp_type` (str, optional): fingerprint type for clustering
    - `cluster_fp_length` (int, optional): fingerprint length for clustering
    - `cluster_fp_radius` (int, optional): fingerprint radius for clustering
    - `selec_check` (bool, optional): whether or not to check for potential selectivity issues
    - `attribute_filter` (list[dict], optional): template attribute filter to apply before template application
    - `priority` (int, optional): set priority for celery task (0 = low, 1 = normal (default), 2 = high)

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = RetroSerializer

    def execute(self, request, data):
        """
        Execute single step retro task and return celery result object.
        """
        args = (data['target'],)
        kwargs = {
            'max_num_templates': data['num_templates'],
            'max_cum_prob': data['max_cum_prob'],
            'fast_filter_threshold': data['filter_threshold'],
            'template_set': data['template_set'],
            'template_prioritizer_version': data['template_prioritizer_version'],
            'precursor_prioritizer': data['precursor_prioritizer'],
            'cluster': data['cluster'],
            'cluster_method': data['cluster_method'],
            'cluster_feature': data['cluster_feature'],
            'cluster_fp_type': data['cluster_fp_type'],
            'cluster_fp_length': data['cluster_fp_length'],
            'cluster_fp_radius': data['cluster_fp_radius'],
            'selec_check': data['selec_check'],
            'attribute_filter': data['attribute_filter'],
            'postprocess': True,
        }

        result = get_top_precursors.apply_async(args, kwargs, priority=data['priority'])

        return result


class TFXRetroModels(GenericAPIView):
    """
    API endpoint for querying available retrosynthetic models for a given template set.

    Method: GET

    Parameters:

    - `template_set` (str): template set name

    Returns:

    - `versions`: List of version numbers that are available
    """

    serializer_class = TFXRetroModelsSerializer

    def get(self, request, *args, **kwargs):
        """
        Handle GET requests for retro models endpoint.
        """
        serializer = self.get_serializer(data=request.query_params)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        url = 'http://template-relevance-{}:8501/v1/models/template_relevance'.format(data['template_set'])
        try:
            api_resp = requests.get(url)
        except requests.exceptions.ConnectionError:
            resp = {'request': data, 'error': 'tensorflow serving model(s) not available for {}'.format(data['template_set'])}
            return Response(resp)

        model_version_status = api_resp.json().get('model_version_status')
        if not model_version_status:
            resp = {'request': data, 'error': 'tensorflow serving model(s) not available for {}'.format(data['template_set'])}
            return Response(resp)

        versions = sorted([
            model.get('version')
            for model in model_version_status
            if model.get('state') == 'AVAILABLE'
        ])

        resp = {
            'request': data,
            'versions': versions
        }

        return Response(resp)


class TemplateRelevanceAPIView(CeleryTaskAPIView):
    """
    API endpoint for a template relevance prediction.

    Method: POST

    Parameters:

    - `smiles` (str): target smiles
    - `num_templates` (int): number of templates for which to return predicted scores
    - `max_cum_prob` (float): maximum cumulative probability of templates to return
    - `template_set` (str): template set name
    - `template_prioritizer_version` (int): template relevance model version number

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = TemplateRelevanceSerializer

    def execute(self, request, data):
        """
        Execute template relevance prediction via celery.
        """
        args = (data['smiles'], data['num_templates'], data['max_cum_prob'])
        kwargs = {
            'template_set': data['template_set'],
            'template_prioritizer_version': data['template_prioritizer_version'],
        }

        result = template_relevance.apply_async(args, kwargs)

        return result


class ApplyOneTemplateByIdxAPIView(CeleryTaskAPIView):
    """
    API endpoint for applying one template by index (and template set).

    Method: POST

    Parameters:

    - `smiles` (str): target smiles
    - `template_idx` (int): template index
    - `template_set` (str): name of template set

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = ApplyOneTemplateByIdxSerializer

    def execute(self, request, data):
        """
        Execute template relevance prediction via celery.
        """
        # First arg is tree builder path ID which is not needed
        args = (0, data['smiles'], data['template_idx'])
        kwargs = {
            'calculate_next_probs': False,
            'template_set': data['template_set'],
        }

        result = apply_one_template_by_idx.apply_async(args, kwargs)

        # Unpack result into dict
        # First item is the tree builder path ID which is not needed
        result = {
            'smiles': result[1],
            'template_idx': result[2],
            'precursors': result[3],
            'ffscore': result[4],
        }

        return result


models = TFXRetroModels.as_view()
singlestep = RetroAPIView.as_view()
temprel = TemplateRelevanceAPIView.as_view()
apply_one_template = ApplyOneTemplateByIdxAPIView.as_view()
