import requests
from rdkit import Chem
from rest_framework import serializers
from rest_framework.generics import GenericAPIView
from rest_framework.response import Response

from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors
from .celery import CeleryTaskAPIView


class RetroSerializer(serializers.Serializer):
    """Serializer for retrosynthesis task parameters."""
    target = serializers.CharField()
    num_templates = serializers.IntegerField(default=100)
    max_cum_prob = serializers.FloatField(min_value=0.0, max_value=1.0, default=0.995)
    filter_threshold = serializers.FloatField(default=0.75)
    template_set = serializers.CharField(default='reaxys')
    template_prioritizer_version = serializers.IntegerField(default=0)

    cluster = serializers.BooleanField(default=True)
    cluster_method = serializers.CharField(default='kmeans')
    cluster_feature = serializers.CharField(default='original')
    cluster_fp_type = serializers.CharField(default='morgan')
    cluster_fp_length = serializers.IntegerField(default=512)
    cluster_fp_radius = serializers.IntegerField(default=1)

    selec_check = serializers.BooleanField(default=True)

    def validate_target(self, value):
        """Verify that the requested target is valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse target smiles with rdkit.')
        return value


class RetroModelsSerializer(serializers.Serializer):
    """Serializer for drawing parameters."""
    template_set = serializers.CharField()


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
    - `template_prioritizer` (str, optional): template prioritization model to use
    - `cluster` (bool, optional): whether or not to cluster results
    - `cluster_method` (str, optional): method for clustering results
    - `cluster_feature` (str, optional): which feature to use for clustering
    - `cluster_fp_type` (str, optional): fingerprint type for clustering
    - `cluster_fp_length` (int, optional): fingerprint length for clustering
    - `cluster_fp_radius` (int, optional): fingerprint radius for clustering
    - `selec_check` (bool, optional): whether or not to check for potential selectivity issues

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = RetroSerializer

    def execute(self, request, data):
        """
        Execute single step retro task and return celery result object.
        """
        target = data['target']
        max_num_templates = data['num_templates']
        max_cum_prob = data['max_cum_prob']
        fast_filter_threshold = data['filter_threshold']
        template_set = data['template_set']
        template_prioritizer_version = data['template_prioritizer_version']

        cluster = data['cluster']
        cluster_method = data['cluster_method']
        cluster_feature = data['cluster_feature']
        cluster_fp_type = data['cluster_fp_type']
        cluster_fp_length = data['cluster_fp_length']
        cluster_fp_radius = data['cluster_fp_radius']

        selec_check = data['selec_check']

        result = get_top_precursors.delay(
            target,
            template_set=template_set,
            template_prioritizer_version=template_prioritizer_version,
            fast_filter_threshold=fast_filter_threshold,
            max_cum_prob=max_cum_prob,
            max_num_templates=max_num_templates,
            cluster=cluster,
            cluster_method=cluster_method,
            cluster_feature=cluster_feature,
            cluster_fp_type=cluster_fp_type,
            cluster_fp_length=cluster_fp_length,
            cluster_fp_radius=cluster_fp_radius,
            selec_check=selec_check,
            postprocess=True,
        )

        return result


class RetroModels(GenericAPIView):
    """
    API endpoint for querying available retrosynthetic models for a given template set.

    Method: GET

    Parameters:

    - `template_set` (str): template set name

    Returns:

    - `versions`: List of version numbers that are available
    """

    serializer_class = RetroModelsSerializer

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

        versions = [
            model.get('version')
            for model in model_version_status
            if model.get('state') == 'AVAILABLE'
        ]

        resp = {
            'request': data,
            'versions': versions
        }

        return Response(resp)


models = RetroModels.as_view()
singlestep = RetroAPIView.as_view()
