from rdkit import Chem
from rest_framework import serializers

from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors as get_top_precursors_c
from askcos_site.askcos_celery.treebuilder.tb_c_worker_preload import get_top_precursors as get_top_precursors_p
from .celery import CeleryTaskAPIView


class RetroSerializer(serializers.Serializer):
    """Serializer for retrosynthesis task parameters."""
    target = serializers.CharField()
    num_templates = serializers.IntegerField(default=100)
    max_cum_prob = serializers.FloatField(min_value=0.0, max_value=1.0, default=0.995)
    filter_threshold = serializers.FloatField(default=0.75)
    template_set = serializers.CharField(default='reaxys')
    template_prioritizer = serializers.CharField(default='reaxys')
    async = serializers.BooleanField(default=True)

    cluster = serializers.BooleanField(default=True)
    cluster_method = serializers.CharField(default='kmeans')
    cluster_feature = serializers.CharField(default='original')
    cluster_fp_type = serializers.CharField(default='morgan')
    cluster_fp_length = serializers.IntegerField(default=512)
    cluster_fp_radius = serializers.IntegerField(default=1)

    def validate_template_set(self, value):
        """Verify that the requested template set is valid."""
        if value not in ['reaxys', 'uspto_50k']:
            raise serializers.ValidationError('Template set {} not available.'.format(value))
        return value

    def validate_template_prioritizer(self, value):
        """Verify that the requested template prioritizer is valid."""
        if value not in ['reaxys', 'uspto_50k']:
            raise serializers.ValidationError('Template prioritizer {} not available.'.format(value))
        return value

    def validate_target(self, value):
        """Verify that the requested target is valid."""
        if not Chem.MolFromSmiles(value):
            raise serializers.ValidationError('Cannot parse target smiles with rdkit.')
        return value


class RetroAPIView(CeleryTaskAPIView):
    """API endpoint for single-step retrosynthesis task."""

    serializer_class = RetroSerializer
    TIMEOUT = 120

    def execute(self, data):
        """
        Execute single step retro task and return celery result object.
        """
        target = data['target']
        max_num_templates = data['num_templates']
        max_cum_prob = data['max_cum_prob']
        fast_filter_threshold = data['filter_threshold']
        template_set = data['template_set']
        template_prioritizer = data['template_prioritizer']

        cluster = data['cluster']
        cluster_method = data['cluster_method']
        cluster_feature = data['cluster_feature']
        cluster_fp_type = data['cluster_fp_type']
        cluster_fp_length = data['cluster_fp_length']
        cluster_fp_radius = data['cluster_fp_radius']

        if max_cum_prob > 0.999 and max_num_templates > 1000:
            result = get_top_precursors_p.delay(
                target,
                template_set=template_set,
                template_prioritizer=template_prioritizer,
                fast_filter_threshold=fast_filter_threshold,
                max_cum_prob=max_cum_prob,
                max_num_templates=max_num_templates,
                cluster=cluster,
                cluster_method=cluster_method,
                cluster_feature=cluster_feature,
                cluster_fp_type=cluster_fp_type,
                cluster_fp_length=cluster_fp_length,
                cluster_fp_radius=cluster_fp_radius
            )
        else:
            result = get_top_precursors_c.delay(
                target,
                template_set=template_set,
                template_prioritizer=template_prioritizer,
                fast_filter_threshold=fast_filter_threshold,
                max_cum_prob=max_cum_prob,
                max_num_templates=max_num_templates,
                cluster=cluster,
                cluster_method=cluster_method,
                cluster_feature=cluster_feature,
                cluster_fp_type=cluster_fp_type,
                cluster_fp_length=cluster_fp_length,
                cluster_fp_radius=cluster_fp_radius
            )

        return result

    def process(self, data, output):
        """
        Post-process output from single step retrosynthesis task.
        """
        smiles, precursors = output

        for precursor in precursors:
            precursor['templates'] = precursor.pop('tforms')

        return precursors


singlestep = RetroAPIView.as_view()
