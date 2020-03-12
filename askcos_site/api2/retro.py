from rdkit import Chem
from celery.exceptions import TimeoutError
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import serializers

from askcos_site.askcos_celery.treebuilder.tb_c_worker import get_top_precursors as get_top_precursors_c
from askcos_site.askcos_celery.treebuilder.tb_c_worker_preload import get_top_precursors as get_top_precursors_p

TIMEOUT = 120


class RetroSerializer(serializers.Serializer):
    """Serializer for retrosynthesis task parameters."""
    async = serializers.BooleanField(default=False)
    target = serializers.CharField()
    num_templates = serializers.IntegerField(default=100)
    max_cum_prob = serializers.FloatField(min_value=0.0, max_value=1.0, default=0.995)
    filter_threshold = serializers.FloatField(default=0.75)
    template_set = serializers.CharField(default='reaxys')
    template_prioritizer = serializers.CharField(default='reaxys')

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


@api_view(['POST'])
def singlestep(request):
    """API endpoint for single-step retrosynthesis task."""
    serializer = RetroSerializer(data=request.data)
    serializer.is_valid(raise_exception=True)
    data = serializer.validated_data

    run_async = data['async']
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

    resp = {'request': data}

    if max_cum_prob > 0.999 and max_num_templates > 1000:
        res = get_top_precursors_p.delay(
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
        res = get_top_precursors_c.delay(
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

    if run_async:
        resp['id'] = res.id
        resp['state'] = res.state
        return Response(resp)

    try:
        (smiles, precursors) = res.get(TIMEOUT)
    except TimeoutError:
        resp['error'] = 'API request timed out (limit {}s)'.format(TIMEOUT)
        res.revoke()
        return Response(resp, status=408)
    except Exception as e:
        resp['error'] = str(e)
        res.revoke()
        return Response(resp, status=400)

    resp['precursors'] = precursors
    for precursor in precursors:
        precursor['templates'] = precursor.pop('tforms')

    return Response(resp)
