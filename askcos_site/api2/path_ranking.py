import json

from rest_framework import serializers

from askcos_site.askcos_celery.treebuilder.path_ranking_worker import get_scores
from .celery import CeleryTaskAPIView


class PathRankingSerializer(serializers.Serializer):
    """Serializer for path ranking task parameters."""
    trees = serializers.CharField()
    cluster = serializers.BooleanField(default=True)
    cluster_method = serializers.CharField(default='hdbscan')
    min_samples = serializers.IntegerField(default=5)
    min_cluster_size = serializers.IntegerField(default=5)


class PathRankingAPIView(CeleryTaskAPIView):
    """
    API endpoint for fast-filter prediction task.

    Method: POST

    Parameters:

    - `trees` (str): list of trees to rank as a json string
    - `cluster` (bool, optional): whether or not to cluster pathways
    - `cluster_method` (str, optional): hdbscan or kmeans
    - `min_samples` (int, optional): min samples for hdbscan
    - `min_cluster_size` (int, optional): min cluster size for hdbscan

    Returns:

    - `task_id`: celery task ID
    """

    serializer_class = PathRankingSerializer

    def execute(self, request, data):
        """
        Execute fast filter task and return celery result object.
        """
        trees = json.loads(data.pop('trees'))
        args = (trees,)
        kwargs = {
            'clustering': data['cluster'],
            'cluster_method': data['cluster_method'],
            'min_samples': data['min_samples'],
            'min_cluster_size': data['min_cluster_size'],
        }
        result = get_scores.apply_async(args, kwargs)
        return result


path_ranker = PathRankingAPIView.as_view()
