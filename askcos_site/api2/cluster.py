from rest_framework import serializers
from rest_framework.generics import GenericAPIView
from rest_framework.response import Response

from makeit.utilities.cluster import group_results


class ClusterSerializer(serializers.Serializer):
    """Serializer for clustering task parameters."""
    original = serializers.CharField()
    outcomes = serializers.ListField(child=serializers.CharField())

    feature = serializers.CharField(default='original')
    fingerprint = serializers.CharField(default='morgan')
    fpradius = serializers.IntegerField(default=1)
    fpnbits = serializers.IntegerField(default=512)
    clustermethod = serializers.CharField(default='kmeans')
    scores = serializers.ListField(child=serializers.FloatField(), required=False)

    def validate_feature(self, value):
        """Check that the feature parameter is valid."""
        if value not in ['original', 'outcomes', 'all']:
            raise serializers.ValidationError("feature should be one of ['original', 'outcomes', 'all'].")
        return value

    def validate_fingerprint(self, value):
        """Check that the fingerprint parameter is valid."""
        if value not in ['morgan']:
            raise serializers.ValidationError("fingerprint should be one of ['morgan'].")
        return value

    def validate_clustermethod(self, value):
        """Check that the clustermethod parameter is valid."""
        if value not in ['hdbscan', 'kmeans']:
            raise serializers.ValidationError("clustermethod should be one of ['hdbscan', 'kmeans'].")
        return value

    def validate(self, attrs):
        """Object level validation for multiple fields."""
        if 'scores' in attrs and len(attrs['scores']) != len(attrs['outcomes']):
            raise serializers.ValidationError('Length of scores ({0}) not equal to length of outcomes ({1}).'
                                              .format(len(attrs['scores']), len(attrs['outcomes'])))
        return attrs


class ClusterAPIView(GenericAPIView):
    """
    API endpoint for clustering similar transformed outcomes

    Method: GET

    Parameters:
        original (str): smiles string
        outcomes (list): list of smiles strings of outcomes

    Optional Parameters:
        feature (str): features to use [original', 'outcomes', 'all']
        fingerprint (str): fingerprint type ['morgan']
        fpradius (int): fingerprint radius, default 1
        fpnbits (int): fingerprint bits, default 512
        clustermethod (str): cluster method ['hdbscan', 'kmeans']
        scores (list): list of scores of precursors

    Returns:
        request: dictionary of request parameters
        output: list of cluster indices for outcomes

    Test case:
        %3B=';'
        curl -k 'https://localhost/api/cluster/' -d 'original=CCOC&outcomes=CCO%3BCC'
    """

    serializer_class = ClusterSerializer

    def post(self, request, *args, **kwargs):
        """
        Handle POST requests for clustering task.
        """
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        idx = group_results(
            data['original'],
            data['outcomes'],
            feature=data['feature'],
            fp_type=data['fingerprint'],
            fp_length=data['fpnbits'],
            fp_radius=data['fpradius'],
            cluster_method=data['clustermethod'],
            scores=data.get('scores')
        )

        resp = {
            'request': data,
            'output': idx,
        }

        return Response(resp)


cluster = ClusterAPIView.as_view()
