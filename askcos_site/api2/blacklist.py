from django.utils import timezone
from rest_framework import serializers
from rest_framework.decorators import action
from rest_framework.permissions import IsAuthenticated
from rest_framework.response import Response
from rest_framework.viewsets import ModelViewSet

from askcos_site.main.models import BlacklistedReactions, BlacklistedChemicals


class BlacklistSerializer(serializers.Serializer):
    """
    Base serializer class for blacklisted items.
    """
    smiles = serializers.CharField()
    description = serializers.CharField(required=False)
    datetime = serializers.CharField(required=False)
    active = serializers.BooleanField(default=True)

    def create(self, validated_data):
        """
        Create a new blacklisted instance.
        """
        now = timezone.now()
        return self.Meta.model.objects.create(
            user=self.context['request'].user,
            description=validated_data.get('description', 'no description'),
            created=now,
            dt=validated_data.get('datetime', now.strftime('%B %d, %Y %H:%M:%S %p UTC')),
            smiles=validated_data['smiles'],
            active=validated_data['active'],
        )

    def to_representation(self, instance):
        """
        Create serialized representation of a blacklisted instance.
        """
        ret = {
            'id': instance.id,
            'smiles': instance.smiles,
            'description': instance.description,
            'created': serializers.DateTimeField().to_representation(instance.created),
            'datetime': instance.dt,
            'active': instance.active,
        }
        return ret


class BlacklistedReactionsSerializer(BlacklistSerializer):
    """
    Serializer for a blacklisted reaction.
    """
    class Meta:
        model = BlacklistedReactions


class BlacklistedChemicalsSerializer(BlacklistSerializer):
    """
    Serializer for a blacklisted reaction.
    """
    class Meta:
        model = BlacklistedChemicals


class BlacklistViewSet(ModelViewSet):
    """
    Base class for API endpoint for accessing blacklisted items.

    Implements custom queryset to return objects belonging to current user
    and custom actions for activating and deactivating a blacklist entry.
    """
    lookup_field = 'id'
    permission_classes = [IsAuthenticated]

    def get_queryset(self):
        """
        Returns blacklisted items associated with the authenticated user.
        """
        return self.model.objects.filter(user=self.request.user)

    def get_serializer_context(self):
        """
        Add request information to context passed to serializer.
        """
        context = super().get_serializer_context()
        context.update({"request": self.request})
        return context

    def destroy(self, request, *args, **kwargs):
        """
        Override the default destroy method to return success status.
        """
        instance = self.get_object()
        data = self.get_serializer(instance).data  # get data before deletion
        self.perform_destroy(instance)
        return Response({'success': True, 'data': data})

    @action(detail=True, methods=['get'])
    def activate(self, request, *args, **kwargs):
        """
        For this blacklisted item, set active to True.
        """
        instance = self.get_object()
        instance.active = True
        instance.save()
        data = self.get_serializer(instance).data  # get data after update
        return Response({'success': True, 'data': data})

    @action(detail=True, methods=['get'])
    def deactivate(self, request, *args, **kwargs):
        """
        For this blacklisted item, set active to False.
        """
        instance = self.get_object()
        instance.active = False
        instance.save()
        data = self.get_serializer(instance).data  # get data after update
        return Response({'success': True, 'data': data})


class BlacklistedReactionsViewSet(BlacklistViewSet):
    """
    API endpoint for accessing blacklisted reactions.
    """
    model = BlacklistedReactions
    serializer_class = BlacklistedReactionsSerializer


class BlacklistedChemicalsViewSet(BlacklistViewSet):
    """
    API endpoint for accessing blacklisted chemicals.
    """
    model = BlacklistedChemicals
    serializer_class = BlacklistedChemicalsSerializer
