from django.utils import timezone
from rdkit import Chem
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


class BlacklistedReactionsSerializer(BlacklistSerializer):
    """
    Serializer for a blacklisted reaction.
    """
    class Meta:
        model = BlacklistedReactions

    def validate_smiles(self, value):
        """
        Verify that the provided SMILES is valid. Returns canonicalized SMILES.
        """
        try:
            reactants, agents, products = value.split('>')
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
        return reactants + '>' + agents + '>' + products


class BlacklistedChemicalsSerializer(BlacklistSerializer):
    """
    Serializer for a blacklisted reaction.
    """
    class Meta:
        model = BlacklistedChemicals

    def validate_smiles(self, value):
        """
        Verify that the provided SMILES is valid. Returns canonicalized SMILES.
        """
        if value:
            mol = Chem.MolFromSmiles(value)
            if not mol:
                raise serializers.ValidationError('Cannot parse smiles with rdkit.')
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        return value


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
        API endpoint to activate specified blacklist entry.

        Method: GET

        Returns:

        - `success`: true if successfully activated
        - `data`: updated entry data
        """
        instance = self.get_object()
        instance.active = True
        instance.save()
        data = self.get_serializer(instance).data  # get data after update
        return Response({'success': True, 'data': data})

    @action(detail=True, methods=['get'])
    def deactivate(self, request, *args, **kwargs):
        """
        API endpoint to deactivate specified blacklist entry.

        Method: GET

        Returns:

        - `success`: true if successfully deactivated
        - `data`: updated entry data
        """
        instance = self.get_object()
        instance.active = False
        instance.save()
        data = self.get_serializer(instance).data  # get data after update
        return Response({'success': True, 'data': data})


class BlacklistedReactionsViewSet(BlacklistViewSet):
    """
    API endpoint for accessing blacklisted reactions.

    Method: GET

    Returns: list of blacklisted reaction entries belonging to the currently authenticated user

    Method: POST

    Parameters:

    - `smiles` (str): reaction SMILES string to be added
    - `description` (str, optional): text description
    - `datetime` (str, optional): timestamp, default current time
    - `active` (bool, optional): whether this entry is active, default true

    Returns: created entry

    ----------
    For a particular entry, specified as URI parameter (`/api/v2/blacklist/reactions/<id>/`):

    Method: GET

    Returns: entry with requested id

    Method: DELETE

    Returns:

    - `success`: true if successfully deleted
    - `data`: data from deleted entry
    """
    model = BlacklistedReactions
    serializer_class = BlacklistedReactionsSerializer


class BlacklistedChemicalsViewSet(BlacklistViewSet):
    """
    API endpoint for accessing blacklisted chemicals.

    Method: GET

    Returns: list of blacklisted chemical entries belonging to the currently authenticated user

    Method: POST

    Parameters:

    - `smiles` (str): chemical SMILES string to be added
    - `description` (str, optional): text description
    - `datetime` (str, optional): timestamp, default current time
    - `active` (bool, optional): whether this entry is active, default true

    Returns: created entry

    ----------
    For a particular entry, specified as URI parameter (`/api/v2/blacklist/chemicals/<id>/`):

    Method: GET

    Returns: entry with requested id

    Method: DELETE

    Returns:

    - `success`: true if successfully deleted
    - `data`: data from deleted entry
    """
    model = BlacklistedChemicals
    serializer_class = BlacklistedChemicalsSerializer
