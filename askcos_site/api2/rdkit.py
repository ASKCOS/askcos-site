from rdkit import Chem
from rest_framework import serializers
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.viewsets import ViewSet


class SmilesSerializer(serializers.Serializer):
    """Serializer for SMILES strings."""
    smiles = serializers.CharField()
    isomericSmiles = serializers.BooleanField(default=True)


class MolfileSerializer(serializers.Serializer):
    """Serializer for Molfiles."""
    molfile = serializers.CharField(trim_whitespace=False)
    isomericSmiles = serializers.BooleanField(default=True)


class SmilesViewSet(ViewSet):
    """
    A ViewSet for exposing various SMILES functionality.
    """

    @action(detail=False, methods=['POST'])
    def canonicalize(self, request):
        """
        Canonicalize the specified SMILES using RDKit.

        Method: POST

        Parameters:

        - `smiles` (str): SMILES string to canonicalize
        - `isomericSmiles` (bool, optional): whether to generate isomeric SMILES

        Returns:

        - `smiles` (str): canonicalized SMILES
        """
        serializer = SmilesSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        smiles = serializer.validated_data['smiles']
        isomericSmiles = serializer.validated_data['isomericSmiles']

        resp = {}

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            resp['error'] = 'Cannot parse smiles with rdkit.'
            return Response(resp, status=400)

        smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)
        if not smiles:
            resp['error'] = 'Cannot canonicalize smiles with rdkit.'
            return Response(resp, status=400)

        resp['smiles'] = smiles

        return Response(resp)

    @action(detail=False, methods=['POST'])
    def validate(self, request):
        """
        Check the syntax and validity of a SMILES string.

        Method: POST

        Parameters:

        - `smiles` (str): SMILES string

        Returns:

        - `correct_syntax` (bool): correctness of the SMILES syntax
        - `valid_chem_name` (bool): validity of the chemical

        RDKit ref: https://github.com/rdkit/rdkit/issues/2430
        """
        serializer = SmilesSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        smiles = serializer.validated_data['smiles']

        mol = Chem.MolFromSmiles(smiles, sanitize=False)

        if mol is None:
            correct_syntax = False
            valid_chem_name = False
        else:
            correct_syntax = True
            try:
                Chem.SanitizeMol(mol)
            except:
                valid_chem_name = False
            else:
                valid_chem_name = True

        resp = {
            'correct_syntax': correct_syntax,
            'valid_chem_name': valid_chem_name,
        }

        return Response(resp)

    @action(detail=False, methods=['POST'])
    def from_molfile(self, request):
        """
        Convert the provided Molfile to a SMILES string.

        Method: POST

        Parameters:

        - `molfile` (str): Molfile input
        - `isomericSmiles` (bool, optional): whether to generate isomeric SMILES

        Returns:

        - `smiles` (str): canonical SMILES
        """
        serializer = MolfileSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        molfile = serializer.validated_data['molfile']
        isomericSmiles = serializer.validated_data['isomericSmiles']

        resp = {}

        mol = Chem.MolFromMolBlock(molfile)
        if not mol:
            resp['error'] = 'Cannot parse sdf molfile with rdkit.'
            return Response(resp, status=400)

        try:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)
        except:
            resp['error'] = 'Cannot parse sdf molfile with rdkit.'
            return Response(resp, status=400)

        resp['smiles'] = smiles

        return Response(resp)

    @action(detail=False, methods=['POST'])
    def to_molfile(self, request):
        """
        Convert the provided SMILES string to a Molfile.

        Method: POST

        Parameters:

        - `smiles` (str): SMILES input

        Returns:

        - `molfile` (str): Molfile output
        """
        serializer = SmilesSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        smiles = serializer.validated_data['smiles']

        resp = {}

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            resp['error'] = 'Cannot parse smiles with rdkit.'
            return Response(resp, status=400)

        try:
            molfile = Chem.MolToMolBlock(mol)
        except:
            resp['error'] = 'Cannot parse smiles with rdkit.'
            return Response(resp, status=400)

        resp['molfile'] = molfile

        return Response(resp)
