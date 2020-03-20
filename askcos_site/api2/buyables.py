import io
import json

import pandas as pd
from bson import ObjectId
from rdkit import Chem
from rest_framework import serializers
from rest_framework.authentication import SessionAuthentication
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.viewsets import ViewSet
from rest_framework_jwt.authentication import JSONWebTokenAuthentication

from askcos_site.globals import buyables_db
from askcos_site.main.views.users import can_modify_buyables


class BuyableSerializer(serializers.Serializer):
    """Serializer for buyable attributes"""
    smiles = serializers.CharField()
    ppg = serializers.FloatField(min_value=0.0)
    source = serializers.CharField(default='')
    allowOverwrite = serializers.BooleanField(required=False)  # no default so the field is empty when not specified

    def validate_smiles(self, value):
        """Verify that the smiles is valid. Returns canonicalized SMILES."""
        mol = Chem.MolFromSmiles(value)
        if not mol:
            raise serializers.ValidationError('Cannot parse target smiles with rdkit')
        return Chem.MolToSmiles(mol, isomericSmiles=True)


class BuyableUploadSerializer(serializers.Serializer):
    """Serializer for buyable upload"""
    file = serializers.FileField()
    format = serializers.CharField()
    returnLimit = serializers.IntegerField(default=1000)
    allowOverwrite = serializers.BooleanField(default=True)


class BuyableQuerySerializer(serializers.Serializer):
    """Serializer for a buyable query"""
    q = serializers.CharField(default='')
    source = serializers.CharField(default='')
    regex = serializers.BooleanField(default=False)
    limit = serializers.IntegerField(default=100)
    canonicalize = serializers.BooleanField(default=True)


class BuyablesViewSet(ViewSet):
    """
    A ViewSet for accessing buyables data.
    """

    authentication_classes = [SessionAuthentication, JSONWebTokenAuthentication]

    def list(self, request):
        """List all buyables in database matching query parameters"""
        serializer = BuyableQuerySerializer(data=request.query_params)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        search = data['q']
        source = data['source']
        regex = data['regex']
        limit = data['limit']
        canon = data['canonicalize']

        query = {}
        if search:
            if regex:
                query['smiles'] = {'$regex': '.*{}.*'.format(search)}
            else:
                if canon:
                    mol = Chem.MolFromSmiles(search)
                    if mol:
                        search = Chem.MolToSmiles(mol, isomericSmiles=True)
                query['smiles'] = search

        if source:
            query['source'] = source

        search_result = list(buyables_db.find(query, {'smiles': 1, 'ppg': 1, 'source': 1}).limit(limit))

        for doc in search_result:
            doc['_id'] = str(doc['_id'])

        resp = {
            'result': search_result,
            'search': search,
        }

        return Response(resp)

    def retrieve(self, request, pk=None):
        """Return single buyables entry by mongo _id"""
        resp = {'error': None, 'result': None}

        try:
            result = buyables_db.find_one({'_id': ObjectId(pk)})
        except Exception as e:
            resp['error'] = 'Could not retrieve item: {0!s}'.format(e)
        else:
            result['_id'] = str(result['_id'])
            resp['result'] = [result]

        return Response(resp)

    def destroy(self, request, pk=None):
        """Delete a specific buyables entry from the database"""
        resp = {'error': None, 'success': False}

        if not can_modify_buyables(request):
            resp['error'] = 'You are not authorized to modify the buyables database. Please ask your site administrator for permission.'
            return Response(resp, status=401)

        delete_result = buyables_db.delete_one({'_id': ObjectId(pk)})
        if delete_result.deleted_count != 1:
            resp['error'] = 'Could not find buyable'
        else:
            resp['success'] = True

        return Response(resp)

    def create(self, request):
        """Add new buyable entry"""
        resp = {'error': None, 'success': False}

        if not can_modify_buyables(request):
            resp['error'] = 'You are not authorized to modify the buyables database. Please ask your site administrator for permission.'
            return Response(resp, status=401)

        serializer = BuyableSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        smiles = data['smiles']
        ppg = data['ppg']
        source = data['source']
        allow_overwrite = data.get('allowOverwrite', True)

        result = self.add_buyable_to_db(
            {'smiles': smiles, 'ppg': ppg, 'source': source},
            allow_overwrite=allow_overwrite,
        )

        if result['error']:
            resp['error'] = result['error']
            return Response(resp)

        resp['success'] = True
        if result.get('duplicate'):
            resp['duplicate'] = True
        if result.get('inserted'):
            resp['inserted'] = result['inserted']
        if result.get('updated'):
            resp['updated'] = result['updated']

        return Response(resp)

    @action(detail=False, methods=['post'])
    def upload(self, request):
        """Upload a list of buyables"""
        resp = {'error': None, 'success': False}

        if not can_modify_buyables(request):
            resp['error'] = 'You are not authorized to modify the buyables database. Please ask your site administrator for permission.'
            return Response(resp, status=401)

        serializer = BuyableUploadSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data

        upload_file = data['file']
        file_format = data['format']
        return_limit = data['returnLimit']
        allow_overwrite = data['allowOverwrite']

        try:
            content = upload_file.read()
        except Exception as e:
            resp['error'] = 'Cannot read file: {0!s}'.format(e)
            return Response(resp)

        if file_format == 'json':
            try:
                upload_json = json.loads(content.decode('utf-8'))
            except json.JSONDecodeError:
                resp['error'] = 'Cannot parse json!'
                return Response(resp)
        elif file_format == 'csv':
            try:
                df = pd.read_csv(io.StringIO(content.decode('utf-8')))
                upload_json = df.to_dict(orient='records')
            except Exception as e:
                resp['error'] = 'Cannot parse csv: {0!s}'.format(e)
                return Response(resp)
        else:
            resp['error'] = 'File format not supported!'
            return Response(resp)

        if not isinstance(upload_json, list) or len(upload_json) == 0 or not isinstance(upload_json[0], dict):
            resp['error'] = 'Improperly formatted json!'
            return Response(resp)

        data_serializer = BuyableSerializer(data=upload_json, many=True)
        data_serializer.is_valid(raise_exception=True)
        buyables_data = data_serializer.validated_data

        result = self.add_buyable_list_to_db(buyables_data, allow_overwrite=allow_overwrite)

        if result.get('error'):
            resp['error'] = result['error']

        resp['success'] = True
        resp['added'] = result['added'][:return_limit]
        resp['updated'] = result['updated'][:return_limit]
        resp['added_count'] = result['added_count']
        resp['updated_count'] = result['updated_count']
        resp['duplicate_count'] = result['duplicate_count']
        resp['count'] = result['count']
        resp['total'] = result['total']

        return Response(resp)

    def add_buyable_list_to_db(self, buyable_list, allow_overwrite=True):
        """Add list of buyable compounds to the database"""
        result = {
            'error': None,
            'count': 0,
            'added': [],
            'updated': [],
            'added_count': 0,
            'updated_count': 0,
            'duplicate_count': 0,
            'total': len(buyable_list)
        }

        for buyable in buyable_list:
            res = self.add_buyable_to_db(buyable, allow_overwrite=allow_overwrite)
            if not res['error']:
                if res.get('duplicate'):
                    result['duplicate_count'] += 1
                if res.get('inserted'):
                    result['added'].append(res['inserted'])
                    result['added_count'] += 1
                if res.get('updated'):
                    result['updated'].append(res['updated'])
                    result['updated_count'] += 1
                result['count'] += 1
            else:
                result['error'] = res['error']

        return result

    def add_buyable_to_db(self, buyable, allow_overwrite=True):
        """Add a single buyable compound to the database"""
        result = {'error': None}

        smiles = buyable['smiles']
        ppg = buyable['ppg']
        source = buyable['source']

        new_doc = {
            'smiles': smiles,
            'ppg': ppg,
            'source': source
        }

        existing_doc = buyables_db.find_one({'smiles': smiles})
        if existing_doc and allow_overwrite:
            buyables_db.update_one(
                {'smiles': smiles},
                {'$set': {'ppg': ppg, 'source': source}}
            )
            new_doc['_id'] = str(existing_doc['_id'])
            result['updated'] = new_doc
        elif existing_doc and not allow_overwrite:
            result['duplicate'] = True
        else:
            insert_result = buyables_db.insert_one(new_doc)
            if not insert_result.inserted_id:
                result['error'] = 'Addition of buyable failed'
                return result
            new_doc['_id'] = str(new_doc['_id'])
            result['inserted'] = new_doc

        return result
