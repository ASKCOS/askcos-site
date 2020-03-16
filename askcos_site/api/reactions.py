import json
from django.http import JsonResponse
import makeit.global_config as gc
from pymongo import MongoClient

client = MongoClient(
    gc.MONGO['path'],
    gc.MONGO['id'],
    connect=gc.MONGO['connect']
)
db_name = gc.REACTIONS['database']
collection = gc.REACTIONS['collection']
reactions_collection = client[db_name][collection]

def reactions(request):
    resp = {}
    if request.method != 'POST':
        resp['error'] = 'Must use POST request'
        return JsonResponse(resp, status=400)
    if not request.body:
        resp['error'] = 'Did not recieve request body'
        return JsonResponse(resp, status=400)
    body = json.loads(request.body.decode('utf-8'))
    _ids = body.get('ids')
    template_set = body.get('template_set')
    query = {'reaction_id': {'$in': _ids}}
    if template_set:
        query['template_set'] = template_set
    reactions_by_ids = list(reactions_collection.find(query))
    resp['reactions'] = reactions_by_ids
    return JsonResponse(resp)
