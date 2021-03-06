from bson.objectid import ObjectId
from django.http import JsonResponse

from askcos_site.globals import retro_templates


def template_sets(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    resp['template_sets'] = retro_templates.distinct('template_set')
    return JsonResponse(resp)


def template(request):
    resp = {}
    resp['request'] = dict(**request.GET)
    _id = request.GET.get('id')
    transform = retro_templates.find_one({'_id': _id})
    if not transform:
        try:
            transform = retro_templates.find_one({'_id': ObjectId(_id)})
        except:
            transform = None
    if not transform:
        resp['error'] = 'Cannot find template with id '+_id
        return JsonResponse(resp)
    transform['_id'] = _id
    transform.pop('product_smiles', None)
    transform.pop('name', None)
    if transform.get('template_set') == 'reaxys':
        refs = transform.pop('references', [''])
        transform['references'] = list(map(lambda x: x.split('-')[0], refs))
    resp['template'] = transform
    return JsonResponse(resp)


def reaxys_export(request):
    resp = {}
    _id = request.GET.get('id')
    transform = retro_templates.find_one({'_id': _id})
    if not transform:
        try:
            transform = retro_templates.find_one({'_id': ObjectId(_id)})
        except:
            transform = None
    if not transform:
        resp['error'] = 'Cannot find template with id '+_id
        return JsonResponse(resp)
    if transform.get('template_set') != 'reaxys':
        resp['error'] = 'Template is not in the reaxys template set'
        return JsonResponse(resp)
    references = '; '.join(map(lambda ref: ref.split('-')[0], transform['references']))
    resp['fileName'] = 'reaxys_query.json'
    resp['version'] = '1.0'
    resp['content'] = {
        'id': 'root',
        'facts': [{
            'id': 'Reaxys487',
            'fields': [{
                'value': references,
                'boundOperator': 'op_num_equal',
                'id': 'RX.ID',
                'displayName': 'Reaction ID'
            }],
            'fieldsLogicOperator': 'AND',
            'exist': False,
            'bio': False
        }]
    }
    resp['exist'] = False
    resp['bio'] = False
    resp['content']['facts'][0].pop('logicOperator', None)
    return JsonResponse(resp)
