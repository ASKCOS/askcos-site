import json

from django.http import JsonResponse
from rdkit import Chem


def canonicalize(request):
    resp = {}
    if request.method != 'POST':
        resp['error'] = 'Must use POST request'
        return JsonResponse(resp, status=400)
    if not request.body:
        resp['error'] = 'Did not recieve request body'
        return JsonResponse(resp, status=400)
    original_smiles = json.loads(request.body.decode('utf-8')).get('smiles')
    if not original_smiles:
        resp['error'] = 'SMILES was not sent or body could not be parsed'
        return JsonResponse(resp, status=400)
    mol = Chem.MolFromSmiles(original_smiles)
    if not mol:
        resp['error'] = 'Cannot parse smiles with rdkit'
        return JsonResponse(resp, status=400)
    smiles = Chem.MolToSmiles(mol)
    if not smiles:
        resp['error'] = 'Cannot canonicalize smiles with rdkit'
        return JsonResponse(resp, status=400)
    resp['smiles'] = smiles
    return JsonResponse(resp)

def molfile_to_smiles(request):
    resp = {}
    molfile = request.GET.get('molfile', '')
    isomericSmiles = request.GET.get('isomericSmiles', 'True') in ['True', 'true']
    mol = Chem.MolFromMolBlock(molfile)
    if not mol:
        resp['error'] = 'Cannot parse sdf molfile with rdkit'
        return JsonResponse(resp, status=400)
    try:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles)
    except:
        resp['error'] = 'Cannot parse sdf molfile with rdkit'
        return JsonResponse(resp, status=400)
    resp['smiles'] = smiles
    return JsonResponse(resp)

def smiles_to_molfile(request):
    resp = {}
    smiles = request.GET.get('smiles', '')
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        resp['error'] = 'Cannot parse smiles with rdkit'
        return JsonResponse(resp, status=400)
    try:
        molfile = Chem.MolToMolBlock(mol)
    except:
        resp['error'] = 'Cannot parse smiles with rdkit'
        return JsonResponse(resp, status=400)
    resp['molfile'] = molfile
    return JsonResponse(resp)
