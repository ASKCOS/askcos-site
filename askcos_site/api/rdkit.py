from rdkit import Chem
from django.http import JsonResponse

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
