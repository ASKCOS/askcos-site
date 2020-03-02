'''
Template-free forward predictor
'''

from celery import shared_task
from celery.signals import celeryd_init
import numpy as np
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
from rdkit import Chem
from rdkit.Chem import Descriptors

tffp = None
CORRESPONDING_QUEUE = 'tffp_worker'

@celeryd_init.connect
def configure_worker(options={}, **kwargs):
    print(options)
    if 'queues' not in options:
        return
    if CORRESPONDING_QUEUE not in options['queues'].split(','):
        return
    print('### STARTING UP A TEMPLATE-FREE FORWARD PREDICTOR WORKER ###')

    global tffp

    # Import as needed
    from makeit.synthetic.evaluation.rexgen_direct.predict import TFFP
    print('Imported TFFP')
    try:
        tffp = TFFP()
    except Exception as e:
        print(e)
        raise(e)
    print('Initialized')
    tffp.predict('CCCO.CCCBr')

    print('Finished configuring TFFP worker')


@shared_task
def get_outcomes(reactants, top_n=10):
    global tffp
    mapped_smiles, results = tffp.predict(reactants, top_n=top_n)
    results_to_return = {}
    original_reactants = reactants.split('.')
    for res in results:
        smiles_list = set(res['smiles'].split('.'))
        smiles_canonical = set()
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if not mol:
                continue
            smiles_canonical.add(Chem.MolToSmiles(mol))
            # Remove unreacted frags
        smiles_canonical = smiles_canonical - set(original_reactants)
        if not smiles_canonical:
            continue # no reaction?
        smiles = max(smiles_canonical, key=len) # NOTE: this is not great...byproducts may be longer
        if not smiles:
            continue
        if smiles in results_to_return:
            results_to_return[smiles]['rank'] = min(results_to_return[smiles]['rank'], res['rank'])
            results_to_return[smiles]['score'] = np.log(np.exp(results_to_return[smiles]['score']) + np.exp(res['score']))
            results_to_return[smiles]['prob'] += res['prob']
        else:
            # Append outcome information
            results_to_return[smiles] = {
                'rank': res['rank'],
                'smiles': smiles,
                'score': float(res['score']),
                'prob': float(res['prob']),
                'mol_wt': float(Descriptors.MolWt(Chem.MolFromSmiles(smiles)))
            }
    results_to_return = sorted(results_to_return.values(), key=lambda x: x['prob'], reverse=True)
    total_prob = sum([outcome['prob'] for outcome in results_to_return])
    for i, outcome in enumerate(results_to_return):
        results_to_return[i]['rank'] = i + 1
        results_to_return[i]['prob'] = outcome['prob'] / total_prob
    return results_to_return
