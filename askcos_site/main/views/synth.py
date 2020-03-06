import time

from django.http import JsonResponse
from django.shortcuts import render
from django.template.loader import render_to_string

from askcos_site.askcos_celery.treeevaluator.scoring_coordinator import evaluate
from askcos_site.main.utils import ajax_error_wrapper
from makeit.utilities.contexts import clean_context


def synth_interactive(request):
    context = {}
    return render(request, 'forward_interactive.html', context)

#@login_required
def synth_interactive_smiles(request, smiles):
    '''Synth interactive initialized w/ reaction smiles'''
    return synth_interactive(request, reactants=smiles.split('>')[0], product=smiles.split('>')[-1])

@ajax_error_wrapper
def ajax_validate_temperature_synth(request):
    '''Checks to see if a temperature is valid'''
    data = {'err': False}

    temperature = request.GET.get('temperature', None)
    print('temperature from Ajax: {}'.format(temperature))

    try:
        temperature = float(temperature)
        data['temperature'] = temperature
    except Exception as e:
        data['err'] = True

    return JsonResponse(data)

@ajax_error_wrapper
def ajax_start_synth(request):
    '''Perform the forward synthesis'''
    data = {'err': False}

    reactants = request.GET.get('reactants', None)
    solvent = request.GET.get('solvent', None)
    temperature = request.GET.get('temperature', None)
    reagents = request.GET.get('reagents', None)
    mincount = int(request.GET.get('mincount', 25))
    maxreturn = int(request.GET.get('maxreturn', 100))
    forward_scorer = request.GET.get('forward_scorer', 'Template_Free')
    print('Conditions for forward synthesis:')
    print('reactants: {}'.format(reactants))
    print('solvent: {}'.format(solvent))
    print('temp: {}'.format(temperature))
    print('reagents: {}'.format(reagents))
    print('mincount: {}'.format(mincount))
    print('max return: {}'.format(maxreturn))
    print('forward scorer: {}'.format(forward_scorer))

    startTime = time.time()


    # context expected is (T1, slvt1, rgt1, cat1, t1, y1)
    if solvent == 'default':
        solvent = ''
        print('reset default solvent')

    res = evaluate.delay(reactants, '',
        contexts=[clean_context((temperature, solvent, reagents, '', -1, -1))],
        forward_scorer=forward_scorer, top_n=maxreturn, return_all_outcomes=True)
    outcomes = res.get(300)[0]['outcomes']

    print('Got top outcomes, length {}'.format(len(outcomes)))
    #print(outcomes)
    data['html_time'] = '{:.3f} seconds elapsed'.format(time.time() - startTime)

    if outcomes:

        data['html'] = render_to_string('synth_outcomes_only.html',
            {'outcomes': outcomes})
    else:
        data['html'] = 'No outcomes found? That is weird...'

    # Save in session in case used wants to print
    request.session['last_synth_interactive'] = {'reactants': reactants,
        'temperature': temperature, 'reagents': reagents, 'solvent': solvent,
        'mincount': mincount, 'outcomes': outcomes, 'forward_scorer': forward_scorer}

    return JsonResponse(data)
