import re

from django.http import JsonResponse
from django.shortcuts import render, HttpResponse
from django.urls import reverse
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from askcos_site.main.utils import ajax_error_wrapper, resolve_smiles


@ajax_error_wrapper
def ajax_smiles_to_image(request):
    '''Takes an Ajax call with a smiles string
    and returns the HTML for embedding an image'''

    smiles = request.GET.get('smiles', None)
    print('SMILES from Ajax: {}'.format(smiles))
    smiles = resolve_smiles(smiles)
    if smiles is None:
        return JsonResponse({'err': True})
    print('Resolved smiles -> {}'.format(smiles))

    url = reverse('draw_smiles', kwargs={'smiles': smiles})
    data = {
        'html': '<img src="' + url + '">',
        'smiles': smiles,
    }
    return JsonResponse(data)


@ajax_error_wrapper
def ajax_rxn_to_image(request):
    '''Takes an Ajax call with a rxn smiles string
    and returns the HTML for embedding an image'''

    reactants = request.GET.get('reactants', '')
    product = request.GET.get('product', '')
    strip = request.GET.get('strip', True)

    reactants = resolve_smiles(reactants)
    product = resolve_smiles(product)
    smiles = reactants + '>>' + product

    if request.method == 'GET' and 'rxnsmiles' in request.GET:
        tmp = request.GET['rxnsmiles']
        if tmp is not None and tmp != '':
            smiles = tmp
    print(smiles)
    print('RXN SMILES from Ajax: {}'.format(smiles))
    url = reverse('draw_reaction', kwargs={'smiles': smiles})
    data = {
        'html': '<img src="' + url + '">',
        'smiles': smiles,
        'reactants': reactants,
        'product': product,
    }
    return JsonResponse(data)


# @login_required
def draw_smiles(request, smiles):
    '''
    Returns a png response for a target smiles
    '''
    from askcos.utilities.io.draw import MolsSmilesToImage, MakeBackgroundTransparent
    isTransparent = request.GET.get('transparent', 'False')
    response = HttpResponse(content_type='image/png')
    if isTransparent.lower() in ['true', 't', 'yes', 'y', '1']:
        MakeBackgroundTransparent(MolsSmilesToImage(str(smiles))).save(response, 'png', quality=70)
    else:
        MolsSmilesToImage(str(smiles)).save(response, 'png', quality=70)
    return response


#@login_required
def draw_template(request, template):
    '''
    Returns a png response for a reaction SMARTS template
    '''
    from askcos.utilities.io.draw import TransformStringToImage
    response = HttpResponse(content_type='image/jpeg')
    TransformStringToImage(str(template)).save(response, 'png', quality=70)
    return response


# @login_required
def draw_reaction(request, smiles):
    '''
    Returns a png response for a SMILES reaction string
    '''
    from askcos.utilities.io.draw import ReactionStringToImage
    response = HttpResponse(content_type='image/jpeg')
    ReactionStringToImage(str(smiles)).save(response, 'png', quality=70)
    return response

def draw_mapped_reaction(request, smiles):
    '''
    Returns a png response for a SMILES reaction string
    '''
    from askcos.utilities.io.draw import ReactionStringToImage
    response = HttpResponse(content_type='image/jpeg')
    print('in views', smiles)
    ReactionStringToImage(str(smiles), strip=False).save(response, 'png', quality=70)
    return response

def draw_highlighted_reaction(request, smiles):
    '''
        Returns a png response for a SMILES reaction string
    '''
    from askcos.utilities.io.draw import MappedReactionToHightlightImage
    response = HttpResponse(content_type='image/jpeg')
    print('in views', smiles)
    MappedReactionToHightlightImage(str(smiles), highlightByReactant=True).save(response, 'png', quality=70)
    return response

def draw_smiles_highlight(request, smiles, reacting_atoms, bonds='False'):
    '''
    Returns a svg xml with atoms highlighted
    '''
    from askcos.utilities.io.draw import MolsSmilesToImageHighlight
    from ast import literal_eval
    reacting_atoms = literal_eval(reacting_atoms)
    #TODO has to be a better way to evaluate string to true or false
    bonds = bonds.lower() in ['true', '1', 't', 'y', 'yes']
    res = MolsSmilesToImageHighlight(smiles, reacting_atoms=reacting_atoms, bonds=bonds, clear_map=True)
    IsTransparent = request.GET.get('transparent', 'False')
    if IsTransparent.lower() in ['true', '1', 't', 'y', 'yes']:
        # <svg:rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300' height='300' x='0' y='0'> </svg:rect>
        # replace #FFFFFF with none
        res = re.sub(r"(.*)<svg:rect(.*)style='([^']*)'(.*)>(.*)</svg:rect>(.*)",
               r"\1<svg:rect\2style='opacity:0.0;fill:none;stroke:none'\4>\5</svg:rect>\6", res)
    response = HttpResponse(res, content_type='image/svg+xml')
  
    return response


def draw(request):
    return render(request, 'draw.html')


def draw_fig(request, fig):
    '''
    Returns a png response for figure object
    '''
    response = HttpResponse(content_type='img/png')
    canvas = FigureCanvas(fig)
    canvas.print_png(response)
    return response
