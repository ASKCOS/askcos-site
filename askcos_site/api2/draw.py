import re

from django.http import HttpResponse
from rest_framework import serializers
from rest_framework.generics import GenericAPIView
from rest_framework.response import Response

from askcos.utilities.io.draw import MolsSmilesToImage, MakeBackgroundTransparent
from askcos.utilities.io.draw import TransformStringToImage
from askcos.utilities.io.draw import ReactionStringToImage
from askcos.utilities.io.draw import MappedReactionToHightlightImage
from askcos.utilities.io.draw import MolsSmilesToImageHighlight


class DrawerSerializer(serializers.Serializer):
    """Serializer for drawing parameters."""
    smiles = serializers.CharField()
    input_type = serializers.CharField(required=False)
    transparent = serializers.BooleanField(default=False)
    draw_map = serializers.BooleanField(default=False)
    highlight = serializers.BooleanField(default=False)
    reacting_atoms = serializers.ListField(child=serializers.FloatField(), required=False)

    def validate_input_type(self, value):
        """Check that input type is accepted value."""
        if value not in ['chemical', 'reaction', 'template']:
            raise serializers.ValidationError("Valid input types: ['chemical', 'reaction', 'template']")


class DrawerAPIView(GenericAPIView):
    """
    API endpoint for drawing molecules, reactions, and templates.

    Notes:

    - Both GET and POST requests are possible. GET requests may be easier
      for simple linking, while POST requests are better for complex data.
    - If `input_type` is not specified, will attempt to determine type.
      Specifying `input_type` can provide faster results.

    Method: GET, POST

    Parameters:

    - `smiles` (str): input SMILES (or SMARTS) string
    - `input_type` (str, optional): one of 'chemical', 'reaction', or 'template'
    - `transparent` (bool, optional): whether background should be transparent (chemical only)
    - `draw_map` (bool, optional): whether atom mapping should be drawn (reaction only)
    - `highlight` (bool, optional): whether to highlight mapped atoms (reaction or chemical)
    - `reacting_atoms` (list, optional): list of atom scores to highlight and label (chemical only)

    Returns: PNG image of input SMILES
    """

    serializer_class = DrawerSerializer

    def get(self, request):
        """
        Handle GET request for drawing.
        """
        serializer = self.get_serializer(data=request.query_params)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        return draw(data)

    def post(self, request):
        """
        Handle POST request for drawing.
        """
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        data = serializer.validated_data
        return draw(data)


def draw(data):
    """
    Return HttpResponse with PNG of requested structure.
    """
    input_type = data.get('input_type')

    if input_type == 'chemical':
        methods = [draw_chemical]
    elif input_type == 'reaction':
        methods = [draw_reaction]
    elif input_type == 'template':
        methods = [draw_template]
    else:
        methods = [draw_chemical, draw_reaction, draw_template]

    for method in methods:
        try:
            response = method(data)
        except:
            continue
        else:
            return response
    else:
        return Response({'error': 'Could not draw requested structure.', 'request': data}, status=400)


def draw_chemical(data):
    """
    Returns HttpResponse containing PNG of chemical SMILES.
    """
    smiles = data.get('smiles')
    transparent = data.get('transparent')
    highlight = data.get('highlight')
    reacting_atoms = data.get('reacting_atoms', [])

    if highlight:
        res = MolsSmilesToImageHighlight(smiles, reacting_atoms=reacting_atoms, bonds=False, clear_map=True)
        if transparent:
            # <svg:rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300' height='300' x='0' y='0'> </svg:rect>
            # replace #FFFFFF with none
            res = re.sub(
                r"(.*)<svg:rect(.*)style='([^']*)'(.*)>(.*)</svg:rect>(.*)",
                r"\1<svg:rect\2style='opacity:0.0;fill:none;stroke:none'\4>\5</svg:rect>\6",
                res
            )
        response = HttpResponse(res, content_type='image/svg+xml')
    else:
        response = HttpResponse(content_type='image/png')
        if transparent:
            MakeBackgroundTransparent(MolsSmilesToImage(smiles)).save(response, 'png', quality=70)
        else:
            MolsSmilesToImage(smiles).save(response, 'png', quality=70)
    return response


def draw_template(data):
    """
    Returns HttpResponse containing PNG of reaction SMARTS.
    """
    template = data.get('smiles')
    response = HttpResponse(content_type='image/jpeg')
    TransformStringToImage(template).save(response, 'png', quality=70)
    return response


def draw_reaction(data):
    """
    Returns HttpResponse containing PNG of reaction SMILES.
    """
    smiles = data.get('smiles')
    strip = not data.get('draw_map')
    highlight = data.get('highlight')
    response = HttpResponse(content_type='image/jpeg')
    if highlight:
        MappedReactionToHightlightImage(smiles, highlightByReactant=True).save(response, 'png', quality=70)
    else:
        ReactionStringToImage(smiles, strip=strip).save(response, 'png', quality=70)
    return response


drawer = DrawerAPIView.as_view()
