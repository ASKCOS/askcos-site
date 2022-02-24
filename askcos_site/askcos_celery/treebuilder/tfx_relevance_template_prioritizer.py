import numpy as np
import requests
import rdkit.Chem as Chem
from askcos.prioritization.templates.relevance import RelevanceTemplatePrioritizer
from rdkit.Chem import AllChem
from scipy.special import softmax

from ..tfx_api_model import TFXAPIModel


class TFXRelevanceTemplatePrioritizer(RelevanceTemplatePrioritizer):
    def __init__(self, model=None, hostname=None, model_name=None, version=None):
        self.model = model or TFXRelevanceTemplateAPIModel(hostname, model_name, version)
        self.fp_length = self.model.get_input_dim()
        self.fp_radius = 2


class TFXRelevanceTemplateAPIModel(TFXAPIModel):
    """Template relevance Tensorflow API Model. 
    
    Overrides input and output transformation methods with template relevance 
    specific methods.

    Attributes:
        hostname (str): hostname of service serving tf model.
        model_name (str): Name of model provided to tf serving.
        version (int): version of the model to use when serving
    """
    def __init__(self, hostname, model_name, version, **kwargs):
        super().__init__(hostname, model_name, version)
        self.metaurl = self.baseurl+'/metadata'
        self.fp_length = self.get_input_dim()

    def get_input_dim(self):
        resp = requests.get(self.metaurl)
        metadata = resp.json()['metadata']['signature_def']['signature_def']['serving_default']
        input_dim = int(list(metadata['inputs'].values())[0]['tensor_shape']['dim'][1]['size'])
        return input_dim

    def transform_input(self, fp, **kwargs):
        """Transforms the input for the API model to a list

        Args:
            fp (np.array): Fingerprint bit vector as numpy array

        Returns:
            list: Fingerprint bit vector as a list
        """
        return fp.tolist()
