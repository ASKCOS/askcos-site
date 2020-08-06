import numpy as np
from askcos.utilities.fingerprinting import create_rxn_Morgan2FP_separately
from ..tfx_api_model import TFXAPIModel

class TFXFastFilter(TFXAPIModel):
    """Fast filter Tensorflow API Model. Overrides input and output transformation methods with fast filter specific methods.

    Attributes:
        hostname (str): hostname of service serving tf model.
        model_name (str): Name of model provided to tf serving.
        version (int): version of the model to use when serving
    """
    def transform_input(self, reactant_smiles, target, rxnfpsize=2048, pfpsize=2048, useFeatures=False):
        """Transforms the input for the API model from SMILES strings to product and reaction fingerprints

        Args:
            reactant_smiles (str): SMILES string of the reactants
            target (str): SMILES string of the target
            rxnfpsize (int): Length of desired fingerprint for the reaction. Must agree with model input shape
            pfpsize (int): Length of desired fingerprint for the product. Must agree with model input shape
            useFeatures (bool): Flag to use features or not when generating fingerprint. Should agree with how model was trained

        Returns:
            list of dict: Input fingerprints, formatted for a call to the tensorflow API model
        """
        pfp, rfp = create_rxn_Morgan2FP_separately(
            reactant_smiles, target, rxnfpsize=rxnfpsize, pfpsize=pfpsize, useFeatures=useFeatures
        )
        pfp = np.asarray(pfp, dtype='float32')
        rfp = np.asarray(rfp, dtype='float32')
        rxnfp = pfp - rfp
        return [{
            'input_1': pfp.tolist(),
            'input_2': rxnfp.tolist()
        }]
    
    def transform_output(self, pred):
        """Transforms the output of the prediction into the fast filter score

        Args:
            pred (np.array): numpy array of the output of the prediction

        Returns
            float: the first element of the prediction is returned
        """
        return pred[0]