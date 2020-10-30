import requests


class TorchserveAPI:
    """
    Base class for Torchserve API queries.
    """

    def __init__(self, hostname, model_name, *args, version=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.prediction_url = 'http://{0}:8080/predictions/{1}'.format(hostname, model_name)
        if version is not None:
            self.prediction_url += '/{0}'.format(version)
        self.management_url = 'http://{0}:8081'.format(hostname)

    def predict(self, data):
        """
        Perform model prediction by sending request to Torchserve API.
        """
        resp = requests.post(self.prediction_url, json=data)
        if resp.status_code == 200:
            return resp.json()
        else:
            return {'error': 'Received {0}: {1}'.format(resp.status_code, resp.reason)}
