import requests


class TorchserveAPI:
    """
    Base class for torchserve API queries.
    """

    def __init__(self, hostname, model_name, *args, version=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.inference_url = 'http://{0}:8080/predictions/{1}'.format(hostname, model_name)
        if version is not None:
            self.inference_url += '/{0}'.format(version)
        self.management_url = 'http://{0}:8081'.format(hostname)

    def inference(self, data):
        """
        Perform inference query.
        """
        resp = requests.post(self.inference_url, json=data)
        if resp.status_code == 200:
            return resp.json()
        else:
            return {'error': 'Received {0}: {1}'.format(resp.status_code, resp.reason)}
