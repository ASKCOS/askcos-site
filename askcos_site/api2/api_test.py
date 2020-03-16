"""
Live test suite for askcos web api

This will test api functionality for a live instance of the askcos site
"""

import unittest

from requests import Session


class TestAPI(unittest.TestCase):
    """Test class for askcos api v2"""

    @classmethod
    def setUpClass(cls):
        """This method is run once before all tests in this class."""
        cls.client = Session()
        cls.client.verify = False

    def test_buyables(self):
        """Test /buyables endpoint"""
        # Get request for main endpoint
        response = self.client.get('https://localhost/api/v2/buyables/')
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(len(result['result']), 100)  # returns 100 results by default

        # Post request to add buyable
        data = {
            'smiles': 'C1CCC1',
            'ppg': '2.0',
            'allowOverwrite': False,
        }
        response = self.client.post('https://localhost/api/v2/buyables/', data=data)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])
        _id = result['inserted']['_id']

        # Post request to upload buyables (both are duplicates)
        filedata = '[{"smiles": "C1CCC1","ppg": "2.0"},{"smiles": "C1CCCC1","ppg": "3.0"}]'
        files = {'file': ('upload.json', filedata)}
        data = {'format': 'json', 'allowOverwrite': False}
        response = self.client.post('https://localhost/api/v2/buyables/upload/', data=data, files=files)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])
        self.assertEqual(result['added'], [])
        self.assertEqual(result['updated'], [])
        self.assertEqual(result['duplicate_count'], 2)
        self.assertEqual(result['total'], 2)

        # Get request with query
        response = self.client.get('https://localhost/api/v2/buyables/?q=C1CCC1')
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(len(result['result']), 1)
        self.assertEqual(result['result'][0]['smiles'], 'C1CCC1')

        # Get request for specific buyable
        response = self.client.get('https://localhost/api/v2/buyables/{0}/'.format(_id))
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(len(result['result']), 1)
        self.assertEqual(result['result'][0]['smiles'], 'C1CCC1')

        # Delete request for specific buyable
        response = self.client.delete('https://localhost/api/v2/buyables/{0}/'.format(_id))
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])

    def test_fast_filter(self):
        """Test /fast-filter endpoint"""
        data = {
            'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
            'products': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
        }
        response = self.client.post('https://localhost/api/v2/fast-filter/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['reactants'], data['reactants'])
        self.assertEqual(request['products'], data['products'])

        # Check result
        self.assertAlmostEqual(result['score'], 0.998, places=2)

        # Test insufficient data
        response = self.client.post('https://localhost/api/v2/fast-filter/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['This field is required.'],
                                           'products': ['This field is required.']})

        # Test unparseable smiles
        response = self.client.post('https://localhost/api/v2/fast-filter/', data={'reactants': 'X', 'products': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['Cannot parse reactants smiles with rdkit.'],
                                           'products': ['Cannot parse products smiles with rdkit.']})

    def test_retro(self):
        """Test /retro endpoint"""
        data = {
            'target': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
        }
        response = self.client.post('https://localhost/api/v2/retro/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['target'], data['target'])
        self.assertEqual(request['num_templates'], 100)
        self.assertEqual(request['max_cum_prob'], 0.995)
        self.assertEqual(request['filter_threshold'], 0.75)
        self.assertEqual(request['template_set'], 'reaxys')
        self.assertEqual(request['template_prioritizer'], 'reaxys')
        self.assertTrue(request['cluster'])
        self.assertEqual(request['cluster_method'], 'kmeans')
        self.assertEqual(request['cluster_feature'], 'original')
        self.assertEqual(request['cluster_fp_type'], 'morgan')
        self.assertEqual(request['cluster_fp_length'], 512)
        self.assertEqual(request['cluster_fp_radius'], 1)

        # Test insufficient data
        response = self.client.post('https://localhost/api/v2/retro/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'target': ['This field is required.']})

        # Test unparseable smiles
        response = self.client.post('https://localhost/api/v2/retro/', data={'target': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'target': ['Cannot parse target smiles with rdkit.']})

    @classmethod
    def tearDownClass(cls):
        """This method is run once after all tests in this class."""
        cls.client.close()


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
