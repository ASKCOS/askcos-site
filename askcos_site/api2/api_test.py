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

    def test_celery_status(self):
        """Test /celery endpoint"""
        response = self.client.get('https://localhost/api/v2/celery/')
        self.assertEqual(response.status_code, 200)

        result = response.json()
        self.assertIsInstance(result['queues'], list)

    def test_celery_task_status(self):
        """Test /celery/task endpoint"""
        response = self.client.get('https://localhost/api/v2/celery/task/')
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'task_id': ['This field is required.']})

        response = self.client.get('https://localhost/api/v2/celery/task/?task_id=abc')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), {'complete': False})

    def test_context(self):
        """Test /context endpoint"""
        data = {
            'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
            'products': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'num_results': 5,
            'return_scores': 'true',
        }
        response = self.client.post('https://localhost/api/v2/context/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['reactants'], data['reactants'])
        self.assertEqual(request['products'], data['products'])
        self.assertEqual(request['num_results'], data['num_results'])
        self.assertTrue(request['return_scores'])

        self.assertEqual(len(result['contexts']), 5)
        c = result['contexts'][0]
        self.assertEqual(c['catalyst'], '')
        self.assertEqual(c['reagent'], 'Cc1ccccc1.[H][N-][H].[Na+]')
        self.assertAlmostEqual(c['score'], 0.339, places=2)
        self.assertEqual(c['solvent'], '')
        self.assertAlmostEqual(c['temperature'], 94.48, places=1)

        # Test insufficient data
        response = self.client.post('https://localhost/api/v2/context/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['This field is required.'],
                                           'products': ['This field is required.']})

        # Test unparseable smiles
        response = self.client.post('https://localhost/api/v2/context/', data={'reactants': 'X', 'products': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['Cannot parse reactants smiles with rdkit.'],
                                           'products': ['Cannot parse products smiles with rdkit.']})

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

    def test_forward(self):
        """Test /forward endpoint"""
        data = {
            'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
            'num_results': 5,
        }
        response = self.client.post('https://localhost/api/v2/forward/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['reactants'], data['reactants'])
        self.assertEqual(request['reagents'], '')
        self.assertEqual(request['solvent'], '')
        self.assertEqual(request['num_results'], data['num_results'])

        self.assertEqual(len(result['outcomes']), 5)
        o = result['outcomes'][0]
        self.assertEqual(o['smiles'], 'CN(C)CCOC(c1ccccc1)c1ccccc1')
        self.assertAlmostEqual(o['mol_wt'], 255.36, places=2)
        self.assertEqual(o['rank'], 1)
        self.assertAlmostEqual(o['score'], -63.3, places=1)
        self.assertAlmostEqual(o['prob'], 0.91, places=2)

        # Test insufficient data
        response = self.client.post('https://localhost/api/v2/forward/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['This field is required.']})

        # Test unparseable smiles
        response = self.client.post('https://localhost/api/v2/forward/', data={'reactants': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['Cannot parse reactants smiles with rdkit.']})

    def test_impurity(self):
        """Test /impurity endpoint"""
        data = {
            'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
        }
        response = self.client.post('https://localhost/api/v2/impurity/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['reactants'], data['reactants'])
        self.assertEqual(request['reagents'], '')
        self.assertEqual(request['products'], '')
        self.assertEqual(request['solvent'], '')
        self.assertEqual(request['top_k'], 3)
        self.assertEqual(request['threshold'], 0.75)
        self.assertEqual(request['predictor'], 'WLN forward predictor')
        self.assertEqual(request['inspector'], 'Reaxys predictor')
        self.assertEqual(request['mapper'], 'WLN atom mapper')
        self.assertTrue(request['check_mapping'])

        self.assertIsNotNone(result['task_id'])

        # Test insufficient data
        response = self.client.post('https://localhost/api/v2/impurity/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['This field is required.']})

        # Test unparseable smiles
        response = self.client.post('https://localhost/api/v2/impurity/', data={'reactants': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['Cannot parse smiles with rdkit.']})

    def test_reactions(self):
        """Test /reactions endpoint"""
        data = {
            'ids': ['1'],
        }
        response = self.client.post('https://localhost/api/v2/reactions/', data=data)
        self.assertEqual(response.status_code, 200)

        result = response.json()
        self.assertEqual(result['reactions'], [])

        # Test insufficient data
        response = self.client.post('https://localhost/api/v2/reactions/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'ids': ['This field is required.']})

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

    def test_scscore(self):
        """Test /scscore endpoint"""
        data = {
            'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
        }
        response = self.client.post('https://localhost/api/v2/scscore/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['smiles'], data['smiles'])

        # Check that we got expected result
        self.assertAlmostEqual(result['score'], 2.159, places=3)

        # Test insufficient data
        response = self.client.post('https://localhost/api/v2/scscore/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['This field is required.']})

        # Test unparseable smiles
        response = self.client.post('https://localhost/api/v2/scscore/', data={'smiles': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['Cannot parse smiles with rdkit.']})

    def test_template(self):
        """Test /template endpoint"""
        # Get request for specific template endpoint
        template_id = '5e1f4b6f63488328509a594e'
        response = self.client.get('https://localhost/api/v2/template/{0}/'.format(template_id))
        self.assertEqual(response.status_code, 200)

        # Check retrieved template
        result = response.json()
        template = result['template']
        self.assertIsNotNone(template)
        self.assertEqual(template['template_set'], 'reaxys')
        self.assertEqual(template['_id'], template_id)
        self.assertEqual(template['references'], ['100336', '100337', '555364', '3948785', '3948834',
                                                  '28127174', '35585623', '38022824', '38022828',
                                                  '38022830', '38022834', '38022833', '38022835',
                                                  '38022845', '41610599', '41610601', '41610620'])

        # Get request for specific template reaxys query export endpoint
        response = self.client.get('https://localhost/api/v2/template/{0}/export'.format(template_id))
        self.assertEqual(response.status_code, 200)

        # Check exported json
        result = response.json()
        self.assertEqual(result['fileName'], 'reaxys_query.json')
        self.assertEqual(len(result['content']['facts']), 1)
        fact = result['content']['facts'][0]
        self.assertEqual(fact['id'], 'Reaxys487')
        self.assertEqual(len(fact['fields']), 1)
        field = fact['fields'][0]
        self.assertEqual(field['id'], 'RX.ID')
        self.assertEqual(field['displayName'], 'Reaction ID')
        self.assertEqual(field['boundOperator'], 'op_num_equal')
        self.assertEqual(field['value'],
                         '100336; 100337; 555364; 3948785; 3948834; 28127174; 35585623; 38022824; 38022828; 38022830; 38022834; 38022833; 38022835; 38022845; 41610599; 41610601; 41610620')

    def test_tree_builder(self):
        """Test /treebuilder endpoint"""
        data = {
            'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'async': False,
            'return_first': True,
        }
        response = self.client.post('https://localhost/api/v2/treebuilder/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['smiles'], data['smiles'])
        self.assertEqual(request['async'], False)
        self.assertEqual(request['return_first'], True)
        self.assertEqual(request['chemical_property_logic'], 'none')
        self.assertEqual(request['chemical_popularity_logic'], 'none')
        self.assertEqual(request['template_set'], 'reaxys')
        self.assertEqual(request['template_prioritizer'], 'reaxys')

        # Check that we got a result (can't check values because it's non-deterministic)
        self.assertIsInstance(result['trees'], list)
        self.assertIsInstance(result['trees'][0], dict)
        self.assertIsInstance(result['trees'][0]['children'], list)

        # Test insufficient data
        response = self.client.post('https://localhost/api/v2/treebuilder/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['This field is required.']})

        # Test unparseable smiles
        response = self.client.post('https://localhost/api/v2/treebuilder/', data={'smiles': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['Cannot parse smiles with rdkit.']})

    @classmethod
    def tearDownClass(cls):
        """This method is run once after all tests in this class."""
        cls.client.close()


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
