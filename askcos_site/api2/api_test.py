"""
Live test suite for askcos web api

This will test api functionality for a live instance of the askcos site
"""

import time
import unittest

from requests import Session

username = ''
password = ''


class TestAPI(unittest.TestCase):
    """Test class for askcos api v2"""

    @classmethod
    def setUpClass(cls):
        """This method is run once before all tests in this class."""
        cls.client = Session()
        cls.client.verify = False

        cls.url = 'https://localhost/api/v2'

    def get(self, endpoint, **kwargs):
        """Process a GET request"""
        return self.client.get(self.url + endpoint, **kwargs)

    def post(self, endpoint, **kwargs):
        """Process a POST request"""
        return self.client.post(self.url + endpoint, **kwargs)

    def delete(self, endpoint, **kwargs):
        """Process a DELETE request"""
        return self.client.delete(self.url + endpoint, **kwargs)

    def get_result(self, task_id):
        """Retrieve celery task output"""
        # Try to get result 10 times in 2 sec intervals
        for _ in range(10):
            response = self.get('/celery/task/{0}/'.format(task_id))
            result = response.json()
            if result.get('complete'):
                return result
            else:
                if result.get('failed'):
                    self.fail('Celery task failed.')
                else:
                    time.sleep(2)

    def authenticate(self):
        """Get authentication token, returns formatted header dict"""
        response = self.post('/token-auth/', data={'username': username, 'password': password})
        self.assertEqual(response.status_code, 200)

        result = response.json()
        self.assertIn('token', result)
        token = result['token']

        return {'Authorization': 'Bearer {0}'.format(token)}

    def test_atom_mapper(self):
        """Test /atom-mapper endpoint"""
        data = {
            'rxnsmiles': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1>>CN(C)CCOC(c1ccccc1)c1ccccc1',
        }
        response = self.post('/atom-mapper/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['rxnsmiles'], data['rxnsmiles'])
        self.assertEqual(request['mapper'], 'WLN atom mapper')
        self.assertEqual(request['priority'], 1)

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])
        self.assertEqual(result['output'], '[CH3:1][N:2]([CH3:3])[CH2:4][CH2:5][Cl:6].[OH:7][CH:8]([c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[c:15]1[cH:16][cH:17][cH:18][cH:19][cH:20]1>>[CH3:1][N:2]([CH3:3])[CH2:4][CH2:5][O:7][CH:8]([c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[c:15]1[cH:16][cH:17][cH:18][cH:19][cH:20]1')

        # Test insufficient data
        response = self.post('/atom-mapper/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'rxnsmiles': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/atom-mapper/', data={'rxnsmiles': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'rxnsmiles': ['Cannot parse reaction smiles.']})

        response = self.post('/atom-mapper/', data={'rxnsmiles': 'X>>Y'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'rxnsmiles': ['Cannot parse reactants using rdkit.']})

        # Test task priority argument
        data = {
            'rxnsmiles': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1>>CN(C)CCOC(c1ccccc1)c1ccccc1',
            'priority': 2,
        }
        response = self.post('/atom-mapper/', data=data)
        self.assertEqual(response.status_code, 200)

        result = response.json()
        request = result['request']
        self.assertEqual(request['priority'], 2)

    @unittest.skipIf(not (username and password), 'Requires login credentials.')
    def test_banlist_chemicals(self):
        """Test /banlist/chemicals endpoint"""
        response = self.get('/banlist/chemicals/')
        self.assertEqual(response.status_code, 401)

        headers = self.authenticate()

        # Post request to add banned chemical
        data = {
            'smiles': 'c1ccccc1',
            'description': 'test',
        }
        response = self.post('/banlist/chemicals/', data=data, headers=headers)
        self.assertEqual(response.status_code, 201)
        result = response.json()
        self.assertEqual(result['smiles'], data['smiles'])
        self.assertEqual(result['description'], data['description'])
        self.assertTrue(result['active'])
        self.assertIn('created', result)
        created_entry = result

        # Get list of banned chemicals
        response = self.get('/banlist/chemicals/', headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertIsInstance(result, list)
        entries = [x for x in result if x['smiles'] == 'c1ccccc1' and x['description'] == 'test']
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0], created_entry)
        entry_id = entries[0]['id']

        # Get detail view of specific entry
        response = self.get('/banlist/chemicals/{0}/'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(result, created_entry)

        # Deactivate entry
        response = self.get('/banlist/chemicals/{0}/deactivate'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])
        self.assertFalse(result['data']['active'])

        # Activate entry
        response = self.get('/banlist/chemicals/{0}/activate'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])
        self.assertTrue(result['data']['active'])

        # Delete entry
        response = self.delete('/banlist/chemicals/{0}/'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])

        # Try to retrieve entry again
        response = self.get('/banlist/chemicals/{0}/'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 404)

        # List entries again
        response = self.get('/banlist/chemicals/', headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertIsInstance(result, list)
        entries = [x for x in result if x['smiles'] == 'c1ccccc1' and x['description'] == 'test']
        self.assertEqual(entries, [])

    @unittest.skipIf(not (username and password), 'Requires login credentials.')
    def test_banlist_reactions(self):
        """Test /banlist/reactions endpoint"""
        response = self.get('/banlist/reactions/')
        self.assertEqual(response.status_code, 401)

        headers = self.authenticate()

        # Post request to add banned chemical
        data = {
            'smiles': 'c1ccccc1>>Brc1ccccc1',
            'description': 'test',
        }
        response = self.post('/banlist/reactions/', data=data, headers=headers)
        self.assertEqual(response.status_code, 201)
        result = response.json()
        self.assertEqual(result['smiles'], data['smiles'])
        self.assertEqual(result['description'], data['description'])
        self.assertTrue(result['active'])
        self.assertIn('created', result)
        created_entry = result

        # Get list of banned chemicals
        response = self.get('/banlist/reactions/', headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertIsInstance(result, list)
        entries = [x for x in result if x['smiles'] == 'c1ccccc1>>Brc1ccccc1' and x['description'] == 'test']
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0], created_entry)
        entry_id = entries[0]['id']

        # Get detail view of specific entry
        response = self.get('/banlist/reactions/{0}/'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(result, created_entry)

        # Deactivate entry
        response = self.get('/banlist/reactions/{0}/deactivate'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])
        self.assertFalse(result['data']['active'])

        # Activate entry
        response = self.get('/banlist/reactions/{0}/activate'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])
        self.assertTrue(result['data']['active'])

        # Delete entry
        response = self.delete('/banlist/reactions/{0}/'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])

        # Try to retrieve entry again
        response = self.get('/banlist/reactions/{0}/'.format(entry_id), headers=headers)
        self.assertEqual(response.status_code, 404)

        # List entries again
        response = self.get('/banlist/reactions/', headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertIsInstance(result, list)
        entries = [x for x in result if x['smiles'] == 'c1ccccc1>>Brc1ccccc1' and x['description'] == 'test']
        self.assertEqual(entries, [])

    def test_buyables(self):
        """Test /buyables endpoint"""
        # Get request for main endpoint
        response = self.get('/buyables/')
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(len(result['result']), 100)  # returns 100 results by default

        # Post request to add buyable
        data = {
            'smiles': 'C1CC2(C1)CCC2',
            'ppg': '2.0',
            'source': 'test',
            'allowOverwrite': False,
        }
        response = self.post('/buyables/', data=data)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])
        _id = result['inserted']['_id']

        # Post request to upload buyables (duplicate entry)
        filedata = '[{"smiles": "C1CC2(C1)CCC2","ppg": "2.0","source": "test"}]'
        files = {'file': ('upload.json', filedata)}
        data = {'format': 'json', 'allowOverwrite': False}
        response = self.post('/buyables/upload/', data=data, files=files)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])
        self.assertEqual(result['inserted'], [])
        self.assertEqual(result['updated'], [])
        self.assertEqual(result['duplicate_count'], 1)
        self.assertEqual(result['total'], 1)

        # Get request with query
        response = self.get('/buyables/?q=C1CC2(C1)CCC2')
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(len(result['result']), 1)
        self.assertEqual(result['result'][0]['smiles'], 'C1CC2(C1)CCC2')

        # Get request with source query
        response = self.get('/buyables/?source=test')
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(len(result['result']), 1)
        self.assertEqual(result['result'][0]['smiles'], 'C1CC2(C1)CCC2')

        # Get request for specific buyable
        response = self.get('/buyables/{0}/'.format(_id))
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(len(result['result']), 1)
        self.assertEqual(result['result'][0]['smiles'], 'C1CC2(C1)CCC2')

        # Delete request for specific buyable
        response = self.delete('/buyables/{0}/'.format(_id))
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertTrue(result['success'])

    def test_celery_status(self):
        """Test /celery endpoint"""
        response = self.get('/celery/')
        self.assertEqual(response.status_code, 200)

        result = response.json()
        self.assertIsInstance(result['queues'], list)

    def test_celery_task_status(self):
        """Test /celery/task endpoint"""
        response = self.get('/celery/task/abc/')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), {'complete': False})

    def test_cluster(self):
        """Test /cluster endpoint"""
        data = {
            'original': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'outcomes': ['OC(c1ccccc1)c1ccccc1.CN(C)CCCl',
                         'ClC(c1ccccc1)c1ccccc1.CN(C)CCO',
                         'O=C(c1ccccc1)c1ccccc1.CN(C)CCO']
        }
        response = self.post('/cluster/', data=data)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(result['request']['original'], data['original'])
        self.assertEqual(result['request']['outcomes'], data['outcomes'])

        self.assertIsInstance(result['output'], list)
        self.assertIn(result['output'], [[0, 0, 1], [1, 1, 0]])

        response = self.post('/cluster/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'original': ['This field is required.'],
                                           'outcomes': ['This field is required.']})

    def test_context(self):
        """Test /context endpoint"""
        data = {
            'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
            'products': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'num_results': 5,
            'return_scores': True,
        }
        response = self.post('/context/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['reactants'], data['reactants'])
        self.assertEqual(request['products'], data['products'])
        self.assertEqual(request['num_results'], data['num_results'])
        self.assertTrue(request['return_scores'])

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])
        self.assertEqual(len(result['output']), 5)
        o = result['output'][0]
        self.assertEqual(o['catalyst'], '')
        self.assertEqual(o['reagent'], 'Cc1ccccc1.[H][N-][H].[Na+]')
        self.assertAlmostEqual(o['score'], 0.339, places=2)
        self.assertEqual(o['solvent'], '')
        self.assertAlmostEqual(o['temperature'], 94.48, places=1)

        # Test insufficient data
        response = self.post('/context/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['This field is required.'],
                                           'products': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/context/', data={'reactants': 'X', 'products': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['Cannot parse reactants smiles with rdkit.'],
                                           'products': ['Cannot parse products smiles with rdkit.']})

    def test_drawing(self):
        """Test /draw endpoint"""
        expected = {'smiles': ['This field is required.']}
        response = self.get('/draw/')
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), expected)
        response = self.post('/draw/')
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), expected)

        data = {
            'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'input_type': 'invalid',
        }
        expected = {'input_type': ["Valid input types: ['chemical', 'reaction', 'template']"]}
        response = self.get('/draw/', params=data)
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), expected)
        response = self.post('/draw/', data=data)
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), expected)

        tests = [
            # Test chemical smiles
            {
                'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            },
            # Test transparent chemical smiles
            {
                'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
                'transparent': True,
            },
            # Test reaction smiles
            {
                'smiles': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1>>CN(C)CCOC(c1ccccc1)c1ccccc1',
            },
            # Test mapped reaction smiles
            {
                'smiles': '[CH3:1][N:2]([CH3:3])[CH2:4][CH2:5][Cl:6].[OH:7][CH:8]([c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[c:15]1[cH:16][cH:17][cH:18][cH:19][cH:20]1>>[CH3:1][N:2]([CH3:3])[CH2:4][CH2:5][O:7][CH:8]([c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[c:15]1[cH:16][cH:17][cH:18][cH:19][cH:20]1',
                'draw_mapped': True,
            },
            # Test highlighted reaction smiles
            {
                'smiles': '[CH3:1][N:2]([CH3:3])[CH2:4][CH2:5][Cl:6].[OH:7][CH:8]([c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[c:15]1[cH:16][cH:17][cH:18][cH:19][cH:20]1>>[CH3:1][N:2]([CH3:3])[CH2:4][CH2:5][O:7][CH:8]([c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[c:15]1[cH:16][cH:17][cH:18][cH:19][cH:20]1',
                'highlight': True,
            },
            # Test template smarts
            {
                'smiles': '[#7:1]-[C:2](=[O;D1;H0:3])-[CH;@;D3;+0:4]1-[CH2;D2;+0:5]-[CH;@;D3;+0:9](-[C:7](-[#8:6])=[O;D1;H0:8])-[NH;D2;+0:10]-[CH;@;D3;+0:11]-1-[c:12]>>[#7:1]-[C:2](=[O;D1;H0:3])-[CH;D2;+0:4]=[CH2;D1;+0:5].[#8:6]-[C:7](=[O;D1;H0:8])-[CH2;D2;+0:9]/[N;H0;D2;+0:10]=[CH;D2;+0:11]/[c:12]',
            },
            # Test site selectivity drawing
            {
                'smiles': 'Cc1ccccc1',
                'highlight': True,
                'reacting_atoms': [
                    1.4901161193847656e-07,
                    2.9802322387695312e-08,
                    0.32867664098739624,
                    0.0002092421054840088,
                    0.9734209775924683,
                    0.00020927190780639648,
                    0.32867667078971863
                ]
            },
        ]

        for data in tests:
            response1 = self.get('/draw/', params=data)
            self.assertEqual(response1.status_code, 200)
            response2 = self.post('/draw/', data=data)
            self.assertEqual(response2.status_code, 200)
            self.assertEqual(response1.content, response2.content)

    def test_fast_filter(self):
        """Test /fast-filter endpoint"""
        data = {
            'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
            'products': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
        }
        response = self.post('/fast-filter/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['reactants'], data['reactants'])
        self.assertEqual(request['products'], data['products'])

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])
        self.assertAlmostEqual(result['output'], 0.998, places=2)

        # Test insufficient data
        response = self.post('/fast-filter/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['This field is required.'],
                                           'products': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/fast-filter/', data={'reactants': 'X', 'products': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['Cannot parse reactants smiles with rdkit.'],
                                           'products': ['Cannot parse products smiles with rdkit.']})

    def test_forward(self):
        """Test /forward endpoint"""
        data = {
            'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
            'num_results': 5,
        }
        response = self.post('/forward/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['reactants'], data['reactants'])
        self.assertEqual(request['reagents'], '')
        self.assertEqual(request['solvent'], '')
        self.assertEqual(request['num_results'], data['num_results'])
        self.assertEqual(request['priority'], 1)

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])
        self.assertEqual(len(result['output']), 5)
        o = result['output'][0]
        self.assertEqual(o['smiles'], 'CN(C)CCOC(c1ccccc1)c1ccccc1')
        self.assertAlmostEqual(o['mol_wt'], 255.36, places=2)
        self.assertEqual(o['rank'], 1)
        self.assertAlmostEqual(o['score'], -63.3, places=1)
        self.assertAlmostEqual(o['prob'], 0.91, places=2)

        # Test insufficient data
        response = self.post('/forward/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/forward/', data={'reactants': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['Cannot parse reactants smiles with rdkit.']})

        # Test task priority argument
        data = {
            'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
            'num_results': 5,
            'priority': 2,
        }
        response = self.post('/forward/', data=data)
        self.assertEqual(response.status_code, 200)

        result = response.json()
        request = result['request']
        self.assertEqual(request['priority'], 2)

    def test_impurity(self):
        """Test /impurity endpoint"""
        data = {
            'reactants': 'CN(C)CCCl.OC(c1ccccc1)c1ccccc1',
        }
        response = self.post('/impurity/', data=data)
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
        self.assertEqual(request['inspector'], 'Reaxys inspector')
        self.assertEqual(request['mapper'], 'WLN atom mapper')
        self.assertTrue(request['check_mapping'])

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Test insufficient data
        response = self.post('/impurity/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/impurity/', data={'reactants': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'reactants': ['Cannot parse smiles with rdkit.']})

    def test_path_ranking(self):
        """Test /path-ranking endpoint"""
        data = {
            'trees': """[{"smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1", "ppg": 6.0, "as_reactant": 0, "as_product": 0, "id": 4, "is_chemical": true, "children": [{"plausibility": 0.9988686442375183, "template_score": 0.07315370440483093, "tforms": ["5e1f4b6e6348832850997243"], "num_examples": 185, "necessary_reagent": "", "id": 10, "is_reaction": true, "children": [{"smiles": "BrC(c1ccccc1)c1ccccc1", "ppg": 1.0, "as_reactant": 0, "as_product": 0, "id": 9, "is_chemical": true, "children": []}, {"smiles": "CN(C)CCO", "ppg": 1.0, "as_reactant": 0, "as_product": 0, "id": 5, "is_chemical": true, "children": []}], "smiles": "BrC(c1ccccc1)c1ccccc1.CN(C)CCO>>CN(C)CCOC(c1ccccc1)c1ccccc1"}]}, {"smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1", "ppg": 6.0, "as_reactant": 0, "as_product": 0, "id": 4, "is_chemical": true, "children": [{"plausibility": 0.9775225520133972, "template_score": 0.17757552862167358, "tforms": ["5e1f4b6e6348832850996b90"], "num_examples": 266, "necessary_reagent": "", "id": 6, "is_reaction": true, "children": [{"smiles": "CN(C)CCO", "ppg": 1.0, "as_reactant": 0, "as_product": 0, "id": 5, "is_chemical": true, "children": []}, {"smiles": "OC(c1ccccc1)c1ccccc1", "ppg": 1.0, "as_reactant": 0, "as_product": 0, "id": 2, "is_chemical": true, "children": []}], "smiles": "CN(C)CCO.OC(c1ccccc1)c1ccccc1>>CN(C)CCOC(c1ccccc1)c1ccccc1"}]}, {"smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1", "ppg": 6.0, "as_reactant": 0, "as_product": 0, "id": 4, "is_chemical": true, "children": [{"plausibility": 0.9868380427360535, "template_score": 0.016103610396385193, "tforms": ["5e1f4b6e634883285099626f"], "num_examples": 697, "necessary_reagent": "", "id": 18, "is_reaction": true, "children": [{"smiles": "BrCCOC(c1ccccc1)c1ccccc1", "ppg": 0.0, "as_reactant": 0, "as_product": 0, "id": 16, "is_chemical": true, "children": [{"plausibility": 0.915161669254303, "template_score": 0.07843249291181564, "tforms": ["5e1f4b6e6348832850995fbf", "5e1f4b6e6348832850995f00"], "num_examples": 3287, "necessary_reagent": "", "id": 20, "is_reaction": true, "children": [{"smiles": "BrCCBr", "ppg": 1.0, "as_reactant": 0, "as_product": 0, "id": 19, "is_chemical": true, "children": []}, {"smiles": "OC(c1ccccc1)c1ccccc1", "ppg": 1.0, "as_reactant": 0, "as_product": 0, "id": 2, "is_chemical": true, "children": []}], "smiles": "BrCCBr.OC(c1ccccc1)c1ccccc1>>BrCCOC(c1ccccc1)c1ccccc1"}]}, {"smiles": "CNC", "ppg": 13.0, "as_reactant": 0, "as_product": 0, "id": 17, "is_chemical": true, "children": []}], "smiles": "BrCCOC(c1ccccc1)c1ccccc1.CNC>>CN(C)CCOC(c1ccccc1)c1ccccc1"}]}, {"smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1", "ppg": 6.0, "as_reactant": 0, "as_product": 0, "id": 4, "is_chemical": true, "children": [{"plausibility": 0.9411429762840271, "template_score": 0.010665317997336388, "tforms": ["5e1f4b6e6348832850995d9c"], "num_examples": 13475, "necessary_reagent": "", "id": 53, "is_reaction": true, "children": [{"smiles": "CN(C)C(=O)COC(c1ccccc1)c1ccccc1", "ppg": 0.0, "as_reactant": 0, "as_product": 0, "id": 52, "is_chemical": true, "children": [{"plausibility": 0.9692589640617371, "template_score": 0.00962438341230154, "tforms": ["5e1f4b6e6348832850995e68"], "num_examples": 3231, "necessary_reagent": "", "id": 61, "is_reaction": true, "children": [{"smiles": "CNC", "ppg": 13.0, "as_reactant": 0, "as_product": 0, "id": 17, "is_chemical": true, "children": []}, {"smiles": "O=C(Cl)COC(c1ccccc1)c1ccccc1", "ppg": 0.0, "as_reactant": 0, "as_product": 0, "id": 60, "is_chemical": true, "children": [{"plausibility": 0.9996992349624634, "template_score": 0.8562169075012207, "tforms": ["5e1f4b6e6348832850995d80"], "num_examples": 26695, "necessary_reagent": "[Cl]", "id": 59, "is_reaction": true, "children": [{"smiles": "O=C(O)COC(c1ccccc1)c1ccccc1", "ppg": 0.0, "as_reactant": 0, "as_product": 0, "id": 58, "is_chemical": true, "children": [{"plausibility": 0.9207471013069153, "template_score": 0.053130947053432465, "tforms": ["5e1f4b6e6348832850995d76"], "num_examples": 112035, "necessary_reagent": "", "id": 71, "is_reaction": true, "children": [{"smiles": "COC(=O)COC(c1ccccc1)c1ccccc1", "ppg": 0.0, "as_reactant": 0, "as_product": 0, "id": 70, "is_chemical": true, "children": [{"plausibility": 0.9830392599105835, "template_score": 0.1137804239988327, "tforms": ["5e1f4b6e634883285099602f"], "num_examples": 1201, "necessary_reagent": "", "id": 69, "is_reaction": true, "children": [{"smiles": "COC(=O)CBr", "ppg": 1.0, "as_reactant": 0, "as_product": 0, "id": 68, "is_chemical": true, "children": []}, {"smiles": "OC(c1ccccc1)c1ccccc1", "ppg": 1.0, "as_reactant": 0, "as_product": 0, "id": 2, "is_chemical": true, "children": []}], "smiles": "COC(=O)CBr.OC(c1ccccc1)c1ccccc1>>COC(=O)COC(c1ccccc1)c1ccccc1"}]}], "smiles": "COC(=O)COC(c1ccccc1)c1ccccc1>>O=C(O)COC(c1ccccc1)c1ccccc1"}]}], "smiles": "O=C(O)COC(c1ccccc1)c1ccccc1>>O=C(Cl)COC(c1ccccc1)c1ccccc1"}]}], "smiles": "CNC.O=C(Cl)COC(c1ccccc1)c1ccccc1>>CN(C)C(=O)COC(c1ccccc1)c1ccccc1"}]}], "smiles": "CN(C)C(=O)COC(c1ccccc1)c1ccccc1>>CN(C)CCOC(c1ccccc1)c1ccccc1"}]}, {"smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1", "ppg": 6.0, "as_reactant": 0, "as_product": 0, "id": 4, "is_chemical": true, "children": [{"plausibility": 0.9868380427360535, "template_score": 0.016103610396385193, "tforms": ["5e1f4b6e634883285099626f"], "num_examples": 697, "necessary_reagent": "", "id": 18, "is_reaction": true, "children": [{"smiles": "BrCCOC(c1ccccc1)c1ccccc1", "ppg": 0.0, "as_reactant": 0, "as_product": 0, "id": 16, "is_chemical": true, "children": [{"plausibility": 0.9061155915260315, "template_score": 0.014876967296004295, "tforms": ["5e1f4b6e6348832850996936"], "num_examples": 319, "necessary_reagent": "", "id": 27, "is_reaction": true, "children": [{"smiles": "BrCCI", "ppg": 0.0, "as_reactant": 0, "as_product": 0, "id": 26, "is_chemical": true, "children": [{"plausibility": 0.9610458612442017, "template_score": 0.05479388311505318, "tforms": ["5e1f4b6e6348832850996027"], "num_examples": 1208, "necessary_reagent": "[I]", "id": 33, "is_reaction": true, "children": [{"smiles": "CS(=O)(=O)OCCBr", "ppg": 0.0, "as_reactant": 0, "as_product": 0, "id": 32, "is_chemical": true, "children": [{"plausibility": 0.9341546893119812, "template_score": 0.011639130301773548, "tforms": ["5e1f4b6e634883285099633a"], "num_examples": 615, "necessary_reagent": "[Br]", "id": 36, "is_reaction": true, "children": [{"smiles": "CS(=O)(=O)OCCOS(C)(=O)=O", "ppg": 21.0, "as_reactant": 0, "as_product": 0, "id": 35, "is_chemical": true, "children": []}], "smiles": "CS(=O)(=O)OCCOS(C)(=O)=O>>CS(=O)(=O)OCCBr"}]}], "smiles": "CS(=O)(=O)OCCBr>>BrCCI"}]}, {"smiles": "OC(c1ccccc1)c1ccccc1", "ppg": 1.0, "as_reactant": 0, "as_product": 0, "id": 2, "is_chemical": true, "children": []}], "smiles": "BrCCI.OC(c1ccccc1)c1ccccc1>>BrCCOC(c1ccccc1)c1ccccc1"}]}, {"smiles": "CNC", "ppg": 13.0, "as_reactant": 0, "as_product": 0, "id": 17, "is_chemical": true, "children": []}], "smiles": "BrCCOC(c1ccccc1)c1ccccc1.CNC>>CN(C)CCOC(c1ccccc1)c1ccccc1"}]}]""",
        }
        response = self.post('/path-ranking/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertTrue(request['cluster'])
        self.assertEqual(request['cluster_method'], 'hdbscan')
        self.assertEqual(request['min_samples'], 5)
        self.assertEqual(request['min_cluster_size'], 5)

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])

        output = result['output']
        self.assertIn('scores', output)
        self.assertEqual(len(output['scores']), 5)
        self.assertEqual(output['scores'][0], -1)
        self.assertEqual(output['scores'][1], -1)

        self.assertIn('encoded_trees', output)
        self.assertEqual(len(output['encoded_trees']), 5)
        self.assertEqual(len(output['encoded_trees'][0]), 0)
        self.assertEqual(len(output['encoded_trees'][1]), 0)
        self.assertEqual(len(output['encoded_trees'][2]), 512)
        self.assertEqual(len(output['encoded_trees'][3]), 512)
        self.assertEqual(len(output['encoded_trees'][4]), 512)

        self.assertIn('clusters', output)
        self.assertEqual(output['clusters'], [-1, -1, 0, 1, 2])

        # Test insufficient data
        response = self.post('/path-ranking/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'trees': ['This field is required.']})

    def test_reactions(self):
        """Test /reactions endpoint"""
        data = {
            'ids': ['1'],
        }
        response = self.post('/reactions/', data=data)
        self.assertEqual(response.status_code, 200)

        result = response.json()
        self.assertEqual(result['reactions'], [])

        # Test insufficient data
        response = self.post('/reactions/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'ids': ['This field is required.']})

    def test_rdkit_smiles(self):
        """Test /rdkit/smiles endpoints"""
        data = {
            'smiles': 'c1ccccc1C(OCCN(C)C)c1ccccc1',
        }

        # Test canonicalization
        response = self.post('/rdkit/smiles/canonicalize/', data=data)
        self.assertEqual(response.status_code, 200)

        # Check that we got expected result
        self.assertEqual(response.json(), {'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1'})

        # Test insufficient data
        response = self.post('/rdkit/smiles/canonicalize/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['This field is required.']})

        # Test validation
        response = self.post('/rdkit/smiles/validate/', data=data)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), {'correct_syntax': True, 'valid_chem_name': True})

        molfile = {
            'molfile': """
     RDKit          2D

 19 20  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    3.8971    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    5.1962    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    6.4952    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    7.7942    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    9.0933    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    7.7942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.2500    3.8971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7500    3.8971    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  1  0
  8  9  1  0
  9 10  1  0
 10 11  1  0
 11 12  1  0
 11 13  1  0
  7 14  1  0
 14 15  2  0
 15 16  1  0
 16 17  2  0
 17 18  1  0
 18 19  2  0
  6  1  1  0
 19 14  1  0
M  END
""",
        }

        # Test molfile generation
        response = self.post('/rdkit/smiles/to_molfile/', data=data)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), molfile)

        # Test molfile parsing
        response = self.post('/rdkit/smiles/from_molfile/', data=molfile)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), {'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1'})

    @unittest.skipIf(not (username and password), 'Requires login credentials.')
    def test_results(self):
        """Test /results endpoint and token authentication"""
        # Test that access is denied without authentication
        response = self.get('/results/')
        self.assertEqual(response.status_code, 401)

        # Test that we can access using the token
        headers = self.authenticate()
        response = self.get('/results/', headers=headers)
        self.assertEqual(response.status_code, 200)

        result = response.json()
        self.assertIsInstance(result['results'], list)
        result_id = result['results'][0]['id']

        # Test checking the status of a result
        response = self.get('/results/{0}/check'.format(result_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(result['id'], result_id)
        self.assertEqual(result['state'], 'completed')
        self.assertIsNone(result['error'])

        # Test retrieving a result
        response = self.get('/results/{0}/'.format(result_id), headers=headers)
        self.assertEqual(response.status_code, 200)
        result = response.json()
        self.assertEqual(result['id'], result_id)
        self.assertIsInstance(result['result'], dict)
        self.assertIsNone(result['error'])

        # Test deleting a non-existent result
        response = self.delete('/results/{0}/'.format('random'), headers=headers)
        self.assertEqual(response.status_code, 404)
        result = response.json()
        self.assertFalse(result['success'])
        self.assertEqual(result['error'], 'Result not found!')

    def test_retro(self):
        """Test /retro endpoint"""
        data = {
            'target': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
        }
        response = self.post('/retro/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['target'], data['target'])
        self.assertEqual(request['num_templates'], 100)
        self.assertEqual(request['max_cum_prob'], 0.995)
        self.assertEqual(request['filter_threshold'], 0.75)
        self.assertEqual(request['template_set'], 'reaxys')
        self.assertEqual(request['template_prioritizer_version'], 0)
        self.assertTrue(request['cluster'])
        self.assertEqual(request['cluster_method'], 'kmeans')
        self.assertEqual(request['cluster_feature'], 'original')
        self.assertEqual(request['cluster_fp_type'], 'morgan')
        self.assertEqual(request['cluster_fp_length'], 512)
        self.assertEqual(request['cluster_fp_radius'], 1)
        self.assertEqual(request['priority'], 1)

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])
        self.assertIsInstance(result['output'], list)

        # Test insufficient data
        response = self.post('/retro/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'target': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/retro/', data={'target': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'target': ['Cannot parse target smiles with rdkit.']})

        # Test task priority argument
        data = {
            'target': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'priority': 2,
        }
        response = self.post('/retro/', data=data)
        self.assertEqual(response.status_code, 200)

        result = response.json()
        request = result['request']
        self.assertEqual(request['priority'], 2)

    def test_retro_models(self):
        """Test /retro/models endpoint"""
        data = {'template_set': 'reaxys'}
        response = self.get('/retro/models/', params=data)
        self.assertEqual(response.status_code, 200)

        result = response.json()
        self.assertEqual(result['request'], data)
        self.assertEqual(result['versions'], ['1'])

        # Test insufficient data
        response = self.get('/retro/models/')
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'template_set': ['This field is required.']})

    def test_root(self):
        """Test / endpoint"""
        response = self.get('/')
        self.assertEqual(response.status_code, 200)

    def test_scscore(self):
        """Test /scscore endpoint"""
        data = {
            'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
        }
        response = self.post('/scscore/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['smiles'], data['smiles'])

        # Check that we got expected result
        self.assertAlmostEqual(result['score'], 2.159, places=3)

        # Test insufficient data
        response = self.post('/scscore/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/scscore/', data={'smiles': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['Cannot parse smiles with rdkit.']})

    def test_selectivity(self):
        """Test /selectivity endpoint"""
        data = {
            'smiles': 'Cc1ccccc1',
        }
        response = self.post('/selectivity/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['smiles'], data['smiles'])

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])
        self.assertIsInstance(result['output'], list)
        self.assertEqual(len(result['output']), 123)

        # Test insufficient data
        response = self.post('/selectivity/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/selectivity/', data={'smiles': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['Cannot parse smiles with rdkit.']})

    def test_selectivity_gen(self):
        """Test /general-selectivity endpoint"""
        data = {
            'rxnsmiles': '[Br:1][Br:2].[NH2:3][c:4]1[n:5][cH:6][n:7][c:8]2[nH:9][cH:10][n:11][c:12]12>O>[Br:2][c:10]1[nH:9][c:8]2[n:7][cH:6][n:5][c:4]([NH2:3])[c:12]2[n:11]1.[Br:2][c:6]1[n:5][c:4]([NH2:3])[c:12]2[c:8]([n:7]1)[nH:9][cH:10][n:11]2',
        }
        response = self.post('/general-selectivity/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['rxnsmiles'], data['rxnsmiles'])

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])
        self.assertIsInstance(result['output'], list)
        self.assertEqual(len(result['output']), 2)
        self.assertAlmostEqual(result['output'][0], 0.99444, places=4)
        self.assertAlmostEqual(result['output'][1], 0.00555, places=4)

        # Test insufficient data
        response = self.post('/general-selectivity/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'rxnsmiles': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/general-selectivity/', data={'rxnsmiles': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'rxnsmiles': ['Cannot parse reaction smiles.']})

    def test_template(self):
        """Test /template endpoint"""
        # Get request for specific template endpoint
        template_id = '5e1f4b6f63488328509a594e'
        response = self.get('/template/{0}/'.format(template_id))
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
        response = self.get('/template/{0}/export'.format(template_id))
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

    def test_template_sets(self):
        """Test /template/sets endpoint"""
        response = self.get('/template/sets/')
        self.assertEqual(response.status_code, 200)

        result = response.json()
        self.assertIn('template_sets', result)
        self.assertIn('reaxys', result['template_sets'])

    def test_tree_builder(self):
        """Test /tree-builder endpoint"""
        data = {
            'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'buyable_logic': 'or',
            'return_first': True,
        }
        response = self.post('/tree-builder/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['smiles'], data['smiles'])
        self.assertEqual(request['return_first'], True)
        self.assertEqual(request['chemical_property_logic'], 'none')
        self.assertEqual(request['chemical_popularity_logic'], 'none')
        self.assertEqual(request['template_set'], 'reaxys')
        self.assertEqual(request['template_prioritizer_version'], 0)
        self.assertEqual(request['priority'], 1)

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])
        self.assertIsInstance(result['output'], list)
        self.assertIsInstance(result['output'][0], dict)
        self.assertIsInstance(result['output'][0]['children'], list)

        # Test insufficient data
        response = self.post('/tree-builder/', data={})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['This field is required.']})

        # Test unparseable smiles
        response = self.post('/tree-builder/', data={'smiles': 'X'})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {'smiles': ['Cannot parse smiles with rdkit.']})

        # Test error when store_results=True without authentication
        data = {
            'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'store_results': True,
        }
        response = self.client.post('https://localhost/api/v2/tree-builder/', data=data)
        self.assertEqual(response.status_code, 401)
        result = response.json()
        self.assertEqual(result['error'], 'You must be authenticated to store tree builder results.')

        # Test task priority argument
        data = {
            'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'buyable_logic': 'or',
            'return_first': True,
            'priority': 2,
        }
        response = self.post('/tree-builder/', data=data)
        self.assertEqual(response.status_code, 200)

        result = response.json()
        request = result['request']
        self.assertEqual(request['priority'], 2)

    def test_tree_builder_v2(self):
        """Test /tree-builder endpoint for tree-builder v2"""
        data = {
            'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'buyable_logic': 'or',
            'return_first': True,
            'version': 2,
        }
        response = self.post('/tree-builder/', data=data)
        self.assertEqual(response.status_code, 200)

        # Confirm that request was interpreted correctly
        result = response.json()
        request = result['request']
        self.assertEqual(request['smiles'], data['smiles'])
        self.assertEqual(request['version'], 2)
        self.assertEqual(request['return_first'], True)
        self.assertEqual(request['chemical_property_logic'], 'none')
        self.assertEqual(request['chemical_popularity_logic'], 'none')
        self.assertEqual(request['template_set'], 'reaxys')
        self.assertEqual(request['template_prioritizer_version'], 0)
        self.assertEqual(request['priority'], 1)

        # Test that we got the celery task id
        self.assertIsInstance(result['task_id'], str)

        # Try retrieving task output
        result = self.get_result(result['task_id'])
        self.assertTrue(result['complete'])
        self.assertIsInstance(result['output'], list)
        self.assertIsInstance(result['output'][0], dict)
        self.assertIsInstance(result['output'][0]['children'], list)

        # Test task priority argument
        data = {
            'smiles': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
            'buyable_logic': 'or',
            'return_first': True,
            'version': 2,
            'priority': 2,
        }
        response = self.post('/tree-builder/', data=data)
        self.assertEqual(response.status_code, 200)

        result = response.json()
        request = result['request']
        self.assertEqual(request['priority'], 2)

    @classmethod
    def tearDownClass(cls):
        """This method is run once after all tests in this class."""
        cls.client.close()


if __name__ == '__main__':
    res = unittest.main(verbosity=3, exit=False)
