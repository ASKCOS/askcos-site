import numpy as np
from pymongo import MongoClient

import askcos.global_config as gc
from askcos.retrosynthetic.transformer import RetroTransformer

# cannot import from askcos_site.globals because RetroTransformerCelery is needed in globals
# that should be fixed, and this should be imported from globals
db_client = MongoClient(gc.MONGO['path'], gc.MONGO['id'], connect=gc.MONGO['connect'])
retro_templates = db_client[gc.RETRO_TEMPLATES['database']][gc.RETRO_TEMPLATES['collection']]

db_comparison_map = {'>': '$gt', '>=': '$gte', '<': '$lt', '<=': '$lte', '==': '$eq'}

class RetroTransformerCelery(RetroTransformer):
    """RetroTransformer with overridden methods for use with web services"""
    def __init__(
        self, load_all=False, template_set='reaxys', template_prioritizer=None,
        precursor_prioritizer='relevanceheuristic', fast_filter=None, 
        cluster='default', cluster_settings=None, template_db=None,
        scscorer=None
    ):
        super().__init__(
            load_all=load_all, template_set=template_set, 
            template_prioritizer=template_prioritizer, 
            precursor_prioritizer=precursor_prioritizer,
            fast_filter=fast_filter, cluster=cluster,
            cluster_settings=cluster_settings, scscorer=scscorer
        )

        self.template_db = template_db
        if self.template_db is None:
            self.template_db = retro_templates

    def get_one_template_by_idx(self, index, template_set=None):
        """Returns one template from given template set with given index.

        Args:
            index (int): index of template to return
            template_set (str): name of template set to return template from

        Returns:
            Template dictionary ready to be applied (i.e. - has 'rxn' object)

        """
        if template_set is None:
            template_set = self.template_set
            
        template = self.template_db.find_one(
            {
                'index': index,
                'template_set': template_set
            }
        )
        
        if not self.load_all:
            template = self.doc_to_template(template)

        if not template:
            raise ValueError('Could not find template from template set "{}" with index "{}"'.format(
                template_set, index
            ))

        return template

    def get_templates_by_indices(self, indices, template_set=None):
        """Returns templates from given template set with given indices

        Args:
            indices (np.array): indices of templates to return
            template_set (str, optional): Template set from which to 
                retrieve templates

        Returns:
            list: templates ready to be applied (with `rxn` attribute)

        """
        if template_set is None:
            template_set = self.template_set

        index_list = indices.tolist()

        cursor = self.template_db.find({
            'index': {'$in': index_list},
            'template_set': template_set
        })

        template_map = {x['index']: x for x in cursor}

        templates = [template_map[i] for i in index_list]

        if not self.load_all:
            # return generator of templates with rchiralReaction if rdchiralReaction initialization was successful
            templates = (x for x in (self.doc_to_template(temp) for temp in templates) if x.get('rxn'))

        return templates

    def filter_by_attributes(self, scores, indices, attribute_filter, template_set=None):
        """Filters template indices by attribute filter(s)

        Args:
            scores (np.array): scores predicted for prioritized templates
            indices (np.array): indices of prioritized templates
            attribute_filter (list[dict]): list of dictionaries defining 
                attribute filters. The format should be {'name': <str>, 
                'logic': <str>, 'value': <int/float>} where `logic` should be
                one of ['>', '>=', '<', '<=', '==].
            template_set (str, optional): NOT USED.

        Returns:
            np.array, np.array: scores, indices, of prioritized templates 
                following application of attribute filters
        """
        if template_set is None:
            template_set = self.template_set

        attribute_query = {
            'template_set': template_set
        }
        for query in attribute_filter:
            attribute_query['attributes.{}'.format(query['name'])] = {
                db_comparison_map[query['logic']]: query['value']
            }
        cursor = self.template_db.find(
            attribute_query,
            {'index': 1, '_id': 0}
        )
        filtered_indices = [doc['index'] for doc in cursor]
        bool_mask = np.isin(indices, filtered_indices)
        indices = indices[bool_mask]
        scores = scores[bool_mask]
        return scores, indices

    def lookup_id(self, template_id):
        """Looks up one template in mognodb by template id"""
        template = self.template_db.find_one({
            '_id': template_id
        })
        return template
