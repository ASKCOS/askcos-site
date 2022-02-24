"""
Tree builder subclass using celery for multiprocessing.
"""

import askcos_site.askcos_celery.treebuilder.tb_c_worker as tb_c_worker
from askcos.retrosynthetic.mcts.tree_builder import MCTS, WAITING
from askcos_site.askcos_celery.treebuilder.retro_transformer_celery import RetroTransformerCelery
from askcos_site.globals import retro_templates


class MCTSCelery(MCTS):
    """
    This is a subclass of MCTS which uses celery for multiprocessing.

    Note regarding model and data loading: This class uses pricing data,
    chemhistorian data, template prioritizer, and retro transformer. The retro
    transformer additionally needs the precursor prioritizer and fast filter.
    If instantiating this class with no arguments, only Pricer and ChemHistorian
    data will be loaded (using database) by default. The RetroTransformer
    is provided via ``tb_c_worker`` which is configured separately. The
    template prioritizer type should be provided to ``get_buyable_paths`` to
    indicate which model to use.

    Attributes:

    """

    def __init__(self, template_prioritizer=None, precursor_prioritizer=None, fast_filter=None, use_db=True, **kwargs):
        super().__init__(
            template_prioritizer=template_prioritizer,
            precursor_prioritizer=precursor_prioritizer,
            fast_filter=fast_filter,
            use_db=use_db,
            **kwargs
        )

        from celery.result import allow_join_result
        self.allow_join_result = allow_join_result
        self.template_prioritizer_version = None

    @staticmethod
    def load_retro_transformer(template_set='reaxys', precursor_prioritizer='relevanceheuristic'):
        """
        Loads retro transformer model.
        """
        retro_transformer = RetroTransformerCelery(
            template_set=template_set,
            template_prioritizer=None,
            precursor_prioritizer=precursor_prioritizer,
            fast_filter=None
        )
        retro_transformer.load(load_templates=False)
        return retro_transformer

    def reset_workers(self, soft_reset=False):
        # general parameters in celery format
        # TODO: anything goes here?
        self.pending_results = []

    def expand(self, _id, smiles, template_idx):  # TODO: make Celery workers
        """Adds pathway to be worked on with Celery.

        Args:
            _id (int): ID of pending pathway.
            smiles (str): SMILES string of molecule to be exanded.
            template_idx (int): ID of template to apply to molecule.
        """
        # Chiral transformation or heuristic prioritization requires
        # same database. _id is _id of active pathway
        self.pending_results.append(tb_c_worker.apply_one_template_by_idx.apply_async(
            args=(_id, smiles, template_idx),
            kwargs={'max_num_templates': self.template_count,
                    'max_cum_prob': self.max_cum_template_prob,
                    'fast_filter_threshold': self.filter_threshold,
                    'template_prioritizer_version': self.template_prioritizer_version,
                    'template_set': self.template_set},
            priority=2,  # Send high priority task
        ))
        self.status[(smiles, template_idx)] = WAITING
        self.active_pathways_pending[_id] += 1

    def prepare(self):
        """Starts parallelization with Celery."""
        try:
            res = tb_c_worker.apply_one_template_by_idx.apply_async(
                args=(1, 'CCOC(=O)[C@H]1C[C@@H](C(=O)N2[C@@H](c3ccccc3)CC[C@@H]2c2ccccc2)[C@@H](c2ccccc2)N1', 1),
                kwargs={'template_set': self.template_set,
                        'template_prioritizer_version': self.template_prioritizer_version},
                priority=2,
            )
            res.get(20)
        except Exception as e:
            res.revoke()
            raise IOError('Did not find any workers? Try again later ({})'.format(e))

    def get_ready_result(self):
        """Yields processed results from Celery.

        Yields:
            list of 5-tuples of (int, string, int, list, float): Results
                from workers after applying a template to a molecule.
        """
        # Update which processes are ready
        self.is_ready = [i for (i, res) in enumerate(self.pending_results) if res.ready()]
        for i in self.is_ready:
            yield self.pending_results[i].get(timeout=0.1)
            self.pending_results[i].forget()
        self.pending_results = [res for (i, res) in enumerate(self.pending_results) if i not in self.is_ready]

    def stop(self, soft_stop=False):
        """Stops work with Celery.

        Args:
            soft_stop (bool, optional): Unused. (default: {false})
        """
        self.running = False
        if self.pending_results != []:  # clear anything left over - might not be necessary
            for i in range(len(self.pending_results)):
                self.pending_results[i].revoke()

    def get_initial_prioritization(self):
        """
        Get template prioritizer predictions to initialize the tree search.
        """
        res = tb_c_worker.template_relevance.apply_async(
            args=(self.smiles, self.template_count, self.max_cum_template_prob),
            kwargs={'template_set': self.template_set,
                    'template_prioritizer_version': self.template_prioritizer_version},
            priority=2,  # Send high priority task
        )
        return res.get(10)

    def work(self, i):
        """
        Explicitly override work method of MCTS since Celery does not use it.
        """
        raise NotImplementedError('MCTSCelery does not support the work method. Did you mean to use MCTS?')

    def wait_until_ready(self):
        """
        No need to wait for Celery workers since they should be pre-initialized.
        """
        pass
