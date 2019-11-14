import os
import numpy as np
import pandas as pd
import loompy as lp
from sklearn.manifold.t_sne import TSNE
from pyscenic.aucell import aucell
from pyscenic.genesig import Regulon, GeneSignature
from typing import List, Mapping, Union, Sequence, Optional
from operator import attrgetter
from multiprocessing import cpu_count
from pyscenic.binarization import binarize
from itertools import chain, repeat, islice
import networkx as nx
import zlib
import base64
import json
import re

from typing import Mapping
from collections import OrderedDict


class Embedding:

    def __init__(self, embedding: pd.DataFrame, embedding_name: str, is_default: bool = False):
        self.embedding = embedding
        self.embedding_name = embedding_name
        self._is_default = is_default

    def get_embedding_name(self):
        return self.embedding_name

    def get_embedding(self):
        return self.embedding

    def is_default(self):
        return self._is_default


class SCopeLoom:

    def __init__(
        self,
        out_fname: str = None,
        ex_mtx: pd.DataFrame = None,
        regulons: List[Regulon] = None,
        cell_annotations: Optional[Mapping[str, str]] = None,
        tree_structure: Sequence[str] = (),
        title: Optional[str] = None,
        nomenclature: str = "Unknown",
        num_workers: int = cpu_count(),
        embeddings: Mapping[str, pd.DataFrame] = {},
        auc_mtx=None,
        auc_regulon_weights_key='gene2weight',
        auc_thresholds=None,
        compress: bool = False,
        save_additional_regulon_meta_data: bool = False,  # Should be set to true only for multi-runs SCENIC
        set_base_loom: bool = False,
        tag: str = None  # Used when merging track and motif-based SCENIC runs
    ):
        self.out_fname = out_fname
        if ex_mtx is not None:
            self.ex_mtx = ex_mtx
        if regulons is not None:
            self.regulons = regulons
        self.cell_annotations = cell_annotations
        self.tree_structure = tree_structure
        self.title = title
        self.nomenclature = nomenclature
        self.num_workers = num_workers
        self.embeddings = embeddings
        self.auc_mtx = auc_mtx
        self.auc_regulon_weights_key = auc_regulon_weights_key
        self.auc_thresholds = auc_thresholds
        self.compress = compress
        self.save_additional_regulon_meta_data = save_additional_regulon_meta_data
        self.tag = tag
        # Utils
        self.id2name = OrderedDict()
        # Loom Data
        self.col_attrs = {}
        self.row_attrs = {}
        self.global_attrs = {}
        # Set Base Loom
        if set_base_loom:
            self.set_base_loom()

    def set_base_loom(self):
        if(self.regulons[0].name.find('_') == -1):
            print(
                "Regulon name does not seem to be compatible with SCOPE. It should include a space to allow selection of the TF.",
                "\nPlease run: \n regulons = [r.rename(r.name.replace('(+)','_('+str(len(r))+'g)')) for r in regulons]",
                "\nor:\n regulons = [r.rename(r.name.replace('(','_(')) for r in regulons]"
            )
        self.set_cell_annotations()
        if self.auc_mtx is None:
            self.auc_mtx = self.calculate_regulon_enrichment()
        if self.auc_thresholds is None:
            self.auc_thresholds = self.binarize_regulon_enrichment()
        if len(self.embeddings.keys()) == 0:
            self.create_loom_default_embedding()
        self.regulon_gene_assignment = self.create_loom_regulon_gene_assignment()
        self.ngenes = self.calculate_nb_genes_per_cell()
        self.clusterings = self.create_loom_clusterings()
        self.regulon_thresholds = self.create_loom_md_regulon_thresholds()
        self.set_generic_col_attrs()
        self.set_generic_row_attrs()
        self.set_generic_global_attrs()
        self.set_tree()

    #######
    # I/O #
    #######

    @staticmethod
    def read_loom(filename: str, tag: str = None):
        with lp.connect(filename, mode='r', validate=False) as loom:
            scope_loom = SCopeLoom(tag=tag)
            # Materialize i.e.: load the loom content otherwise storage is empty
            lp.to_html(loom)
            # Set the main matrix
            scope_loom.ex_mtx = pd.DataFrame(loom[:, :], index=loom.ra.Gene, columns=loom.ca.CellID).T
            # Set the column, row and global attribute using the underlying Dict of the AttributeManager
            scope_loom.col_attrs = loom.ca.__dict__["storage"]
            scope_loom.row_attrs = loom.ra.__dict__["storage"]
            scope_loom.global_attrs = loom.attrs.__dict__["storage"]
            # Decompress and decode the MetaData global attribute
            try:
                scope_loom.global_attrs["MetaData"] = SCopeLoom.decompress_decode(value=scope_loom.global_attrs["MetaData"])
            except Exception:
                # MetaData is uncompressed
                scope_loom.global_attrs["MetaData"] = json.loads(scope_loom.global_attrs["MetaData"])
        return scope_loom

    ##############
    # Embeddings #
    ##############

    def add_embedding(self, embedding: np.ndarray, embedding_name, is_default: bool = False):
        df_embedding = pd.DataFrame(embedding, columns=['_X', '_Y'], index=self.ex_mtx.index)
        _embedding = Embedding(embedding=df_embedding, embedding_name=embedding_name, is_default=is_default)
        if is_default:
            self.default_embedding = _embedding
        self.embeddings[embedding_name] = _embedding

    @staticmethod
    def get_embedding_id(embedding: Embedding, _list):
        if embedding.is_default():
            return '-1'
        elif len(_list) == '0':
            return '0'
        else:
            return str(len(_list) - 1)

    def create_loom_md_embeddings(self):
        md_embeddings = []
        for _, embedding in self.embeddings.items():
            if embedding.is_default():
                md_embeddings = md_embeddings + [{'id': '-1', 'name': embedding.get_embedding_name()}]
            else:
                md_embeddings = md_embeddings + [
                    {
                        'id': SCopeLoom.get_embedding_id(embedding=embedding, _list=md_embeddings),
                        'name': embedding.get_embedding_name()
                    }
                ]
        return {"embeddings": md_embeddings}

    def create_loom_ca_embeddings(self):
        """ Returns a Dictionary (Dict) with preformated data to be stored in Loom file format.

        Parameters:
            None

        Returns:
            dict: A Dictionary (Dict) with preformated data to be stored in Loom file format.

        """
        default_embedding = None
        embeddings_X = pd.DataFrame(index=self.ex_mtx.index)
        embeddings_Y = pd.DataFrame(index=self.ex_mtx.index)

        for _, embedding in self.embeddings.items():
            if(embedding.get_embedding().shape[1] != 2):
                raise Exception('The embedding should have two columns.')
            # Set the default embedding
            if embedding.is_default():
                default_embedding = embedding.get_embedding()
            # Update the Embeddings_[X|Y]
            embedding_id = str(SCopeLoom.get_embedding_id(embedding, embeddings_X.columns))
            embedding = embedding.get_embedding().copy()
            embedding.columns = ['_X', '_Y']
            embeddings_X = pd.merge(
                embeddings_X,
                embedding['_X'].to_frame().rename(
                    columns={'_X': embedding_id}
                ).astype('float32'), left_index=True, right_index=True)
            embeddings_Y = pd.merge(
                embeddings_Y,
                embedding['_Y'].to_frame().rename(
                    columns={'_Y': embedding_id}
                ).astype('float32'), left_index=True, right_index=True)
        return {
            "Embedding": SCopeLoom.df_to_named_matrix(df=default_embedding),
            "Embeddings_X": SCopeLoom.df_to_named_matrix(df=embeddings_X),
            "Embeddings_Y": SCopeLoom.df_to_named_matrix(df=embeddings_Y)
        }

    def create_loom_default_embedding(self):
        # Create an embedding based on tSNE.
        # Name of columns should be "_X" and "_Y".
        self.add_embedding(embedding=TSNE().fit_transform(self.auc_mtx), embedding_name="tSNE (default)", is_default=True)

    ###############
    # Clusterings #
    ###############

    def create_loom_clusterings(self):
        # Encode cell type clusters.
        # The name of the column should match the identifier of the clustering.
        return pd.DataFrame(
            data=self.ex_mtx.index.values,
            index=self.ex_mtx.index,
            columns=['0']
        ).replace(self.cell_annotations).replace(self.name2idx())

    ###############
    # Annotations #
    ###############

    def set_cell_annotations(self):
        if self.cell_annotations is None:
            self.cell_annotations = dict(zip(self.ex_mtx.index.astype(str), ['-'] * self.ex_mtx.shape[0]))

    ###########
    # Metrics #
    ###########

    def calculate_nb_genes_per_cell(self):
        # Calculate the number of genes per cell.
        binary_mtx = self.ex_mtx.copy()
        binary_mtx[binary_mtx != 0] = 1.0
        return binary_mtx.sum(axis=1).astype(int)

    def add_metrics(self, metrics: List[str]):
        md_metrics = []
        for metric in metrics:
            md_metrics.append({"name": metric})
        self.global_attrs["MetaData"].update({'metrics': md_metrics})

    ##########
    # SCENIC #
    ##########

    def fix_loom_md_regulon_data_for_scope(self):
        # Fix regulon objects to display properly in SCope:
        # Rename regulons in the thresholds object, motif
        md_regulon_thresholds = self.global_attrs['MetaData']['regulonThresholds']
        for _, x in enumerate(md_regulon_thresholds):
            tmp = re.sub(r"(_?)\(", '_(', x.get('regulon'))  # + '-motif'
            x.update({'regulon': tmp})
        return {
            'regulonThresholds': md_regulon_thresholds
        }

    def fix_loom_ca_regulon_data_for_scope(self):
        regulons_auc_col_attrs = list(filter(lambda col_attrs_key: 'RegulonsAUC' in col_attrs_key, self.col_attrs.keys()))

        def fix(col_attrs_key):
            regulons = pd.DataFrame(self.col_attrs[col_attrs_key], index=self.col_attrs['CellID'])
            # Add underscore for SCope compatibility:
            regulons.columns = regulons.columns.str.replace('_?\\(', '_(')
            return {
                col_attrs_key: SCopeLoom.df_to_named_matrix(regulons)
            }
        # Update the keys
        regulons_auc_col_attrs_update = map(fix, regulons_auc_col_attrs)
        # Convert list of dict to dict
        return {next(iter(x)): x.get(next(iter(x))) for x in regulons_auc_col_attrs_update}

    def fix_loom_ra_regulon_data_for_scope(self):
        regulons_row_attrs = list(filter(lambda row_attrs_key: 'Regulons' in row_attrs_key, self.row_attrs.keys()))

        def fix(row_attrs_key):
            regulons = pd.DataFrame(self.row_attrs[row_attrs_key], index=self.row_attrs['Gene'])
            # Add underscore for SCope compatibility:
            regulons.columns = regulons.columns.str.replace('_?\\(', '_(')
            return {
                row_attrs_key: SCopeLoom.df_to_named_matrix(regulons)
            }
        # Update the keys
        regulons_row_attrs_update = map(fix, regulons_row_attrs)
        # Convert list of dict to dict
        return {next(iter(x)): x.get(next(iter(x))) for x in regulons_row_attrs_update}

    def calculate_regulon_enrichment(self):
        # Calculate regulon enrichment per cell using AUCell.
        # Create regulons with weight based on given key
        print("Using {} to weight the genes when running AUCell.".format(self.auc_regulon_weights_key))
        regulon_signatures = list(map(lambda x: GeneSignature(name=x.name, gene2weight=self.get_regulon_gene_data(x, self.auc_regulon_weights_key)), self.regulons))
        auc_mtx = aucell(self.ex_mtx, regulon_signatures, num_workers=self.num_workers)  # (n_cells x n_regulons)
        auc_mtx = auc_mtx.loc[self.ex_mtx.index]
        return auc_mtx

    def binarize_regulon_enrichment(self):
        _, auc_thresholds = binarize(self.auc_mtx)
        return auc_thresholds

    @staticmethod
    def is_SCENIC_multi_runs_mode(scopeloom):
        return 'RegulonGeneOccurrences' in scopeloom.row_attrs.keys() and 'RegulonGeneWeights' in scopeloom.row_attrs.keys()

    def merge_regulon_data(self, scope_loom):
        # Check if SCENIC has been run in multi-runs mode
        is_multi_runs_mode = SCopeLoom.is_SCENIC_multi_runs_mode(self) and SCopeLoom.is_SCENIC_multi_runs_mode(scope_loom)

        # RegulonsAUC
        # Relabel columns with suffix indicating the regulon source
        auc_mtx = pd.DataFrame(data=self.col_attrs['RegulonsAUC'], index=self.col_attrs['CellID'])
        auc_mtx.columns = auc_mtx.columns + '-' + self.tag

        scope_loom_auc_mtx = pd.DataFrame(data=scope_loom.col_attrs['RegulonsAUC'], index=scope_loom.col_attrs['CellID'])
        scope_loom_auc_mtx.columns = scope_loom_auc_mtx.columns + '-' + scope_loom.tag

        # Regulons (regulon assignment matrices)
        regulons = pd.DataFrame(self.row_attrs['Regulons'], index=self.row_attrs['Gene'])
        regulons.columns = regulons.columns + '-' + self.tag

        scope_loom_regulons = pd.DataFrame(scope_loom.row_attrs['Regulons'], index=scope_loom.row_attrs['Gene'])
        scope_loom_regulons.columns = scope_loom_regulons.columns + '-' + scope_loom.tag

        # If multi-runs SCENIC
        # Rename Regulons Gene Occurrences
        if is_multi_runs_mode:
            regulon_gene_occurrences = pd.DataFrame(self.row_attrs['RegulonGeneOccurrences'], index=self.row_attrs['Gene'])
            regulon_gene_occurrences.columns = regulon_gene_occurrences.columns + '-' + self.tag

            scope_loom_regulon_gene_occurrences = pd.DataFrame(scope_loom.row_attrs['RegulonGeneOccurrences'], index=scope_loom.row_attrs['Gene'])
            scope_loom_regulon_gene_occurrences.columns = scope_loom_regulon_gene_occurrences.columns + '-' + scope_loom.tag

        # If multi-runs SCENIC
        # Rename Regulons Gene Weights
        if is_multi_runs_mode:
            regulon_gene_weights = pd.DataFrame(self.row_attrs['RegulonGeneWeights'], index=self.row_attrs['Gene'])
            regulon_gene_weights.columns = regulon_gene_weights.columns + '-' + self.tag

            scope_loom_regulon_gene_weights = pd.DataFrame(scope_loom.row_attrs['RegulonGeneWeights'], index=scope_loom.row_attrs['Gene'])
            scope_loom_regulon_gene_weights.columns = scope_loom_regulon_gene_weights.columns + '-' + scope_loom.tag

        # Combine meta data regulons
        # Rename regulons in the thresholds object, motif
        rt = self.global_attrs["MetaData"]["regulonThresholds"]
        for _, x in enumerate(rt):
            tmp = x.get('regulon') + '-' + self.tag
            x.update({'regulon': tmp})

        # Rename regulons in the thresholds object, track
        scope_loom_rt = scope_loom.global_attrs["MetaData"]["regulonThresholds"]
        for _, x in enumerate(scope_loom_rt):
            tmp = x.get('regulon') + '-' + scope_loom.tag
            x.update({'regulon': tmp})
            # blank out the "motifData" field for track-based regulons:
            x.update({'mofitData': 'NA.png'})

        # merge regulon threshold dictionaries:
        rt_merged = rt + scope_loom_rt

        # Remove because we will save them separately
        del self.row_attrs["Regulons"]
        del self.col_attrs["RegulonsAUC"]

        if is_multi_runs_mode:
            del self.row_attrs["RegulonGeneOccurrences"]
            del self.row_attrs["RegulonGeneWeights"]

        # Update the attributes
        self.row_attrs.update({f'{self.tag.capitalize()}Regulons': SCopeLoom.df_to_named_matrix(regulons)})
        self.col_attrs.update({f'{self.tag.capitalize()}RegulonsAUC': SCopeLoom.df_to_named_matrix(auc_mtx)})

        if is_multi_runs_mode:
            self.row_attrs.update({f'{self.tag.capitalize()}RegulonGeneOccurrences': SCopeLoom.df_to_named_matrix(regulon_gene_occurrences)})
            self.row_attrs.update({f'{self.tag.capitalize()}RegulonGeneWeights': SCopeLoom.df_to_named_matrix(regulon_gene_weights)})

        self.row_attrs.update({f'{scope_loom.tag.capitalize()}Regulons': SCopeLoom.df_to_named_matrix(scope_loom_regulons)})
        self.col_attrs.update({f'{scope_loom.tag.capitalize()}RegulonsAUC': SCopeLoom.df_to_named_matrix(scope_loom_auc_mtx)})

        if is_multi_runs_mode:
            self.row_attrs.update({f'{scope_loom.tag.capitalize()}RegulonGeneOccurrences': SCopeLoom.df_to_named_matrix(scope_loom_regulon_gene_occurrences)})
            self.row_attrs.update({f'{scope_loom.tag.capitalize()}RegulonGeneWeights': SCopeLoom.df_to_named_matrix(scope_loom_regulon_gene_weights)})

        self.global_attrs["MetaData"].update({'regulonThresholds': rt_merged})

    def get_regulon_gene_data(self, regulon, key):
        if key == 'gene2occurrence':
            return regulon.gene2occurrence
        if key == 'gene2weight':
            return regulon.gene2weight
        raise Exception('Cannot retrieve {} from given regulon. Not implemented.'.format(key))

    def get_regulon_meta_data(self, name, threshold):

        def fetch_logo(context):
            for elem in context:
                if elem.endswith('.png'):
                    return elem
            return ""

        name2logo = {reg.name: fetch_logo(reg.context) for reg in self.regulons}

        regulon = list(filter(lambda x: x.name == name, self.regulons))[0]
        regulon_meta_data = {
            "regulon": name,
            "defaultThresholdValue": (threshold if isinstance(threshold, float) else threshold[0]),
            "defaultThresholdName": "guassian_mixture_split",
            "allThresholds": {
                "guassian_mixture_split": (threshold if isinstance(threshold, float) else threshold[0])
            },
            "motifData": name2logo.get(name, "")
        }
        if self.save_additional_regulon_meta_data:
            regulon_meta_data.update({
                "tf": regulon.transcription_factor,
                "score": regulon.score if hasattr(regulon, 'score') else 0.0,
                "orthologousIdentity": regulon.orthologous_identity if hasattr(regulon, 'orthologous_identity') else 0.0,
                "annotation": regulon.annotation if hasattr(regulon, 'annotation') else '',
                "similarityQValue": regulon.score if hasattr(regulon, 'similarity_qvalue') else 0.0
            })
        return regulon_meta_data

    def create_loom_md_regulon_thresholds(self):
        return [self.get_regulon_meta_data(name, threshold) for name, threshold in self.auc_thresholds.iteritems()]

    def create_loom_regulon_gene_assignment(self):
        # Encode genes in regulons as "binary" membership matrix.
        genes = np.array(self.ex_mtx.columns)
        n_genes = len(genes)
        n_regulons = len(self.regulons)
        data = np.zeros(shape=(n_genes, n_regulons), dtype=int)
        for idx, regulon in enumerate(self.regulons):
            data[:, idx] = np.isin(genes, regulon.genes).astype(int)
        regulon_gene_assignment = pd.DataFrame(
            data=data,
            index=self.ex_mtx.columns,
            columns=list(map(attrgetter('name'), self.regulons))
        )
        return regulon_gene_assignment

    # Multi-runs SCENIC additional data

    def encode_regulon_gene_data(self, key):
        genes = np.array(self.ex_mtx.columns)
        n_genes = len(genes)
        n_regulons = len(self.regulons)
        data = np.zeros(shape=(n_genes, n_regulons), dtype=float)
        for idx, regulon in enumerate(self.regulons):
            regulon_data = pd.DataFrame.from_dict(self.get_regulon_gene_data(regulon, key), orient='index')
            data[np.isin(genes, regulon_data.index), idx] = regulon_data[0][genes[np.isin(genes, regulon_data.index)]]
        return data

    def create_loom_regulon_gene_data(self, key):
        return pd.DataFrame(
            data=self.encode_regulon_gene_data(key=key),
            index=self.ex_mtx.columns,
            columns=list(map(attrgetter('name'), self.regulons))
        )

    def add_row_attr_regulon_gene_weights(self):
        regulon_gene_weights = self.create_loom_regulon_gene_data(key='gene2weight')
        if 'RegulonGeneWeights' not in self.row_attrs.keys():
            print("Added 'RegulonGeneWeights' to the row attributes.")
            self.row_attrs['RegulonGeneWeights'] = SCopeLoom.df_to_named_matrix(regulon_gene_weights)

    def add_row_attr_regulon_gene_occurrences(self):
        regulon_gene_occurrences = self.create_loom_regulon_gene_data(key='gene2occurrence')
        if 'RegulonGeneOccurrences' not in self.row_attrs.keys():
            print("Added 'RegulonGeneOccurrences' to the row attributes.")
            self.row_attrs['RegulonGeneOccurrences'] = SCopeLoom.df_to_named_matrix(regulon_gene_occurrences)

    ####

    def name2idx(self):
        return dict(map(reversed, enumerate(sorted(set(self.cell_annotations.values())))))

    def set_generic_col_attrs(self):
        self.col_attrs = {
            "CellID": self.ex_mtx.index.values.astype('str'),
            "nGene": self.ngenes.values,
            "Embedding": SCopeLoom.df_to_named_matrix(self.default_embedding.get_embedding()),
            "RegulonsAUC": SCopeLoom.df_to_named_matrix(self.auc_mtx),
            "Clusterings": SCopeLoom.df_to_named_matrix(self.clusterings),
            "ClusterID": self.clusterings.values
        }

    def set_generic_row_attrs(self):
        self.row_attrs = {
            "Gene": self.ex_mtx.columns.values.astype('str'),
            "Regulons": SCopeLoom.df_to_named_matrix(self.regulon_gene_assignment)
        }

    def set_generic_global_attrs(self):
        self.global_attrs = {
            "title": os.path.splitext(os.path.basename(self.out_fname))[0] if self.title is None else self.title,
            "MetaData": json.dumps({
                "embeddings": [{'id': identifier, 'name': name} for identifier, name in self.id2name.items()],
                "annotations": [{
                    "name": "",
                    "values": []
                }],
                "clusterings": [{
                    "id": 0,
                    "group": "celltype",
                    "name": "Cell Type",
                    "clusters": [{"id": idx, "description": name} for name, idx in self.name2idx().items()]
                }],
                "regulonThresholds": self.regulon_thresholds
            }),
            "Genome": self.nomenclature
        }

    def set_tree(self):
        assert len(self.tree_structure) <= 3, ""
        self.global_attrs.update(("SCopeTreeL{}".format(idx + 1), category) for idx, category in enumerate(list(islice(chain(self.tree_structure, repeat("")), 3))))

    def add_meta_data(self, _dict):
        md = self.global_attrs["MetaData"]
        meta_data = json.loads(md)
        meta_data.update(_dict)
        self.global_attrs["MetaData"] = meta_data

    def export(self, out_fname: str, save_embeddings: bool = True, compress_meta_data: bool = False):
        if out_fname is None:
            raise ValueError("The given out_fname cannot be None.")

        # Embeddings
        if save_embeddings:
            self.col_attrs.update(self.create_loom_ca_embeddings())
            self.global_attrs["MetaData"].update(self.create_loom_md_embeddings())

        # SCENIC
        if any('Regulons' in s for s in self.row_attrs.keys()):
            self.row_attrs.update(self.fix_loom_ra_regulon_data_for_scope())
        if any('RegulonsAUC' in s for s in self.col_attrs.keys()):
            self.col_attrs.update(self.fix_loom_ca_regulon_data_for_scope())
        if 'regulonThresholds' in self.global_attrs["MetaData"].keys():
            self.global_attrs["MetaData"].update(self.fix_loom_md_regulon_data_for_scope())

        # Compress MetaData global attribute
        # Should be compressed if Loompy version 2
        self.global_attrs["MetaData"] = json.dumps(self.global_attrs["MetaData"])
        if compress_meta_data:
            self.global_attrs["MetaData"] = SCopeLoom.compress_encode(value=self.global_attrs["MetaData"])

        # Create loom file for use with the SCope tool.
        # The loom file format opted for rows as genes to facilitate growth along the column axis (i.e add more cells)
        # PySCENIC chose a different orientation because of limitation set by the feather format: selectively reading
        # information from disk can only be achieved via column selection. For the ranking databases this is of utmost
        # importance.
        lp.create(
            filename=out_fname,
            layers=self.ex_mtx.T.values,
            row_attrs=self.row_attrs,
            col_attrs=self.col_attrs,
            file_attrs=self.global_attrs
        )

    # Utility functions

    @staticmethod
    def df_to_named_matrix(df: pd.DataFrame):
        # Create meta-data structure.
        # Create a numpy structured array
        return np.array([tuple(row) for row in df.values],
                        dtype=np.dtype(list(zip(df.columns, df.dtypes))))

    @staticmethod
    def compress_encode(value):
        '''
        Compress using ZLIB algorithm and encode the given value in base64.
        Taken from: https://github.com/aertslab/SCopeLoomPy/blob/5438da52c4bcf48f483a1cf378b1eaa788adefcb/src/scopeloompy/utils/__init__.py#L7
        '''
        return base64.b64encode(zlib.compress(value.encode('ascii'))).decode('ascii')

    @staticmethod
    def decompress_decode(value):
        try:
            value = value.decode('ascii')
            return json.loads(zlib.decompress(base64.b64decode(value)))
        except AttributeError:
            return json.loads(zlib.decompress(base64.b64decode(value.encode('ascii'))).decode('ascii'))
