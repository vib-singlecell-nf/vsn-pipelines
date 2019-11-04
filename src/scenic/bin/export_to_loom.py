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

from typing import Mapping
from collections import OrderedDict


class SCopeLoom:

    def __init__(
        self,
        ex_mtx: pd.DataFrame, regulons: List[Regulon], out_fname: str,
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
        save_additional_regulon_meta_data: bool = False
    ):
        self.ex_mtx = ex_mtx
        self.regulons = regulons
        self.out_fname = out_fname
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
        # Utils
        self.id2name = OrderedDict()
        # Loom Data
        self.col_attrs = {}
        self.row_attrs = {}
        self.general_attrs = {}
        # Set
        if(regulons[0].name.find('_') == -1):
            print(
                "Regulon name does not seem to be compatible with SCOPE. It should include a space to allow selection of the TF.",
                "\nPlease run: \n regulons = [r.rename(r.name.replace('(+)','_('+str(len(r))+'g)')) for r in regulons]",
                "\nor:\n regulons = [r.rename(r.name.replace('(','_(')) for r in regulons]"
            )
        self.set_cell_annotations()
        if auc_mtx is None:
            self.auc_mtx = self.calculate_regulon_enrichment()
        if self.auc_thresholds is None:
            self.auc_thresholds = self.binarize_regulon_enrichment()
        if len(self.embeddings) == 0:
            self.embeddings = self.create_loom_default_embedding()
        loom_embeddings = self.create_loom_embeddings()
        self.embeddings_X = loom_embeddings["embeddings_X"]
        self.embeddings_Y = loom_embeddings["embeddings_Y"]
        self.regulon_gene_assignment = self.create_loom_regulon_gene_assignment()
        self.ngenes = self.calculate_nb_genes_per_cell()
        self.default_embedding = self.get_default_embedding()
        self.clusterings = self.create_loom_clusterings()
        self.regulon_thresholds = self.create_loom_md_regulon_thresholds()
        self.set_generic_columns_attrs()
        self.set_generic_row_attrs()
        self.set_generic_general_attrs()
        self.set_tree()

    def set_cell_annotations(self):
        if self.cell_annotations is None:
            self.cell_annotations = dict(zip(self.ex_mtx.index.astype(str), ['-'] * self.ex_mtx.shape[0]))

    def get_regulon_gene_data(self, regulon, key):
        if key == 'gene2occurrence':
            return regulon.gene2occurrence
        if key == 'gene2weight':
            return regulon.gene2weight
        raise Exception('Cannot retrieve {} from given regulon. Not implemented.'.format(key))

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

    def create_loom_default_embedding(self):
        # Create an embedding based on tSNE.
        # Name of columns should be "_X" and "_Y".
        embedding = {
            "tSNE (default)": pd.DataFrame(
                data=TSNE().fit_transform(self.auc_mtx),
                index=self.ex_mtx.index, columns=['_X', '_Y']
            )
        }  # (n_cells, 2)
        return embedding

    def create_loom_embeddings(self):
        self.embeddings_X = pd.DataFrame(index=self.ex_mtx.index)
        self.embeddings_Y = pd.DataFrame(index=self.ex_mtx.index)

        for idx, (name, df_embedding) in enumerate(self.embeddings.items()):
            if(len(df_embedding.columns) != 2):
                raise Exception('The embedding should have two columns.')

            embedding_id = idx - 1  # Default embedding must have id == -1 for SCope.
            self.id2name[embedding_id] = name

            embedding = df_embedding.copy()
            embedding.columns = ['_X', '_Y']
            return {
                "embeddings_X": pd.merge(self.embeddings_X, embedding['_X'].to_frame().rename(columns={'_X': str(embedding_id)}), left_index=True, right_index=True),
                "embeddings_Y": pd.merge(self.embeddings_Y, embedding['_Y'].to_frame().rename(columns={'_Y': str(embedding_id)}), left_index=True, right_index=True)
            }

    def calculate_nb_genes_per_cell(self):
        # Calculate the number of genes per cell.
        binary_mtx = self.ex_mtx.copy()
        binary_mtx[binary_mtx != 0] = 1.0
        return binary_mtx.sum(axis=1).astype(int)

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

    def get_default_embedding(self):
        default_embedding = next(iter(self.embeddings.values())).copy()
        default_embedding.columns = ['_X', '_Y']
        return default_embedding

    def name2idx(self):
        return dict(map(reversed, enumerate(sorted(set(self.cell_annotations.values())))))

    def create_loom_clusterings(self):
        # Encode cell type clusters.
        # The name of the column should match the identifier of the clustering.
        return pd.DataFrame(
            data=self.ex_mtx.index.values,
            index=self.ex_mtx.index,
            columns=['0']
        ).replace(self.cell_annotations).replace(self.name2idx())

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

    def set_generic_columns_attrs(self):
        self.column_attrs = {
            "CellID": self.ex_mtx.index.values.astype('str'),
            "nGene": self.ngenes.values,
            "Embedding": SCopeLoom.create_structure_array(self.default_embedding),
            "RegulonsAUC": SCopeLoom.create_structure_array(self.auc_mtx),
            "Clusterings": SCopeLoom.create_structure_array(self.clusterings),
            "ClusterID": self.clusterings.values,
            'Embeddings_X': SCopeLoom.create_structure_array(self.embeddings_X),
            'Embeddings_Y': SCopeLoom.create_structure_array(self.embeddings_Y),
        }

    def set_generic_row_attrs(self):
        self.row_attrs = {
            "Gene": self.ex_mtx.columns.values.astype('str'),
            "Regulons": SCopeLoom.create_structure_array(self.regulon_gene_assignment)
        }

    def set_generic_general_attrs(self):
        self.general_attrs = {
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
        self.general_attrs.update(("SCopeTreeL{}".format(idx + 1), category) for idx, category in enumerate(list(islice(chain(self.tree_structure, repeat("")), 3))))

    def add_meta_data(self, _dict):
        md = self.general_attrs["MetaData"]
        meta_data = json.loads(md)
        meta_data.update(_dict)
        self.general_attrs["MetaData"] = json.dumps(meta_data)

    def export(self):
        # Compress MetaData global attribute
        if self.compress:
            self.general_attrs["MetaData"] = SCopeLoom.compress_encode(value=self.general_attrs["MetaData"])

        # Create loom file for use with the SCope tool.
        # The loom file format opted for rows as genes to facilitate growth along the column axis (i.e add more cells)
        # PySCENIC chose a different orientation because of limitation set by the feather format: selectively reading
        # information from disk can only be achieved via column selection. For the ranking databases this is of utmost
        # importance.
        lp.create(
            filename=self.out_fname,
            layers=self.ex_mtx.T.values,
            row_attrs=self.row_attrs,
            col_attrs=self.column_attrs,
            file_attrs=self.general_attrs
        )

    # Multi-runs SCENIC additional data

    def encode_regulon_gene_data(self, key):
        genes = np.array(self.ex_mtx.columns)
        n_genes = len(genes)
        n_regulons = len(self.regulons)
        data = np.zeros(shape=(n_genes, n_regulons), dtype=float)
        for idx, regulon in enumerate(self.regulons):
            regulon_data = self.get_regulon_gene_data(regulon, key)
            regulon_genes = list(dict(map(lambda x: x, regulon_data.items())).keys())
            data[np.isin(genes, regulon_genes), idx] = list(dict(map(lambda x: x, regulon_data.items())).values())
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
            self.row_attrs['RegulonGeneWeights'] = SCopeLoom.create_structure_array(regulon_gene_weights)

    def add_row_attr_regulon_gene_occurrences(self):
        regulon_gene_occurrences = self.create_loom_regulon_gene_data(key='gene2occurrence')
        if 'RegulonGeneOccurrences' not in self.row_attrs.keys():
            print("Added 'RegulonGeneOccurrences' to the row attributes.")
            self.row_attrs['RegulonGeneOccurrences'] = SCopeLoom.create_structure_array(regulon_gene_occurrences)

    # Utility functions

    @staticmethod
    def create_structure_array(df):
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
