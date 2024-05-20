from _types import *
import matplotlib.pyplot as plotter

import os
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Structure import Structure
from typing import List, Tuple,Iterable
import numpy as np
from scipy.spatial import distance_matrix
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA
import matplotlib.pyplot as plotter
import warnings

from Bio.Blast import NCBIWWW, NCBIXML
import requests

import plotly
import plotly.graph_objects as go
import plotly.figure_factory as ff

import warnings
class AAComposition():
    ZERO_AMINO_ACID_PROFILE = {
        "Alanine": 0,
        "Arginine": 0,
        "Asparagine": 0,
        "Aspartate": 0,
        "Cysteine": 0,
        "Glutamate": 0,
        "Glutamine": 0,
        "Glycine": 0,
        "Histidine": 0,
        "Isoleucine": 0,
        "Leucine": 0,
        "Lysine": 0,
        "Methionine": 0,
        "Phenylalanine": 0,
        "Proline": 0,
        "Serine": 0,
        "Threonine": 0,
        "Tryptophan": 0,
        "Tyrosine": 0,
        "Valine": 0
        }
    
    def __init__(self, 
            sequence: str = None
            ) -> None:
        self._sequence = sequence
        self._amino_acid_histogram_profile = None
        self._amino_acid_composition_profile = None
        self._total_amino_acid_count = 0
    
    @property
    def sequence(self) -> str:
        return self._sequence
    
    @sequence.setter
    def sequence(self, value: str = None) -> None:
        self._sequence = value
        self._amino_acid_histogram_profile = None
        self._amino_acid_composition_profile = None
        self._total_amino_acid_count = 0
    
    @property
    def amino_acid_histogram_profile(self) -> dict:
        if self._amino_acid_histogram_profile is None:
            self._compute_amino_acid_histogram_profile()
        return self._amino_acid_histogram_profile
    
    @property
    def amino_acid_composition_profile(self) -> dict:
        if self._amino_acid_histogram_profile is None:
            self._compute_amino_acid_histogram_profile()
        if self._amino_acid_composition_profile is None:
            self._compute_amino_acid_composition_profile()
        return self._amino_acid_composition_profile
    
    def _compute_amino_acid_histogram_profile(self) -> None:
        self._amino_acid_histogram_profile = self.ZERO_AMINO_ACID_PROFILE
        for one_letter_code in self.sequence:
            amino_acid = AminoAcids.amino_acids_1_letter_code[one_letter_code]
            self._amino_acid_histogram_profile[amino_acid.name] += 1
            self._total_amino_acid_count += 1
    
    def _compute_amino_acid_composition_profile(self) -> None:
        self._amino_acid_composition_profile = {}
        for amino_acid in self._amino_acid_histogram_profile.keys():
            amino_acid_count = self._amino_acid_histogram_profile[amino_acid]
            if not amino_acid_count == 0:
                raw_value = (amino_acid_count/self._total_amino_acid_count)
                composition_value = round(raw_value, 3)
            else:
                composition_value = 0
            percentage_value = 100 * composition_value
            self._amino_acid_composition_profile[amino_acid] = percentage_value
    
    def display_amino_acid_composition_profile(self, 
            order_type: str = None, 
            anotation_state: bool = True
            ) -> None:
        displayable_profile = self._get_displayable_profile(
            profile= self.amino_acid_composition_profile,
            order_type= order_type
            )
        self._plot_barh_chart(
            labels= list(displayable_profile.keys()),
            values= list(displayable_profile.values()),
            title= "Amino Acid Composition Profile",
            x_label= "Rate (%)",
            y_label= "Amino Acids",
            annotation_state= anotation_state,
            annotation_label_padding_size= (0.25, -0.25)
            )
    
    def _get_sorted_profile(self, profile, order_type):
        if order_type == 'ascending':
            return dict(sorted(profile.items(), key=lambda item: item[1]))
        elif order_type == 'descending':
            return dict(sorted(profile.items(), key=lambda item: item[1], reverse=True))
        else:
            return profile

    def _get_displayable_profile(self, profile=None, order_type=None):
        if profile is None:
            profile = self.amino_acid_histogram_profile
        return self._get_sorted_profile(profile, order_type)

    def display_amino_acid_histogram_profile(self, order_type=None, annotation_state=True):
        displayable_profile = self._get_displayable_profile(order_type=order_type)
        labels = list(displayable_profile.keys())
        values = list(displayable_profile.values())

        fig = go.Figure(go.Bar(
            x=values,
            y=labels,
            orientation='h'
        ))

        if annotation_state:
            for i, value in enumerate(values):
                fig.add_annotation(
                    x=value,
                    y=labels[i],
                    text=str(value),
                    showarrow=False,
                    xshift=0.75 if order_type != 'descending' else -0.25
                )

        fig.update_layout(
            title="Amino Acid Histogram Profile",
            xaxis_title="Count",
            yaxis_title="Amino Acids"
        )

        fig.show()

    
    def _plot_barh_chart(
        self,
        labels: list, 
        values: list,
        title: str, 
        x_label: str, 
        y_label: str, 
        annotation_state: bool = True,
        annotation_label_padding_size: tuple = (0, 0)
        ) -> None:
    # Create the bar chart
        fig = go.Figure(go.Bar(
            x=values,
            y=labels,
            orientation='h'
        ))

        # Add annotations if required
        if annotation_state:
            annotations = []
            for index, value in enumerate(values):
                annotations.append(dict(
                    x=value + annotation_label_padding_size[0],
                    y=index + annotation_label_padding_size[1],
                    xref='x',
                    yref='y',
                    text=str(round(value, 1)),
                    showarrow=False,
                    font=dict(
                        family='Arial, sans-serif',
                        size=12,
                        color='black'
                    ),
                    align='center',
                    ax=20,
                    ay=-30
                ))
            fig.update_layout(annotations=annotations)

        # Update layout for labels and title
        fig.update_layout(
            title=title,
            xaxis_title=x_label,
            yaxis_title=y_label,
            yaxis=dict(
                autorange='reversed'  # This is to match the 'barh' orientation
            )
        )

        # Show the plot
        fig.show()

    
    def _sort_by_value(
            self,
            profile: dict = None, 
            descending_order: bool = True
            ) -> dict:
        sorted_amino_acids_by_value = sorted(
            profile, 
            key= profile.get, 
            reverse= descending_order
            )
        sorted_profile = {}
        for amino_acid in sorted_amino_acids_by_value:
            sorted_profile[amino_acid] = profile[amino_acid]
        return sorted_profile



class ENMA():
    
    def __init__(self, 
            reference_file_path: str, 
            ensembles_files_paths: Tuple[str, ...]
            ) -> None:
        warnings.filterwarnings("ignore")
        self._reference_file_path = reference_file_path
        self._ensembles_files_paths = ensembles_files_paths
        self._reference = None
        self._ensemble = None
        self._analyse_results = None
        self._ensemble_analyse_records = None
    
    @property
    def reference(self) -> Reference:
        return self._reference
    
    @property
    def ensemble(self) -> Ensemble:
        return self._ensemble
    
    @property
    def analyse_results(self) -> AnalyseResults:
        return self._analyse_results
    
    @property
    def ensemble_analyse_records(self) -> EnsembleAnalyseRecords:
        return self._ensemble_analyse_records
    
    def perform_ensemble_normal_mode_analyse(self) -> None:
        ENMA._print_information_about_analyse_start()
        self._parse_the_structures()
        self._perform_analyse_and_parse_results()
        self._ensemble_analyse_records = EnsembleAnalyseRecords(
            reference= self._reference,
            ensemble= self._ensemble,
            analyse_results= self._analyse_results
            )
        ENMA._print_information_about_analyse_end()
        
    def display_eigenvalue_distribution(self) -> None:
        ENMA._plot_distribution_histogram(
            data= self.analyse_results.eigenvalues,
            x_label= "\nEigenvalue",
            y_label= "Frequency\n",
            title= f"Eigenvalue Distribution of Ensembles\n"
            )
    
    def display_eigenvalues_mode_by_mode(self) -> None:
        ENMA._plot_scree_plot(
            x_axis= np.array(range(1,len(self.analyse_results.eigenvalues)+1)),
            data= self.analyse_results.eigenvalues,
            x_label= "\nMode Number",
            y_label= "Eigenvalue\n",
            title= "Eigenvalues Mode by Mode\n"
            )
    
    def display_pca(self) -> None:
        ENMA._plot_pca(
            data= self.analyse_results.pca,
            labels= self.ensemble.names,
            x_label= "\nPrincipal Component 1",
            y_label= "Principal Component 2\n",
            title= "Principal Component Analysis of Ensembles Coordinates\n"
            )
    
    def _parse_the_structures(self) -> None:
        attributes = self._get_structures_attributes()
        self._reference = Reference(
            name= attributes["ids"]["reference"],
            structure= attributes["structures"]["reference"],
            structure_coordinate= attributes["coordinates"]["reference"]
            )
        self._ensemble = Ensemble(
            names= attributes["ids"]["ensembles"],
            structures= attributes["structures"]["ensembles"],
            structures_coordinates= attributes["coordinates"]["ensembles"]
            )
    
    def _get_structures_attributes(self) -> dict:
        structures_ids = self._get_structure_ids()
        structures = self._get_structures(
            reference_id= structures_ids["reference"],
            ensemble_ids= structures_ids["ensembles"]
            )
        aligned_structures = ENMA._get_aligned_structures(
            reference_structure= structures["reference"],
            ensemble_structures= structures["ensembles"]
            )
        coordinates = ENMA._get_coordinates(
            reference_structure= aligned_structures["reference"],
            ensemble_structures= aligned_structures["ensembles"]
            )
        return {
            "ids": structures_ids,
            "structures": aligned_structures,
            "coordinates": coordinates
            }
    
    def _get_structures(self, 
            reference_id: str, 
            ensemble_ids: List[str]
            ) -> dict:
        reference = ENMA._get_structure_from_pdb_file(
            structure_id= reference_id,
            file_path= self._reference_file_path
            )
        ensembles = []
        for id, file_path in zip(ensemble_ids, self._ensembles_files_paths):
            ensembles.append(
                ENMA._get_structure_from_pdb_file(
                    structure_id= id,
                    file_path= file_path
                    )
                )
        return {
            "reference": reference,
            "ensembles": ensembles
            }
    
    @staticmethod
    def _get_structure_from_pdb_file(
            structure_id: str, 
            file_path: str
            ) -> Structure:
        parser = PDBParser()
        return parser.get_structure(
            id= structure_id, 
            file= file_path
            )
    
    def _get_structure_ids(self) -> dict:
        reference = ENMA._extract_the_structure_id(
            file_path= self._reference_file_path
            )
        ensembles = []
        for file_path in self._ensembles_files_paths:
            ensembles.append(
                ENMA._extract_the_structure_id(
                    file_path= file_path
                    )
                )
        return {
            "reference": reference,
            "ensembles": ensembles
            }  
    
    @staticmethod
    def _extract_the_structure_id(
            file_path: str
            ) -> str:
        base_name = os.path.basename(file_path)
        return os.path.splitext(base_name)[0]
    
    @staticmethod
    def _get_aligned_structures(
            reference_structure: Structure, 
            ensemble_structures: List[Structure]
            ) -> dict:
        structures = {
            "reference": reference_structure,
            "ensembles": ensemble_structures
            }
        aligner = Superimposer()
        aligner.set_atoms(
            list(structures["reference"].get_atoms()), 
            list(structures["reference"].get_atoms())
            )
        for structure in structures["ensembles"]:
            aligner.apply(
                list(structure.get_atoms())
                )
        return structures
    
    @staticmethod
    def _get_coordinates(
            reference_structure: Structure, 
            ensemble_structures: List[Structure]
            ) -> dict:
        reference = ENMA._extract_coordinates(
            structure= reference_structure
            )
        ensembles = []
        for structure in ensemble_structures:
            ensembles.append(
                ENMA._extract_coordinates(
                    structure= structure
                    )
                )
        trimmed_ensembles = ENMA._trim_coordinates(
            structures_coordinates= ensembles
            )
        return {
            "reference": reference,
            "ensembles": trimmed_ensembles
            }
        
    @staticmethod
    def _extract_coordinates(
            structure: Structure
            ) -> np.ndarray:
        coordinates = []
        for atom in structure.get_atoms():
            coordinates.append(atom.get_coord())
        return np.array(coordinates)
    
    @staticmethod
    def _trim_coordinates(
            structures_coordinates: List[np.ndarray]
            ) -> list:
        minimum_length = ENMA._find_minimum_coordinate(
            structures_coordinates= structures_coordinates
            )
        trimmed_coordinates = []
        for structure_coordinates in structures_coordinates:
            trimmed_coordinates.append(
                structure_coordinates[:minimum_length]
                )
        return trimmed_coordinates
            
    @staticmethod
    def _find_minimum_coordinate(
            structures_coordinates: List[np.ndarray]
            ) -> int:
        lengths = []
        for structure_coordinates in structures_coordinates:
            lengths.append(
                len(structure_coordinates)
                )
        return min(lengths)
    
    def _perform_analyse_and_parse_results(self) -> None:
        results = ENMA._perform_analyse_and_get_results(
            reference_coordinates= self.reference.structure_coordinate,
            ensemble_coordinates= self.ensemble.structures_coordinates
            )
        self._analyse_results = AnalyseResults(
            mean= results["mean"],
            covariance_matrix= results["covariance_matrix"],
            distance_matrix= results["distance_matrix"],
            hessian_matrix= results["hessian_matrix"],
            eigenvalues= results["eigenvalues"],
            eigenvectors= results["eigenvectors"],
            pca= results["pca"]
            )
    
    @staticmethod
    def _perform_analyse_and_get_results(
            reference_coordinates: np.ndarray, 
            ensemble_coordinates: List[np.ndarray]
            ) -> dict:
        mean = ENMA._calculate_mean(
            ensemble_coordinates= np.array(ensemble_coordinates)
            )
        covariance_matrix = ENMA._calculate_covariance(
            ensemble_coordinates= np.array(ensemble_coordinates),
            ensemble_mean= mean
            )
        distance_matrix = ENMA._calculate_distance_matrix(
            reference_coordinates= reference_coordinates,
            ensemble_coordinates= ensemble_coordinates
            )
        hessian_matrix = ENMA._calculate_hessian_matrix(
            distance_matrix= distance_matrix,
            covariance_matrix= covariance_matrix
            )
        eigens = ENMA._calculate_eigens(
            hessian_matrix= hessian_matrix
            )
        pca = ENMA._perform_pca(
            ensemble_coordinates= ensemble_coordinates
            )
        return {
            "mean": mean,
            "covariance_matrix": covariance_matrix,
            "distance_matrix": distance_matrix,
            "hessian_matrix": hessian_matrix,
            "eigenvalues": eigens["eigenvalues"],
            "eigenvectors": eigens["eigenvectors"],
            "pca": pca
            }
    
    @staticmethod
    def _calculate_mean(
            ensemble_coordinates: np.ndarray
            ) -> np.ndarray:
        return np.mean(
            a= ensemble_coordinates,
            axis= 0
            )
    
    @staticmethod
    def _calculate_covariance(
            ensemble_coordinates: np.ndarray, 
            ensemble_mean: np.ndarray
            ) -> np.ndarray:
        fluctuations = ensemble_coordinates - ensemble_mean[np.newaxis,:,:]
        outers = []
        for fluctuation in fluctuations:
            outers.append(
                np.outer(
                    a=fluctuation.flatten(), 
                    b=fluctuation.flatten()
                    )
                )
        return np.mean(
            a= outers,
            axis= 0
            )
    
    @staticmethod
    def _calculate_distance_matrix(
            reference_coordinates: np.ndarray, 
            ensemble_coordinates: List[np.ndarray]
            ) -> np.ndarray:
        size = ensemble_coordinates[0].shape[0]
        return distance_matrix(
            reference_coordinates[:size],
            reference_coordinates[:size]
            )
    
    @staticmethod
    def _calculate_hessian_matrix(
            distance_matrix: np.ndarray, 
            covariance_matrix: np.ndarray
            ) -> np.ndarray:
        hessian_matrix = np.zeros_like(covariance_matrix)
        size = len(distance_matrix)
        for row in range(size):
            for column in range(size):
                if row != column and distance_matrix[row,column] != 0:
                    covariance = covariance_matrix[row,column]
                    distance = distance_matrix[row,column]**3
                    hessian_matrix[row,column] = covariance / distance
        return hessian_matrix
    
    @staticmethod
    def _calculate_eigens(
            hessian_matrix: np.ndarray
            ) -> dict:
        hessian_sparse = csr_matrix(hessian_matrix)
        eigenvalues, eigenvectors = eigh(hessian_sparse.todense())
        return {
            "eigenvalues": eigenvalues,
            "eigenvectors": eigenvectors
            }
        
    @staticmethod
    def _perform_pca(
            ensemble_coordinates: List[np.ndarray]
            ) -> np.ndarray:
        coordinates_array = np.array(ensemble_coordinates)
        coordinates_size = len(ensemble_coordinates)
        pca_data = coordinates_array.reshape(coordinates_size, -1)
        pca = PCA(n_components=2)
        return pca.fit_transform(pca_data)
    
    @staticmethod
    def _plot_scree_plot(
            x_axis: np.ndarray, 
            data: np.ndarray, 
            x_label: str, 
            y_label: str, 
            title: str
            ) -> None:
        plotter.figure(figsize=(12,8))
        plotter.plot(
            x_axis,
            data,
            "-o"
            )
        plotter.xlabel(x_label)
        plotter.ylabel(y_label)
        plotter.title(title)
        plotter.show()
        
    @staticmethod
    def _plot_distribution_histogram(data: np.ndarray, x_label: str, y_label: str, title: str) -> None:
        # Create the histogram
        fig = go.Figure(data=[go.Histogram(x=data, nbinsx=50)])
        
        # Update layout for labels and title
        fig.update_layout(
            title=title,
            xaxis_title=x_label,
            yaxis_title=y_label,
            bargap=0.2,  # Gap between bars of adjacent location coordinates
        )
        
        # Show the plot
        fig.show()
    @staticmethod
    def _plot_pca(
            data: np.ndarray, 
            labels: List[str], 
            x_label: str, 
            y_label: str, 
            title: str
            ) -> None:
        plotter.figure(figsize=(12,8))
        colors = plotter.get_cmap("tab10")
        for index, label in enumerate(labels):
            pc1 = data[index, 0]
            pc2 = data[index, 1]
            plotter.scatter(
                pc1,
                pc2,
                color= colors(index),
                label= label
                )
        plotter.xlabel(x_label)
        plotter.ylabel(y_label)
        plotter.title(title)
        plotter.legend(loc="best")
        plotter.show()
        
    @staticmethod
    def _print_information_about_analyse_start() -> None:
        print("Ensemble normal mode analyse is started...")
        print("It can be take several seconds...")
    
    @staticmethod
    def _print_information_about_analyse_end() -> None:
        print("Ensemble normal mode analyse is ended.")
        print("Anlyse results are displayable.")





class Blast():
    DOWNLOAD_PDB_FILE_URL: str = "https://files.rcsb.org/download"
    
    def __init__(self, 
            sequence: str, 
            accession: str
            ) -> None:
        warnings.filterwarnings("ignore")
        self._sequence = sequence
        self._accession = accession
        self._similar_structures = None    
        
    @property
    def sequence(self) -> str:
        return self._sequence
    
    @sequence.setter 
    def sequence(self, value) -> None:
        self._sequence = value
        
    @property
    def accession(self) -> str:
        return self._accession
    
    @accession.setter
    def accession(self, value) -> None:
        self._accession = value
        
    @property
    def similar_strucutures(self) -> BlastRecords:
        return self._similar_structures
    
    def download_files(self, 
            output_folder_path: str = None, 
            similar_structures_download_status: bool = True, 
            reference_structure_download_status: bool = True
            ) -> None:
        if output_folder_path is None:
            output_folder_path = os.getcwd()
        if similar_structures_download_status:
            for accession in self._get_accession_list_of_similar_structures():
                Blast._download_a_pdb_file(
                    output_folder_path= output_folder_path,
                    accession= accession
                    )
        if reference_structure_download_status:
            Blast._download_a_pdb_file(
                output_folder_path= output_folder_path,
                accession= self.accession
                )
    
    def show_alignment(self) -> None:
        if not self.similar_strucutures is None:
            self._plot_alignments(
                data= self._get_similar_structure_alignment_data(),
                x_label= "\nPositions",
                y_label= "Accession of Proteins\n",
                title= f"Sequence Alignment Plot of {self.accession}\n"
                )
    
    
    def make_blast_search(self, 
            size: int = 20, 
            program: str = "blastp", 
            database: str = "pdb"
            ) -> None:
    
        blast_datas = self._get_blast_datas_from_database(
            program= program,
            database= database,
            hitlist_size= size + 1
            )
        
        self._similar_structures = self._parse(
            datas= blast_datas
            )
                
    @staticmethod
    def _download_a_pdb_file(
            output_folder_path: str, 
            accession: str
            ) -> None:
        output_file_name = f"{accession}.pdb"
        url = os.path.join(
            Blast.DOWNLOAD_PDB_FILE_URL,
            output_file_name
            )
        output_path = os.path.join(
            output_folder_path,
            output_file_name
            )
        Blast._download_a_file(
            url= url,
            output_path= output_path
            )
            
    def _get_blast_datas_from_database(self, 
            program: str, 
            database: str, 
            hitlist_size: int
            ) -> object:
        datas = NCBIWWW.qblast(
            program, 
            database, 
            self.sequence, 
            hitlist_size = hitlist_size
            )
        return NCBIXML.read(datas)

    
    def _get_accession_list_of_similar_structures(self) -> list:
        accession_list = []
        if not self.similar_strucutures is None:
            for record in self.similar_strucutures.records:
                accession_list.append(str(record.accession))
            print(accession_list)
            return accession_list                
    
    
    @staticmethod
    def _download_a_file(
            url: str, 
            output_path: str
            ) -> None:
        response = requests.get(url= url)
        if response.status_code == 200:
            with open(output_path, "wb") as file:
                file.write(response.content)
    
    
    
    def _get_parsed_blast_records(self, 
            alignments: Iterable
            ) -> list:
        blast_records = []
        for alignment in alignments:
            blast_record = self._get_parsed_blast_record(
                alignment= alignment
                )
            if not blast_record is None:
                blast_records.append(blast_record)
        return blast_records

    def _parse(self, 
            datas: object
            ) -> BlastRecords:
        records = self._get_parsed_blast_records(
            alignments= datas.alignments
            )
        data = {
            "reference_sequence": {
                "accession": self.accession,
                "length": len(self.sequence),
                "sequence": self.sequence
                },
            "records": records
            }
        return BlastRecords(**data)
    def _is_the_reference_protein(self, 
            accession: str
            ) -> bool:
        return accession.startswith(self.accession)
    
    def _get_parsed_blast_record(self, 
            alignment: object
            ) -> BlastRecord:
        if not self._is_the_reference_protein(alignment.accession):
            hsps = Blast._get_parsed_high_scoring_pairs(
                high_scoring_pairs= alignment.hsps
                )
            underscore_index = str(alignment.accession).find("_")
            record = {
                "accession": alignment.accession[:underscore_index],
                "chain_code": alignment.accession[underscore_index+1:],
                "length": alignment.length,
                "high_scoring_pairs": hsps
                }
            return BlastRecord(**record)
    
    
    
    
    
    def _get_similar_structure_alignment_data(self) -> dict:
        alignments = {}
        minimum_position = float("inf")
        maximum_position = float("-inf")
        for blast_record in self.similar_strucutures.records:
            fragments = []
            for high_scoring_pair in blast_record.high_scoring_pairs:
                positions = {
                    "start_position": high_scoring_pair.query_start_position,
                    "end_position": high_scoring_pair.query_end_position
                    }
                if positions["end_position"] > maximum_position:
                    maximum_position = positions["end_position"]
                if positions["start_position"] < minimum_position:
                    minimum_position = positions["start_position"]
                fragments.append(positions)
            accession = blast_record.accession
            chain_code = blast_record.chain_code
            alignments[f"{accession} (Chain {chain_code})"] = fragments
        reference = {
            "accession": self.accession,
            "start_position": minimum_position,
            "end_position": maximum_position
            }
        return {
            "reference": reference,
            "alignments": alignments
            }
    @staticmethod
    def _get_parsed_high_scoring_pairs(
            high_scoring_pairs: Iterable
            ) -> list:
        parsed_high_scoring_pairs = []
        for high_scoring_pair in high_scoring_pairs:
            parsed_high_scoring_pairs.append(
                HighScoringPair(
                    score= high_scoring_pair.score,
                    bits= high_scoring_pair.bits,
                    expect= high_scoring_pair.expect,
                    positives= high_scoring_pair.positives,
                    gaps= high_scoring_pair.gaps,
                    alignment_length= high_scoring_pair.align_length,
                    query_sequence= high_scoring_pair.query,
                    match_sequence= high_scoring_pair.match,
                    subject_sequence= high_scoring_pair.sbjct,
                    query_start_position= high_scoring_pair.query_start,
                    query_end_position= high_scoring_pair.query_end,
                    subject_start_position= high_scoring_pair.sbjct_start,
                    subject_end_position= high_scoring_pair.sbjct_end
                    )
                )
        return parsed_high_scoring_pairs

    def _plot_alignments(self,data: dict, x_label: str, y_label: str, title: str):
        fig = go.Figure()

        # Add horizontal lines for each alignment
        for alignment in data["alignments"].keys():
            fragments = data["alignments"][alignment]
            for fragment in fragments:
                fig.add_trace(go.Scatter(
                    x=[fragment["start_position"], fragment["end_position"]],
                    y=[alignment, alignment],
                    mode='lines',
                    line=dict(color='blue', width=5),
                    name=alignment
                ))

        # Add horizontal line for the reference
        fig.add_trace(go.Scatter(
            x=[data["reference"]["start_position"], data["reference"]["end_position"]],
            y=[data["reference"]["accession"], data["reference"]["accession"]],
            mode='lines',
            line=dict(color='red', width=8),
            name=data["reference"]["accession"]
        ))

        # Set labels and title
        fig.update_layout(
            xaxis_title=x_label,
            yaxis_title=y_label,
            title=title
        )

        # Show plot
        fig.show()

    
    
    
    