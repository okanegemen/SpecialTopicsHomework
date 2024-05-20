from dataclasses import dataclass
from typing import List
@dataclass
class AminoAcidInfo():
    name: str
    one_letter_code: str
    three_letter_code: str


class AminoAcid():
    alanine: AminoAcidInfo = AminoAcidInfo(
        name= "Alanine", 
        one_letter_code= "A", 
        three_letter_code= "Ala"
        )
    arginine: AminoAcidInfo = AminoAcidInfo(
        name= "Arginine", 
        one_letter_code= "R", 
        three_letter_code= "Arg"
        )
    asparagine: AminoAcidInfo = AminoAcidInfo(
        name= "Asparagine",
        one_letter_code= "N",
        three_letter_code= "Asn"
        )
    aspartate: AminoAcidInfo = AminoAcidInfo(
        name= "Aspartate",
        one_letter_code= "D",
        three_letter_code= "Asp"
        )
    cysteine: AminoAcidInfo = AminoAcidInfo(
        name= "Cysteine",
        one_letter_code= "C",
        three_letter_code= "Cys"
        )
    glutamate: AminoAcidInfo = AminoAcidInfo(
        name= "Glutamate",
        one_letter_code= "E",
        three_letter_code= "Glu"
        )
    glutamine: AminoAcidInfo = AminoAcidInfo(
        name= "Glutamine",
        one_letter_code= "Q",
        three_letter_code= "Gln"
        )
    glycine: AminoAcidInfo = AminoAcidInfo(
        name= "Glycine",
        one_letter_code= "G",
        three_letter_code= "Gly"
        )
    histidine: AminoAcidInfo = AminoAcidInfo(
        name= "Histidine",
        one_letter_code= "H",
        three_letter_code= "His"
        )
    isoleucine: AminoAcidInfo = AminoAcidInfo(
        name= "Isoleucine",
        one_letter_code= "I",
        three_letter_code= "Ile"
        )
    leucine: AminoAcidInfo = AminoAcidInfo(
        name= "Leucine",
        one_letter_code= "L",
        three_letter_code= "Leu"
        )
    lysine: AminoAcidInfo = AminoAcidInfo(
        name= "Lysine",
        one_letter_code= "K",
        three_letter_code= "Lys"
        )
    methionine: AminoAcidInfo = AminoAcidInfo(
        name= "Methionine",
        one_letter_code= "M",
        three_letter_code= "Met"
        )
    phenylalanine: AminoAcidInfo = AminoAcidInfo(
        name= "Phenylalanine",
        one_letter_code= "F",
        three_letter_code= "Phe"
        )
    proline: AminoAcidInfo = AminoAcidInfo(
        name= "Proline",
        one_letter_code= "P",
        three_letter_code= "Pro"
        )
    serine: AminoAcidInfo = AminoAcidInfo(
        name= "Serine",
        one_letter_code= "S",
        three_letter_code= "Ser"
        )
    threonine: AminoAcidInfo = AminoAcidInfo(
        name= "Threonine",
        one_letter_code= "T",
        three_letter_code= "Thr"
        )
    tryptophan: AminoAcidInfo = AminoAcidInfo(
        name= "Tryptophan",
        one_letter_code= "W",
        three_letter_code= "Trp"
        )
    tyrosine: AminoAcidInfo = AminoAcidInfo(
        name= "Tyrosine",
        one_letter_code= "Y",
        three_letter_code= "Tyr"
        )
    valine: AminoAcidInfo = AminoAcidInfo(
        name= "Valine",
        one_letter_code= "V",
        three_letter_code= "Val"
        )
    
class AminoAcids():
    amino_acids: dict = {
        "Alanine": AminoAcid.alanine,
        "Arginine": AminoAcid.arginine,
        "Asparagine": AminoAcid.asparagine,
        "Aspartate": AminoAcid.aspartate,
        "Cysteine": AminoAcid.cysteine,
        "Glutamate": AminoAcid.glutamine,
        "Glutamine": AminoAcid.glutamine,
        "Glycine": AminoAcid.glycine,
        "Histidine": AminoAcid.histidine,
        "Isoleucine": AminoAcid.isoleucine,
        "Leucine": AminoAcid.leucine,
        "Lysine": AminoAcid.lysine,
        "Methionine": AminoAcid.methionine,
        "Phenylalanine": AminoAcid.phenylalanine,
        "Proline": AminoAcid.proline,
        "Serine": AminoAcid.serine,
        "Threonine": AminoAcid.threonine,
        "Tryptophan": AminoAcid.tryptophan,
        "Tyrosine": AminoAcid.tyrosine,
        "Valine": AminoAcid.valine
        }
    amino_acids_1_letter_code = {
        "A": AminoAcid.alanine,
        "R": AminoAcid.arginine,
        "N": AminoAcid.asparagine,
        "D": AminoAcid.aspartate,
        "C": AminoAcid.cysteine,
        "E": AminoAcid.glutamine,
        "Q": AminoAcid.glutamine,
        "G": AminoAcid.glycine,
        "H": AminoAcid.histidine,
        "I": AminoAcid.isoleucine,
        "L": AminoAcid.leucine,
        "K": AminoAcid.lysine,
        "M": AminoAcid.methionine,
        "F": AminoAcid.phenylalanine,
        "P": AminoAcid.proline,
        "S": AminoAcid.serine,
        "T": AminoAcid.threonine,
        "W": AminoAcid.tryptophan,
        "Y": AminoAcid.tyrosine,
        "V": AminoAcid.valine
        }
    amino_acids_3_letter_code = {
        "Ala": AminoAcid.alanine,
        "Arg": AminoAcid.arginine,
        "Asn": AminoAcid.asparagine,
        "Asp": AminoAcid.aspartate,
        "Cys": AminoAcid.cysteine,
        "Glu": AminoAcid.glutamine,
        "Gln": AminoAcid.glutamine,
        "Gly": AminoAcid.glycine,
        "His": AminoAcid.histidine,
        "Ile": AminoAcid.isoleucine,
        "Lue": AminoAcid.leucine,
        "Lys": AminoAcid.lysine,
        "Met": AminoAcid.methionine,
        "Phe": AminoAcid.phenylalanine,
        "Pro": AminoAcid.proline,
        "Ser": AminoAcid.serine,
        "Thr": AminoAcid.threonine,
        "Trp": AminoAcid.tryptophan,
        "Tyr": AminoAcid.tyrosine,
        "Val": AminoAcid.valine
        }


@dataclass
class HighScoringPair():
    score: float
    bits: float
    expect: float
    positives: int
    gaps: int
    alignment_length: int
    query_sequence: str
    match_sequence: str
    subject_sequence: str
    query_start_position: int
    query_end_position: int
    subject_start_position : int
    subject_end_position: int

@dataclass
class BlastRecord():
    accession: str
    chain_code: str
    length: int
    high_scoring_pairs : List[HighScoringPair]
    
@dataclass
class ReferenceSequence():
    accession: str
    length: str
    sequence: str
    
@dataclass
class BlastRecords():
    reference_sequence: ReferenceSequence
    records: List[BlastRecord]



from dataclasses import dataclass
from typing import List
from numpy import ndarray
from Bio.PDB.Structure import Structure

@dataclass
class Reference():
    name: str
    structure: Structure
    structure_coordinate: ndarray

@dataclass
class Ensemble():
    names: List[str]
    structures: List[Structure]
    structures_coordinates: List[ndarray]

@dataclass
class AnalyseResults():
    mean: ndarray
    covariance_matrix: ndarray
    distance_matrix: ndarray
    hessian_matrix: ndarray
    eigenvalues: ndarray
    eigenvectors: ndarray
    pca: ndarray

@dataclass 
class EnsembleAnalyseRecords():
    reference: Reference
    ensemble: Ensemble
    analyse_results: AnalyseResults