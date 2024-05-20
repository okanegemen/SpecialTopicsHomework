from implementation import *
from Bio import SeqIO

cur_dir = os.getcwd()

sequence_3ASW = ""
fasta_file_path = f"{cur_dir}/HomeWork/_fasta_files/3ASW.fasta"


size = 10 if fasta_file_path == f"{cur_dir}/HomeWork/_fasta_files/3ASW.fasta" else 9
with open(fasta_file_path, "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_3ASW += str(record.seq)

processor = AAComposition(
    sequence= sequence_3ASW
    )

processor.display_amino_acid_composition_profile(
    order_type= "a", 
    anotation_state= True
    )

processor.display_amino_acid_histogram_profile()

processor_blast = Blast(
    sequence= sequence_3ASW,
    accession= "3ASW"
    )
print("Beginned!")
processor_blast.make_blast_search(
    size = size,
    program= "blastp",
    database= "pdb"
    )

processor_blast.show_alignment()

processor_blast.download_files(
    output_folder_path= "_pdb_files",
    similar_structures_download_status= True,
    reference_structure_download_status= True
    )

reference_file_path = f"{cur_dir}/HomeWork/_pdb_files"
ensemble_files_paths = (
    "_pdb_files/3AT0.pdb",
    "_pdb_files/3AU0.pdb",
    "_pdb_files/4F1Z.pdb",
    "_pdb_files/4F27.pdb",
    "_pdb_files/4JE0.pdb",
    "_pdb_files/4JDZ.pdb",
    "_pdb_files/5CF3.pdb",
    "_pdb_files/1N67.pdb",
    "_pdb_files/2VR3.pdb"
    )

processor_enma = ENMA(
    reference_file_path= reference_file_path,
    ensembles_files_paths= ensemble_files_paths
    )

processor_enma.perform_ensemble_normal_mode_analyse()

processor_enma.display_eigenvalue_distribution()

processor_enma.display_eigenvalues_mode_by_mode()




