import warnings
warnings.filterwarnings("ignore")

from Bio import Entrez, SeqIO

class ProteinSynthesis():
  def __init__(self):
    self.genetic_code = [[["F", "F", "L", "L"], ["S", "S", "S", "S"], ["Y", "Y", " ", " "], ["C", "C", " ", "W"],],
                          [["L", "L", "L", "L"], ["P", "P", "P", "P"], ["H", "H", "Q", "Q"], ["R", "R", "R", "R"],],
                          [["I", "I", "I", "M"], ["T", "T", "T", "T"], ["N", "N", "K", "K"], ["S", "S", "R", "R"],],
                          [["V", "V", "V", "V"], ["A", "A", "A", "A"], ["D", "D", "E", "E"], ["G", "G", "G", "G"],]]

    self.base_encode = {
        "U": 0,
        "C": 1,
        "A": 2,
        "G": 3
    }

  def transcribed(self, sequence):
    new_sequence = ""
    for base in sequence:
      if base == "A":
        new_sequence += "U" # A is paired with T but later T change to U
      elif base == "T":
        new_sequence += "A"
      elif base == "C":
        new_sequence += "G"
      elif base == "G":
        new_sequence += "C"
    return new_sequence

  def translation(self, sequence):
    new_sequence = ""
    seq_length = len(sequence)
    throw_away = seq_length % 3
    for ind in range(0, seq_length - throw_away, 3):
      new_sequence += self.genetic_code[self.base_encode[sequence[ind]]][self.base_encode[sequence[ind + 1]]][self.base_encode[sequence[ind + 2]]]
    return new_sequence

  def sintesis(self, dna_sequence):
    transcribed_sequence = self.transcribed(dna_sequence)
    translation_sequence = self.translation(transcribed_sequence)
    return translation_sequence

def get_data_from_ncbi():
  Entrez.email = "venus.angela.kurniawan@mail.ugm.ac.id"
  accession_number = "NC_045512.2"
  # accession_number = "MW400961.1"
  handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
  genome_data = SeqIO.read(handle, "genbank")
  handle.close()
  return genome_data
 
protein_sintesis = ProteinSynthesis()
first_sequence = "GGATGCCAATG"
second_sequence = "TCGGTGAATCTGTTTGAT"
original_sequence = get_data_from_ncbi().seq

first_transcribed = protein_sintesis.transcribed(first_sequence)
first_translation = protein_sintesis.translation(first_transcribed)

print("\nSequence:", first_sequence)
print("Transcibed to:", first_transcribed)
print("Translate to:", first_translation)

second_transcribed = protein_sintesis.transcribed(second_sequence)
second_translation = protein_sintesis.translation(second_transcribed)

print("\nSequence:", second_sequence)
print("Transribed to:", second_transcribed)
print("Translate to:", second_translation)

amino_acid = protein_sintesis.sintesis(original_sequence)
amino_acid_check = original_sequence.complement().transcribe().translate(stop_symbol=" ")
print("\nOriginal Sequence:", original_sequence[:50])
print("Amino Acid (Code):", amino_acid[:50])
print("Amino Acid (Check):", amino_acid_check[:50])
if amino_acid == amino_acid_check:
  print("Two sequences are equal")
else:
  print("Two sequences are not equal")