import os, glob
from Bio import SeqIO
from Bio import Entrez
import pandas as pd 
    
name2 = []
trna_num = []
trna_bp_pos = []
intron_start = []     # Declare lists to store information of output files
intron_end = []
score2 = [] 
note = []
tax = []
org = []

path = 'candidatus_complete_trnascan_out' # Fill in directory name with tRNAScan-Se output and structure files
for filename in glob.glob(os.path.join(path, '*.out')):
   with open(os.path.join(os.getcwd(), filename), 'r') as f: # Open files that end with .out in directory in readonly mode
        print(filename) # prints names of files being processed
        count = 0
        for lines in open(filename): # Loops through the files that end with .out 
            if lines.startswith("Sequence") or lines.startswith("Name") or lines.startswith("--------"):
                next
            else: # Takes information out of files and puts them in lists
                name2.append(lines.split()[0])
                trna_num.append(lines.split()[1])
                trna_bp_pos.append(lines.split()[2] + "-" + lines.split()[3])
                intron_start.append(lines.split()[6])
                intron_end.append(lines.split()[7])
                score2.append(lines.split()[8])
                try:
                    note.append(lines.split()[9])
                except:
                    note.append("")
                if count < 1:
                    count+=1
                    try: # connecting with Entrez to search for taxonomy in NCBI databases
                        Entrez.email = "Your.Name.Here@example.org"
                        handle = Entrez.efetch(db="nucleotide", id=lines.split()[0], rettype="gb", retmode="text")
                        x = SeqIO.read(handle, 'genbank')
                        nameO = x.annotations['organism']
                        nameT = x.annotations['taxonomy']
                        org.append(x.annotations['organism'])
                        tax.append(x.annotations['taxonomy'])
                    except:
                        org.append('')
                        tax.append('')
                else: 
                    org.append(nameO)
                    tax.append(nameT)
                    
                    
typ = []
anticodon = []
anticodon_pos = []     # Declare lists to store information of output files
bp_pos = []
score = []
trna_length = []
name = []

path = 'candidatus_complete_trnascan_out' # Fill in directory name with tRNAScan-Se output and structure files
for filename in glob.glob(os.path.join(path, '*.struct')):
   with open(os.path.join(os.getcwd(), filename), 'r') as f: # Open files that end with .out in directory in readonly mode
        print(filename) # prints names of files being processed
        for lines in open(filename): # Loops through the files that end with .struct
            if lines.startswith("Type:"): # Takes information out of files and puts them in lists
                typ.append(lines.split()[1])
                anticodon.append(lines.split()[3])
                anticodon_pos.append(lines.split()[5])
                bp_pos.append(lines.split()[6])
                score.append(lines.split()[8])
            if lines.endswith("bp\n"): # Takes information out of files and puts them in lists
                name.append(lines.split(".t")[0])
                trna_length.append(lines.split()[3])               
                
                
                
arm = []
extra = []
for c, a, t in zip(typ, anticodon, org): # for loop to classify the tRNAs with extra arms
    if "Leu" in c:
        if "TAA" or "CAA" in a and "Frankia datiscae" in t:
            extra.append("No")
            arm.append("Class 2")
        else:
            extra.append("Yes")
            arm.append("Class 2")
    elif "Ser" in c:
        arm.append("Class 1")
        extra.append("Yes")
    elif "SeC" in c:
        extra.append("Yes")
        arm.append("")
    elif "Met" in c or "Val" in c or "Ile" in c or "Cys" in c or "Arg" in c or "Glu" in c or "Gln" in c or "Lys" in c or "Tyr" in c or "Trp" in c:
        extra.append("No")
        arm.append("Class 1")
    else:
        arm.append("Class 2")
        extra.append("No")
        
 
mydict = {'Name': name2, 'TRNA number': trna_num, 'TRNA genome position': trna_bp_pos, 'Intron begin': intron_start,
          'Intron end': intron_end, 'Score': score2, 'Note': note, 'Organism scientific name': org,
          'Taxonomic classification': tax,
          'TRNA type': typ, 'TRNA length': trna_length, 'Anticodon': anticodon,
          'Anticodon position': anticodon_pos, 'Anticodon genome position': bp_pos, 
          'Score': score,"tRNA extra arm": extra, "tRNA class": arm}

df = pd.DataFrame(mydict)

df.to_excel('tRNAScan-Se_output_tabulated_with_extra_arm_classifier_complete_bac.xlsx')
