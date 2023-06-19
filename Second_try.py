import subprocess
import sys

def recupererKmers(nomfichier,k) :
    nomSortie = nomfichier.split(".")[0] + "_kmers.fasta"
    result = subprocess.run(['gkampi','-P','-k',str(k),'-o',nomSortie,'-f','-Q',nomfichier],capture_output=True, text=True,universal_newlines=True)
    if result.returncode != 0 :
        print("Erreur dans la récuperation du fichier d'entrée")
        return
    #else :
        #print("Les kmers sont stockées dans le fichier "+str(nomSortie))

    kmer = {}
    sequences = []
    position = []
    with open (nomSortie) as fasta_out :
        for sequence in fasta_out :
            if sequence[0] != ">" :
                sequences += [sequence[:(len(sequence)-1)]]
            else :
                position += [sequence[2:(len(sequence)-1)]]


    for entier in range(len(sequences)) :
        if sequences[entier] not in kmer :
            kmer[sequences[entier]] = []
            kmer[sequences[entier]] = [int(position[entier])]
        else :
            kmer[sequences[entier]] += [int(position[entier])]
            

    
    #print("Dictionnaire realise correctement") 
    return kmer

# Function that forms a dictionary with the sequences of k-mers and their different positions, also creates a fasta file with the k-mers from the assemblies
# Input: String of the filename used to extract the k-mers, int the required K to form the k-mers
# Output: Dictionary of all the generated k-mers and their relative positions
# Side effect: Creates a fasta file containing all the k-mers in the current directory
# Comment: The position currently provided by gkampi is the absolute position of the k-mer in all the sequences (does not take into account scaffold changes and considers everything as a single sequence, so position 0 is for scaffold 1 only, and no other scaffolds. Additionally, it does not consider 'N' characters in its k-mer and position calculations). To add the positions relative to the scaffold, additional modifications are required.
# Comment: If the output files already exist, they will be overwritten.

kmer_ref = recupererKmers("Test_sequences/rice_ref.fasta",100)
kmer_ass = recupererKmers("Test_sequences/rice_ass.fasta",100)

def extract_scaffold(file_name : str) -> dict:
    liste_scaffold = []
    scaffold_courant = ""

    with open(file_name) as fasta_in:
        for ligne in fasta_in:
            if ligne[0] == ">":  
                if scaffold_courant != "":
                    liste_scaffold.append(scaffold_courant)
                    scaffold_courant = ""
            else:
                scaffold_courant += ligne.strip()

        if scaffold_courant != "":
            liste_scaffold.append(scaffold_courant)

    return liste_scaffold

# Function that extract the scaffold of a fasta file 
# Input : string : name of the fasta 
# Output : list : list of the scaffold present in the fasta file 

scaffold_assemblage = extract_scaffold("Test_sequences/rice_ass.fasta")
scaffold_reference = extract_scaffold("Test_sequences/rice_ref.fasta")

def unique (dictionnaire_kmer,kmer) :
    return len(dictionnaire_kmer[kmer]) == 1
# Function that test if a sequence is unique in the dictionnary ( form = {"sequence" : [position(s)], ...})
# Input : dict : dictionnary of all the sequences to compare  ( form = {"sequence" : [position(s)], ...}) , string : sequence to campare the dictionnary with 
# Output : bool : True  = the sequence is alone in the dictionnary 

def trouver_ancre (dico_kmer1 : dict,dico_kmer2 : dict) -> tuple[str,int,int] :

    for kmer1 in dico_kmer1.keys() :
        for kmer2 in dico_kmer2.keys() :
            if kmer1 == kmer2 and unique(dico_kmer1,kmer1) and unique(dico_kmer2,kmer2) :
                return kmer1,dico_kmer1[kmer1],dico_kmer2[kmer2]
    
    print("Aucune ancre n'a ete trouvee ",file = sys.stderr)
    return None, None, None
# Function that finds the first unique sequence present in both dictionnary ( form = {"sequence" : [position(s)],...})
# Input: dict : first dictionnary to compare , dict : second dictionnary to compare, int : k-mer length
# Output: If a unique sequence anchor is found: returns the k-mer, position of the first nucleotide in the first sequence, position of the first nucleotide in the second sequence
# Output: Otherwise, returns None, None, None

ancre = trouver_ancre(recupererKmers("Test_sequences/rice_ref.fasta",100),recupererKmers("Test_sequences/rice_ass.fasta",100))

def extraire_entete(filename : str,position : int) -> str :
    id = []
    numero_id= 0
    pos =0
    with open(filename) as fichier :
        for scaffold in fichier :
            if scaffold[0] != ">" and pos <= position :
                numero_id += 1
                pos += len(scaffold)-1
            if scaffold[0] == ">" :
                id += [scaffold[1:len(scaffold)-1]]
    if len(id) != 1 :
        return id[numero_id-1]
    else :
        return id[0]            
# Function that extracts the header of the scaffold containing the requested position
# Input: string: Name of the fasta file, int: Position of the nucleotide in the assembly
# Output: string: Header of the scaffold containing the requested position

tete = extraire_entete("Test_sequences/rice_ref.fasta",int(ancre[1][0]))

def trouve_scaffold(liste_scaffold : list, position : int) -> str:
    pos = 0
    
    for scaffold in liste_scaffold:
        pos += len(scaffold)
        if pos > position:
            return scaffold
    return None  
# Function that find the scaffold containing the given position ()
# Input : list : List of scaffold , int : Position searched 
# Output : Return the scaffold or None if the position isn't in the scaffold's list

scaffold_ancre_assemblage = trouve_scaffold(scaffold_assemblage,int(ancre[2][0])) 
scaffold_ancre_reference = trouve_scaffold(scaffold_reference,int(ancre[1][0]))

def position_in_scaffold (scaffold : str,kmer : str) -> int :
    itterateur =0
    
    while itterateur + len(kmer) < len(scaffold) :
        if scaffold[itterateur : itterateur+len(kmer)] == kmer :
            return itterateur
        itterateur+=1
    print("Position non trouvé",file = sys.stderr)
    return None
# Function that returns the position of the kmer in the given scaffold if it is present 
# Input: string : scaffold that contains the kmer , string  : kmer searched 
# Output: int position of the start of the kmer in the scaffold, or None

position_assemblage = position_in_scaffold(scaffold_ancre_assemblage,ancre[0])
position_reference = position_in_scaffold(scaffold_ancre_reference,ancre[0])

def next_kmer_overlap (position : int ,value_overlap : int, kmer_dictionnary : dict ) -> tuple[str,int]:
    for kmer in kmer_dictionnary.keys() : 
        if position + value_overlap == kmer_dictionnary[kmer][0] :
            return  kmer,kmer_dictionnary[kmer][0]

    print("Pas de prochain kmer trouvé", file = sys.stderr)
    return None,None


def next_kmer_no_overlap (position : int, kmer_length : int , kmer_dictionnary : dict) -> str :
    for kmer in kmer_dictionnary.keys() :
        if position + kmer_length == kmer_dictionnary[kmer][0] : 
            return kmer
    
    print("Pas de prochain kmer trouvé" , file=sys.stderr)
    return None

def extend_ancor (ancor : str,ancor_position_in_scaffold1 : int,ancor_position_in_scaffold2 : int ,kmer_dictionnary_scaffold1 : dict,kmer_dictionnary_scaffold2 : dict,scaffold1 : str,scaffold2 : str)  : 
    kmer_length = len(ancor)
    position2 = ancor_position_in_scaffold2 + kmer_length
    position1 = ancor_position_in_scaffold1 + kmer_length
    print(position1,position2)
    kmer_of_scaffolds = {}
    while position1 + kmer_length < len(scaffold1) and position2 + kmer_length < len(scaffold2) : # Tant que les kmer non chevauchant sont possible on regarde ceux-ci
        kmer_scaffold1 = next_kmer_no_overlap (position1,kmer_length,kmer_dictionnary_scaffold1) 
        kmer_scaffold2= next_kmer_no_overlap (position2,kmer_length,kmer_dictionnary_scaffold2) 
        print("zizi")
    
        if kmer_scaffold1 == kmer_scaffold2 and kmer_scaffold1 == scaffold1: #Si les kmer non chevauchant sont égaux on les enregiste et on passe au suivant
            if kmer_scaffold1 in kmer_of_scaffolds :
                kmer_of_scaffolds[kmer_scaffold1] += [[position1,position2]]
                position1 += kmer_length
                position2 += kmer_length
            else :
                kmer_of_scaffolds[kmer_scaffold1] = []
                kmer_of_scaffolds[kmer_scaffold1] += [[position1,position2]]
                position1 += kmer_length
                position2 += kmer_length
        else :
            overlap = kmer_length-1
            while overlap != 0 :
                
                kmer_scaffold1 = next_kmer_overlap (position1,overlap,kmer_dictionnary_scaffold1) 
                kmer_scaffold2= next_kmer_overlap (position2,overlap,kmer_dictionnary_scaffold2) 
                if kmer_scaffold1 == kmer_scaffold2 :
                    if kmer_scaffold1 in kmer_of_scaffolds :
                        kmer_of_scaffolds[kmer_scaffold1] += [[position1,position2]]
                        position1 += overlap
                        position2 += overlap
                    else :
                        kmer_of_scaffolds[kmer_scaffold1] = []
                        kmer_of_scaffolds[kmer_scaffold1] += [[position1,position2]]
                        position1 += overlap
                        position2 += overlap
                else :
                    overlap -=1
            
        print("position : "+ str(position1)+"/" +str(len(scaffold1)))
        print("position : "+ str(position2)+"/" +str(len(scaffold2)))


    #On arrive a la fin de l'un des scaffold au moin, on regarde si il est possible de trouver d'autre kmer chevauchant
    overlap = kmer_length -1
    while position1 + overlap >  len(scaffold1) or position2+overlap > len(scaffold2) : #Tant que le overlap depasse une des deux sequences on le diminue
        overlap -= 1
    
    kmer_scaffold1 = next_kmer_overlap (position1,kmer_length,kmer_dictionnary_scaffold1) 
    kmer_scaffold2= next_kmer_overlap (position2,kmer_length,kmer_dictionnary_scaffold2) 
    
    if kmer_scaffold1 == kmer_scaffold2 : 
        if kmer_scaffold1 in kmer_of_scaffolds :
            kmer_of_scaffolds[kmer_scaffold1] += [[position1,position2]]
        else :
            kmer_of_scaffolds[kmer_scaffold1] = []
            kmer_of_scaffolds[kmer_scaffold1] += [[position1,position2]]

    
    return kmer_of_scaffolds


print(extend_ancor(ancre[0],position_assemblage,position_reference,kmer_ass,kmer_ref,scaffold_ancre_assemblage,scaffold_ancre_reference))
