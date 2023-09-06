import subprocess

#A faire : 
# trouver une ancre unique sur les deux alignement ( sequence + position + n° scaffold) 
# Extraire les identifiants des sequences ( fonction ou gkampi ou regeps)
# Trouver les kmers qui suivent :pas besoin qu'ils soit unique seulement l'ancre doit être unique
# Arranger la fonction proochain kmer pour ces besoins  sikhr


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

# def chevauchement (sequence1,sequence2) :
#     chevauchement = -1
#     for itterateur in range(len(sequence1)) :
#         if sequence1[len(sequence1)-itterateur:] == sequence2[:itterateur] and itterateur != 0:
#             chevauchement = itterateur

#     if chevauchement == -1 :
#         return(False,None)          
#     else :
#         return (True,chevauchement)
def chevauchement(sequence1, sequence2):
    chevauchement = 0
    for itterateur in range(len(sequence1)):
        
        if sequence1[:len(sequence1)-itterateur] == sequence2[itterateur:]:
            chevauchement = len(sequence1)-itterateur
        
    return chevauchement

# Function that is testing if one of the suffix of the sequence 1 is corresponding to one of the prefix of the sequence 2
# Input : string : first sequence , string : second sequence 
# Output :  int : Lenght of the overlapping sequence 

def chevauchement_inverse (sequence1,sequence2) :
    chevauchement = -1
    for iterateur in range(len(sequence1)) :
        if sequence1[:iterateur] == sequence2[len(sequence2)-iterateur:] and iterateur != 0 :
            chevauchement = iterateur

    if chevauchement == -1 :
        return(False,None)  
    else :
        return True, chevauchement

# Function that is testing if one of the prefix of the sequence 1 is corresponding to one of the suffix of the sequence 2
# Input : string : first sequence , string : second sequence 
# Output :  int : Lenght of the overlapping sequence 

def kmers_chevauchant(dictionnaire_sequence,kmer) :
    sequence_chevauchement = {}
    for sequences in dictionnaire_sequence.keys() :
        chev = chevauchement(sequences,kmer)
        if chev != 0  and sequences in sequence_chevauchement:
            sequence_chevauchement[sequences] += [chev]    
        elif chev != 0 and sequences not in sequence_chevauchement:
            sequence_chevauchement[sequences]= []
            sequence_chevauchement[sequences] = [chev]
    return sequence_chevauchement
# Function making a list of all the overlaping sequences between the sequence the dictionnary of sequences (suffix of kmer =  prefix of dictionnary sequence)
# Input :  dict : dictionnary of the sequences (form = {"sequence" : [position(s)] ,...,"sequence" : [position(s)]}) , string : sequence to compare the dictionnary with
# Output : dict : dictionnary of the sequences that overlap the sequences previously given ( form = {"sequence" : [length of the overlapping sequence]})

def kmers_chevauchant_inverse(dictionnaire_sequence,kmer) :
    sequence_chevauchement = {}
    for sequences in dictionnaire_sequence.keys() :
        chev = chevauchement_inverse(sequences,kmer)
        if chev[0] and sequences not in sequence_chevauchement:
            sequence_chevauchement[sequences]= []
            sequence_chevauchement[sequences] = [chev[1]]
        elif chev[0] and sequences in sequence_chevauchement:
            sequence_chevauchement[sequences] += [chev[1]]    
   
    return sequence_chevauchement
# Function making a list of all the overlaping sequences between the sequence the dictionnary of sequences (prefix of kmer =  suffix of dictionnary sequence)
# Input :  dict : dictionnary of the sequences (form = {"sequence" : [position(s)] ,...,"sequence" : [position(s)]}) , string : sequence to compare the dictionnary with
# Output : dict : dictionnary of the sequences that overlap the sequences previously given ( form = {"sequence" : [length of the overlapping sequence]})

def extract_scaffold(file_name):
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



def trouver_ancre (dico_kmer1,dico_kmer2) :

    for kmer1 in dico_kmer1.keys() :
        for kmer2 in dico_kmer2.keys() :
            if kmer1 == kmer2 and unique(dico_kmer1,kmer1) and unique(dico_kmer2,kmer2) :
                return kmer1,dico_kmer1[kmer1],dico_kmer2[kmer2]
    
    print("Aucune ancre n'a ete trouvee ")
    return None, None, None
# Function that finds the first unique sequence present in both dictionnary ( form = {"sequence" : [position(s)],...})
# Input: dict : first dictionnary to compare , dict : second dictionnary to compare, int : k-mer length
# Output: If a unique sequence anchor is found: returns the k-mer, position of the first nucleotide in the first sequence, position of the first nucleotide in the second sequence
# Output: Otherwise, returns None, None, None

ancre = trouver_ancre(kmer_ref,kmer_ass)
#a revoir : marche pas si dico dans l'autre sens???

def extraire_entete(filename,position) :
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


def trouve_scaffold(liste_scaffold, position):
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




def position_in_scaffold (scaffold,kmer) :
    itterateur =0
    
    while itterateur + len(kmer) < len(scaffold) :
        if scaffold[itterateur : itterateur+len(kmer)] == kmer :
            return itterateur
        itterateur+=1
    print("Position non trouvee")
    return None
# Function that returns the position of the kmer in the given scaffold if it is present 
# Input: string : scaffold that contains the kmer , string  : kmer searched 
# Output: int position of the start of the kmer in the scaffold, or None

position_assemblage = position_in_scaffold(scaffold_ancre_assemblage,ancre[0])
position_reference = position_in_scaffold(scaffold_ancre_reference,ancre[0])



def next_kmer_non_chevauchant(scaffold,kmer,position) :
    longueur_kmer = len(kmer)
    if len(scaffold) > position+2*longueur_kmer :
        prochain_kmer = scaffold[position: position+longueur_kmer]
        return prochain_kmer
    else :
        print("Scaffold finis")
        return None
# Function that find the next kmer non overlapping kmer in the scaffold
# Input : string : Scaffold containing the kmer , string : kmer in the scaffold, int : position of the first nucleotide of the kmer in the scaffold  
# Output :  string : Next K-mer of the scaffold ( no overlap )


def precedent_kmer_non_chevauchant(scaffold,kmer,position) :
    longueur_kmer = len(kmer)
    if position-longueur_kmer > 0 :
        prochain_kmer = scaffold[position-longueur_kmer : position]
        return prochain_kmer
    else :
        print("Scaffold trop court")
    return None
# Function that find the previous K-mer non overlapping kmer in the scaffold
# Input : string : Scaffold containing the kmer , string : kmer in the scaffold, int : position of the first nucleotide of the kmer in the scaffold  
# Output :  string : Previous K-mer of the scaffold ( no overlap )


def next_kmer_chevauchant(scaffold,kmer,position,dico_kmer) :
    chevauchement = kmers_chevauchant(dico_kmer,kmer)
    for sequence_chevauchante in chevauchement.keys() :
        if kmer == scaffold[3544:3644] :
            if sequence_chevauchante == scaffold[position+(len(kmer)-chevauchement[sequence_chevauchante][0]): position +((2*len(kmer))-chevauchement[sequence_chevauchante][0])]:
                print(position + 2*(len(kmer)-chevauchement[sequence_chevauchante][0]))
                print(chevauchement[sequence_chevauchante][0])
                print(scaffold[3544:3644] in chevauchement)
        if len(scaffold)-1 > position + 2*(len(kmer)-chevauchement[sequence_chevauchante][0]) and sequence_chevauchante == scaffold[position+(len(kmer)-chevauchement[sequence_chevauchante][0]): position +((2*len(kmer))-chevauchement[sequence_chevauchante][0])] and chevauchement[sequence_chevauchante][0] != len(kmer):
            
            return sequence_chevauchante,chevauchement[sequence_chevauchante][0]
        
    #print("Pas de prochain kmer, scaffold trop court")
    return None,None
# Function that find the next overlapping K-mer in the scaffold
# Input : string : Scaffold containing the K-mer , string : K-mer in the scaffold, int : position of the first nucleotide of the K-mer in the scaffold, dict : dictionnary of all the K-mer (form = {"sequence" : [position(s)]})  
# Output :  string : Next K-mer of the scaffold 

def precedent_kmer_chevauchant (scaffold,kmer,position,dico_kmer) :
    chevauchement = kmers_chevauchant_inverse(dico_kmer,kmer)
    keys = chevauchement.keys()
    for sequence_chevauchante in keys :
        if position - (len(kmer)-chevauchement[sequence_chevauchante][0])+len(kmer) > 0 :
            if sequence_chevauchante == scaffold[position -(len(kmer)+chevauchement[sequence_chevauchante][0]): position] :
                return sequence_chevauchante
        else :
            print("Pas de precedent kmer, scaffold trop court")
            return None
        
# Function that find the previous overlapping K-mer in the scaffold
# Input : string : Scaffold containing the K-mer , string : K-mer in the scaffold, int : position of the first nucleotide of the K-mer in the scaffold, dict : dictionnary of all the K-mer (form = {"sequence" : [position(s)]})  
# Output :  string : Previous K-mer of the scaffold 


def prochain_kmer_in_scaffold_non_chevauchant (scaffold1,kmer,position,dico_kmer) :

    
    apres_ancre = {}
    non_correspondant = {}
    pos = position + len(kmer)
    while len(scaffold1) > pos :
        kmer_actuelle = scaffold1[pos : pos + len(kmer)]
        prochain = next_kmer_non_chevauchant(scaffold1,kmer_actuelle,pos)
        if prochain == scaffold1[pos : pos+len(kmer)] :
            if prochain in dico_kmer :
        
                if prochain in apres_ancre :
                    apres_ancre[prochain] += [pos]
                else :
                    apres_ancre[prochain] = []
                    apres_ancre[prochain] += [pos]
        else :
            if prochain != None :
                if prochain in non_correspondant :
                    non_correspondant[prochain] += [pos]
                else :
                    non_correspondant[prochain] = []
                    non_correspondant[prochain] += [pos]

        pos+= len(kmer)
    
    return apres_ancre,non_correspondant

# Function that find  all the next K-mer in the scaffold present in the dictionnary (without overlap)
# Input : string : The scaffold , string : The K-mer you already now in the sequence, String : the position of the known K-mer, dict : The dictionnary of all the kmer you have (form = {"sequence" : [position(s)], ...})
# Output : dict : All the K-mer present in the scaffold and the dictionnary after the known K-mer (without overlap)

def precedent_kmer_in_scaffold_non_chevauchant(scaffold,kmer,position,dico_kmer) : 
    avant_ancre = {}
    non_correspondant = {}
    pos = position 

    while pos > 0 :
        precedent = precedent_kmer_non_chevauchant(scaffold,kmer,pos)
        if precedent == scaffold[pos-len(kmer) : pos] :
            if precedent in dico_kmer :

                if precedent in avant_ancre :
                    avant_ancre[precedent] += [pos]
                else :
                    avant_ancre[precedent] = []
                    avant_ancre[precedent] += [pos]
        else :
            if precedent != None :
                if precedent in non_correspondant :
                    non_correspondant[precedent] += [pos]
                else :
                    non_correspondant[precedent] = []
                    non_correspondant[precedent] += [pos]

        pos-= len(kmer)

    return avant_ancre ,non_correspondant
# Function that find  all the previous K-mer in the scaffold present in the dictionnary (without overlap)
# Input : string : The scaffold , string : The K-mer you already now in the sequence, String : the position of the known K-mer, dict : The dictionnary of all the kmer you have (form = {"sequence" : [position(s)], ...})
# Output : dict : All the K-mer present in the scaffold and the dictionnary before the known K-mer (without overlap)

#kmer = precedent_kmer_in_scaffold_non_chevauchant(scaffold_ancre_reference,ancre[0],position_reference,kmer_ref)
 #print(kmer)


def prochain_kmer_in_scaffold_chevauchant (scaffold, kmer, position,dico_kmer) :
    apres_ancre = {}
    non_correspondant = {}
    pos = position 
    
    while pos < len (scaffold)-1 :
        kmer_actuelle = scaffold[pos : pos+len(kmer)]
        if pos +len(kmer) < len(scaffold) :
            prochain = next_kmer_chevauchant(scaffold,kmer_actuelle,pos,dico_kmer)
            
            chevauchement = prochain [1]
        
        
            if prochain[0] == scaffold[pos+len(kmer)-chevauchement : pos+2*len(kmer)-chevauchement] :
                if prochain[0] in dico_kmer :
                    if prochain[0] in apres_ancre :
                        apres_ancre[prochain[0]] += [pos]
                    else :
                        
                        apres_ancre[prochain[0]] = []
                        apres_ancre[prochain[0]] += [pos]
                else :
                    if prochain[0] != None :
                        if prochain[0] in non_correspondant :
                            non_correspondant[prochain[0]] += [pos]
                    else :
                        non_correspondant[prochain[0]] = []
                        non_correspondant[prochain[0]] += [pos]
                
        pos+= len(kmer)-chevauchement
        
        print("progression en cours :" + str((pos/len(scaffold))*100)+ " % du scaffold effectue" )
        print("Position = " + str(pos) + "/" + str(len(scaffold)))
        
    return apres_ancre ,non_correspondant
# Function that find  all the next overlaping K-mer in the scaffold present in the dictionnary 
# Input : string : The scaffold , string : The K-mer you already now in the sequence, String : the position of the known K-mer, dict : The dictionnary of all the kmer you have (form = {"sequence" : [position(s)], ...})
# Output : dict : All the overlaping K-mer present in the scaffold and the dictionnary after the known K-mer

print(prochain_kmer_in_scaffold_chevauchant(scaffold_ancre_assemblage,ancre[0],position_assemblage,kmer_ass))

def precedent_kmer_in_scaffold_chevauchant(scaffold,kmer,position,dico_kmer) : 
    avant_ancre = {}
    non_correspondant = {}
    pos = position 

    while pos > 0 :
        precedent = precedent_kmer_chevauchant(scaffold,kmer,pos,dico_kmer)
        if precedent == scaffold[pos-len(kmer) : pos] :
            if precedent in dico_kmer :

                if precedent in avant_ancre :
                    avant_ancre[precedent] += [pos]
                else :
                    avant_ancre[precedent] = []
                    avant_ancre[precedent] += [pos]
        else :
            if precedent != None :
                if precedent in non_correspondant :
                    non_correspondant[precedent] += [pos]
                else :
                    non_correspondant[precedent] = []
                    non_correspondant[precedent] += [pos]

        pos-= len(kmer)

    return avant_ancre ,non_correspondant
# Function that find  all the previous K-mer from the ancor in the scaffold present in the dictionnary (without overlap)
# Input : string : The scaffold , string : The K-mer you already now in the sequence, String : the position of the known K-mer, dict : The dictionnary of all the kmer you have (form = {"sequence" : [position(s)], ...})
# Output : dict : All the K-mer present in the scaffold and the dictionnary before the known K-mer (without overlap)

#kmer = precedent_kmer_in_scaffold_non_chevauchant(scaffold_ancre_reference,ancre[0],position_reference,kmer_ref)
#print(kmer)

def scaffold_identique_non_chevauchant (scaffold1,scaffold2,ancre,dico_kmer1,dico_kmer2,position_ancre1,position_ancre2) :
    total_identique  = 0
    identiques = {}
    precedent1 = precedent_kmer_in_scaffold_non_chevauchant(scaffold1,ancre[0],position_ancre1,dico_kmer1)[0]
    precedent2 = precedent_kmer_in_scaffold_non_chevauchant(scaffold2,ancre[0],position_ancre2,dico_kmer2)[0]

    for kmer1 in precedent1.keys() :
        for kmer2 in precedent2.keys() :
            if kmer1 == kmer2 :
                total_identique += 1
                if kmer1 in identiques :
                    identiques[kmer1] += [precedent1[kmer1][0],precedent2[kmer1][0]] #Position des kmers
                else :
                    identiques[kmer1] = []
                    identiques[kmer1] += [precedent1[kmer1][0],precedent2[kmer1][0]] 
    
    suivant1 = prochain_kmer_in_scaffold_non_chevauchant(scaffold1,ancre[0],position_ancre1,dico_kmer1)[0]
    suivant2 = prochain_kmer_in_scaffold_non_chevauchant(scaffold2,ancre[0],position_ancre2,dico_kmer2)[0]

    for kmer1 in suivant1.keys() :
        for kmer2 in suivant2.keys() :
            if kmer1 == kmer2 :
                total_identique += 1
                if kmer1 in identiques :
                    identiques[kmer1] += [suivant1[kmer1][0],suivant2[kmer1][0]] #Position des kmers
                else :
                    identiques[kmer1] = []
                    identiques[kmer1] += [suivant1[kmer1][0],suivant2[kmer1][0]] 
    
    total = len(suivant1)+len(suivant2)+len(precedent1)+len(precedent2)
    pourcentage_identique = total_identique / total

    return identiques,pourcentage_identique
# Function that say if the scaffolds contains the same K-mer from the ancor to the end of the scaffold ( no overlap )
# Input : string : first scaffold, string : second scaffold, string : ancor, dict : dictionnary of all the K-mer of the first scaffold (form = {"sequence" : [position(s),...]}),  dict : dictionnary of all the K-mer of the second scaffold (form = {"sequence" : [position(s),...]}), position of the ancor in the scaffold 1, position of the ancor in the scaffold 2
# Output : list : List of the K-mer that are the same in the two scaffold, Percentage of identicity

# identique = scaffold_identique_non_chevauchant(scaffold_ancre_reference,scaffold_ancre_assemblage,ancre,kmer_ref,kmer_ass,position_reference,position_assemblage)
# print(identique)

def scaffold_identique_chevauchant (scaffold1,scaffold2,ancre,dico_kmer1,dico_kmer2,position_ancre1,position_ancre2) :
    total_identique  = 0
    identiques = {}
    precedent1 = precedent_kmer_in_scaffold_chevauchant(scaffold1,ancre[0],position_ancre1,dico_kmer1)[0]
    precedent2 = precedent_kmer_in_scaffold_chevauchant(scaffold2,ancre[0],position_ancre2,dico_kmer2)[0]

    for kmer1 in precedent1.keys() :
        for kmer2 in precedent2.keys() :
            if kmer1 == kmer2 :
                total_identique += 1
                if kmer1 in identiques :
                    identiques[kmer1] += [precedent1[kmer1][0],precedent2[kmer1][0]] #Position des kmers
                else :
                    identiques[kmer1] = []
                    identiques[kmer1] += [precedent1[kmer1][0],precedent2[kmer1][0]] 
    
    suivant1 = prochain_kmer_in_scaffold_chevauchant(scaffold1,ancre[0],position_ancre1,dico_kmer1)[0]
    suivant2 = prochain_kmer_in_scaffold_chevauchant(scaffold2,ancre[0],position_ancre2,dico_kmer2)[0]

    for kmer1 in suivant1.keys() :
        for kmer2 in suivant2.keys() :
            if kmer1 == kmer2 :
                total_identique += 1
                if kmer1 in identiques :
                    identiques[kmer1] += [suivant1[kmer1][0],suivant2[kmer1][0]] #Position des kmers
                else :
                    identiques[kmer1] = []
                    identiques[kmer1] += [suivant1[kmer1][0],suivant2[kmer1][0]] 
    
    total = len(suivant1)+len(suivant2)+len(precedent1)+len(precedent2)
    print(total)
    print(total_identique)
    pourcentage_identique = total_identique / total

    return identiques,pourcentage_identique
# Function that say if the scaffolds contains the same K-mer from the ancor to the end of the scaffold 
# Input : string : first scaffold, string : second scaffold, string : ancor, dict : dictionnary of all the K-mer of the first scaffold (form = {"sequence" : [position(s),...]}),  dict : dictionnary of all the K-mer of the second scaffold (form = {"sequence" : [position(s),...]}), position of the ancor in the scaffold 1, position of the ancor in the scaffold 2
# Output : list : List of the K-mer that are the same in the two scaffold, Percentage of identicity

# identique = scaffold_identique_chevauchant(scaffold_ancre_reference,scaffold_ancre_assemblage,ancre,kmer_ref,kmer_ass,position_reference,position_assemblage)
# print(identique)