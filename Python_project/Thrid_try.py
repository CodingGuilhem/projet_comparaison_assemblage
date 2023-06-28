import subprocess
import sys
import os

def into_dictionary ( dictionary : dict, key,input_arg) -> dict : 
    """
    Function to put the input_arg in the dictionary, if the key already exist the input arg will be combined with the other arguments, else the key will be created
    Input : dict : the dictionary were you need to put the arg, a key of the dictionary, the arguments you want to add in your dictionary
    Output : dict : The same dictionary but ordered by the values of the dictionary
    """
    if key in dictionary :
        dictionary[key] += [input_arg]
    else :
        dictionary[key] = []
        dictionary[key] += [input_arg]


def recupererKmers(nomfichier,k) :
    """
    Function that forms a dictionary with the sequences of k-mers and their different positions, also creates a fasta file with the k-mers from the assemblies
    Input: String of the filename used to extract the k-mers, int the required K to form the k-mers
    Output: Dictionary of all the generated k-mers and their relative positions
    Side effect: Creates a fasta file containing all the k-mers in the current directory
    Comment: The position currently provided by gkampi is the absolute position of the k-mer in all the sequences (does not take into account scaffold changes and considers everything as a single sequence, so position 0 is for scaffold 1 only, and no other scaffolds. Additionally, it does not consider 'N' characters in its k-mer and position calculations). To add the positions relative to the scaffold, additional modifications are required.
    Comment: If the output files already exist, they will be overwritten.
    """
    nomSortie = nomfichier.split(".")[0] + "_kmers.fasta"
    result = subprocess.run(['gkampi','-P','-k',str(k),'-o',nomSortie,'-f','-Q',nomfichier],capture_output=True, text=True,universal_newlines=True)
    if result.returncode != 0 :
        print("Erreur dans la récuperation du fichier d'entrée")
        return
    #else :
        #print("K-mer stored in : "+str(nomSortie))

    kmer= {}

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
            

   
    #print("Dictionary done successfully") 
    return kmer




kmer_ref = recupererKmers("Test_sequences/rice_ref.fasta",100)
kmer_ass = recupererKmers("Test_sequences/rice_ass.fasta",100)

def extract_scaffold(file_name : str) -> dict:
    """
    Function that extract the scaffold of a fasta file 
    Input : string : name of the fasta 
    Output : dict: dictionary of the scaffold present in the fasta file and their absolute position ( form = {"scaffold" : [position of the first nucleotide of the scaffold]})
    Output : dict: dictionary of the scaffold present in the fasta file and their absolute position ( form = {"scaffold" : [position of the first nucleotide of the scaffold]})

    """
    scaffold_dictionary = {}
    scaffold_courant = ""
    absolute_position = 0
    nbr_scaffold = 0
    depart_scaffold = 0
    with open(file_name) as fasta_in:
        for line in fasta_in:
            if line[0] == ">":  
                if scaffold_courant != "":
                    into_dictionary(scaffold_dictionary,scaffold_courant,depart_scaffold)
                    scaffold_courant = ""
                    nbr_scaffold += 1
                    depart_scaffold = absolute_position
            else:
                scaffold_courant += line.strip()
                absolute_position += len(line)

          
                
  
        if scaffold_courant != "" and nbr_scaffold != 0:
            into_dictionary(scaffold_dictionary,scaffold_courant,depart_scaffold)
        else :
            into_dictionary(scaffold_dictionary,scaffold_courant,0)

    return scaffold_dictionary


scaffold_assemblage = extract_scaffold("Test_sequences/rice_ass.fasta")
scaffold_reference = extract_scaffold("Test_sequences/rice_ref.fasta")

def unique (dictionnaire_kmer,kmer) :
    """
    Function that test if a sequence is unique in the dictionary ( form = {"sequence" : [position(s)], ...})
    Input : dict : dictionary of all the sequences to compare  ( form = {"sequence" : [position(s)], ...}) , string : sequence to compare the dictionary with 
    Output : bool : True  = the sequence is alone in the dictionary 
    """
    return len(dictionnaire_kmer[kmer]) == 1


def trouver_ancre (dico_kmer1 : dict,dico_kmer2 : dict) -> tuple[str,int,int] :
    """
    Function that finds the first unique sequence present in both dictionary ( form = {"sequence" : [position(s)],...})
    Input: dict : first dictionary to compare , dict : second dictionary to compare, int : k-mer length
    Output: If a unique sequence anchor is found: returns the k-mer, position of the first nucleotide in the first sequence, position of the first nucleotide in the second sequence
    Output: Otherwise, returns None, None, None
    """

    for kmer1 in dico_kmer1.keys() :
        for kmer2 in dico_kmer2.keys() :
            if kmer1 == kmer2 and unique(dico_kmer1,kmer1) and unique(dico_kmer2,kmer2) :
                return kmer1,dico_kmer1[kmer1],dico_kmer2[kmer2]
    
    print("No anchor found",file = sys.stderr)
    return None, None, None


ancre = trouver_ancre(recupererKmers("Test_sequences/rice_ref.fasta",100),recupererKmers("Test_sequences/rice_ass.fasta",100))

def extraire_entete(filename : str,position : int) -> str :
    """
    Function that extracts the header of the scaffold containing the requested position
    Input: string: Name of the fasta file, int: Position of the nucleotide in the assembly
    Output: string: Header of the scaffold containing the requested position

    """
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

head = extraire_entete("Test_sequences/rice_ref.fasta",int(ancre[1][0]))

def trouve_scaffold(scaffold_dict : dict, position : int) -> str:
    """
    Function that find the scaffold containing the given position ()
    git Input : list : List of scaffold , int : Position searched 
    Output : Return the scaffold or None if the position isn't in the scaffold's list
    """
    pos = 0
    
    for scaffold in scaffold_dict.keys():
        pos += len(scaffold)
        if pos > position:
            return scaffold
    return None  


scaffold_ancre_assemblage = trouve_scaffold(scaffold_assemblage,int(ancre[2][0])) 
scaffold_ancre_reference = trouve_scaffold(scaffold_reference,int(ancre[1][0]))

def next_kmer_overlap (position : int ,overlap_value : int, kmer_dictionary : dict ) -> tuple[str,int]:
    """
    Function to find the next overlapping K-mer of the sequence using his absolute position
    Input : int : The absolute position of the K-mer, int : The number of nucleotide that overlap the K-mer, dict : A K-mer dictionary of the form : {"sequence" : [Absolute position], ...}
    Output : str : The K-mer overlapping the inputted K-mer, int : The absolute position of the K-mer overlapping
    """

    for kmer in kmer_dictionary.keys() : 
        if position + overlap_value == kmer_dictionary[kmer][0] :
            return  kmer,kmer_dictionary[kmer][0]

    #print("Next K-mer not found", file = sys.stderr)
    return None,None


def next_kmer_no_overlap (position : int, kmer_length : int , kmer_dictionary : dict) -> str :
    """
    Function to find the next K-mer without overlap of the sequence using his absolute position
    Input : int : The absolute position of the K-mer, int : The number of nucleotide that overlap the K-mer, dict : A K-mer dictionary of the form : {"sequence" : [Absolute position], ...}
    Output : str : The next K-mer of the inputted K-mer, int : The absolute position of the K-mer find

    """

    for kmer in kmer_dictionary.keys() :
        if position + kmer_length == kmer_dictionary[kmer][0] : 
            return kmer
    
    print("Next K-mer not found" , file=sys.stderr)
    return None

def previous_kmer_overlap(position :  int ,kmer_length : int, overlap_value : int , kmer_dictionary : dict) -> tuple[str,int] :
    """
    Function to find the previous overlapping K-mer of the sequence using his absolute position
    Input : int : The absolute position of the K-mer, int : the length of the K-mer generated, int : The number of nucleotide that overlap the K-mer, dict : A K-mer dictionary of the form : {"sequence" : [Absolute position], ...}
    Output : str : The previous K-mer overlapping the inputted K-mer, int : The absolute position of the K-mer find

    """


    for kmer in kmer_dictionary.keys() : 
        if position -kmer_length + overlap_value == kmer_dictionary[kmer][0] :
            return  kmer,kmer_dictionary[kmer][0]

    #print("Next K-mer not found", file = sys.stderr)
    return None,None

def previous_kmer_no_overlap (position : int, kmer_length : int , kmer_dictionary : dict) -> str :
    """
    Function to find the previous K-mer without overlap of the sequence using his absolute position
    Input : int : The absolute position of the K-mer, int : The number of nucleotide that overlap the K-mer, dict : A K-mer dictionary of the form : {"sequence" : [Absolute position], ...}
    Output : str : The K-mer overlapping the inputted K-mer
    """

    for kmer in kmer_dictionary.keys() :
        if position - kmer_length == kmer_dictionary[kmer][0] :
            return kmer

def display_charging_bar(percentage):
    """
    Display a charging bar in the 
    Input :  Percentage of the work done
    Output : Display the bar
    """
    bar_length = 50  # Length of the charging bar
    filled_length = int(bar_length * percentage / 100)
    bar = '#' * filled_length + '-' * (bar_length - filled_length)
    percentage_display = round(percentage, 1)  
    # Printing the loading bar
    print(f"Finding all the anchors : [{bar}] {percentage_display}%", end='\r')


def find_all_anchor(dico_kmer1 : dict,dico_kmer2 : dict) -> dict:
    """
    Function who find all the unique sequences present in the two dictionary 
    Input : dict : First dictionary of K-mers, dict : Second dictionary of K-mers (for of the dict : {"sequence" : [position(s) of the sequence]})
    Output : dict :  All the K-mers present only one time in each dictionary (form : {"sequence" : [position in ditc1 , position in dict2]})

    """

    anchor = {}
    pos = 0
    for kmer1 in dico_kmer1 : 
        pos +=1
        for kmer2 in dico_kmer2 :
            pos += 1
            if kmer1 == kmer2 and unique(dico_kmer1,kmer1) and unique(dico_kmer2,kmer2) :
                
                if kmer1 in anchor :
                    anchor[kmer1] = [dico_kmer1[kmer1][0],dico_kmer2[kmer2][0]]
                else : 
                    anchor[kmer1]= [] 
                    anchor[kmer1] = [dico_kmer1[kmer1][0],dico_kmer2[kmer2][0]]
        
        display_charging_bar((pos /(len(dico_kmer1)*len(dico_kmer2)))*100)
        #print("Finding all the anchors : " + str( (pos /(len(dico_kmer1)*len(dico_kmer2)))*100 ) + " %")
    
    return anchor
    

# Winning time by stocking the anchors of the sequences :
if not(os.path.exists("sequences_ancre.txt")) :
    all_anchors = find_all_anchor(kmer_ass, kmer_ref)
    with open("sequences_ancre.txt", "w") as file:
        for kmer, positions in all_anchors.items():
            line = f"{kmer}: {positions}\n"
            file.write(line)
else :
    all_anchors = {}
    with open("sequences_ancre.txt", "r") as file:
        for line in file:
            kmer, positions = line.strip().split(": ")
            positions = eval(positions)  
            all_anchors[kmer] = positions
# Ici all anchor = {sequence : [position_assemblage,position_reference]}
def delete_key(dictionnaire, clef):
    if clef in dictionnaire:
        del dictionnaire[clef]

def dictionary_in_dictionary(dictionary : dict, dictionary_to_add : dict) -> dict :

    """
    Function that add a dictionary into another dictionary
    Input : dict : The dictionary that will be added to the other dictionary, dict : The dictionary that will be added
    Output : dict : The first dictionary with the second dictionary added to it

    """
    for key in dictionary_to_add.keys() :
        into_dictionary(dictionary,key,dictionary_to_add[key][0])
    return dictionary

def extract_position_from_dictionary (dictionary : dict) -> tuple[int,int] :
    """
    Function that makes all the K-mer of the dictionary as one single string 
    Input :  dict : A dictionary of all the K-mer needed to be concatenated (form = {"sequence" : [position(s)]}) , str : The anchor K-mer : will be added into the dictionary, int : The anchor position
    Output : str : The concatenation of all the K-mer in the dictionary and the anchors in the order of their potions , int : the position of the first K-mer of the full anchor, int : the last position of the full anchor

    """
    position = list(dictionary.values())
    if len(position) == 1 :
        return min(position)[0], len(list(dictionary.keys())[0])
    else :
        return min(position)[0], max(position)[0]


def extend_anchor_front (start_position1 : int,start_position2 : int,kmer_length : int,kmer_dictionary1 : dict, kmer_dictionary2 : dict ,absolute_position_of_scaffold1 : int,absolute_position_of_scaffold2 : int, scaffold_length1 : int, scaffold_length2 : int) -> dict:
    """
    function that extend the back part of the anchor until the end of the scaffold
    Input : int : The position of the first K-mer of the anchor, dict : The dictionary of the K-mer of the scaffold, str : The scaffold, int : The absolute position of the scaffold
    Output : dict : The dictionary of the K-mer of the scaffold extended (form : {"sequence" : [position(s)]})
    """
    
    actual_position1 = start_position1
    actual_position2 = start_position2
    max_position1 = absolute_position_of_scaffold1 + scaffold_length1
    max_position2 = absolute_position_of_scaffold2 + scaffold_length2 
    kmer_of_anchors1 = {}
    kmer_of_anchors2 = {}

    while actual_position1 < max_position1 or actual_position2 < max_position2 :
        
        next_kmer1 = next_kmer_no_overlap(actual_position1,kmer_length,kmer_dictionary1)
        next_kmer2 = next_kmer_no_overlap(actual_position2,kmer_length,kmer_dictionary2)
        if actual_position1 + kmer_length > max_position1 or actual_position2 + kmer_length > max_position2 or next_kmer1 == None or next_kmer2 == None :
           
            overlap = kmer_length - 1
            next_kmer1 = next_kmer_overlap(actual_position1,overlap,kmer_dictionary1)
            next_kmer2 = next_kmer_overlap(actual_position2,overlap,kmer_dictionary2)
            
            while overlap > 0  and (next_kmer1 == (None,None) or next_kmer2 == (None,None)) :
                
                overlap -= 1
                next_kmer1 = next_kmer_overlap(actual_position1,overlap,kmer_dictionary1)
                next_kmer2 = next_kmer_overlap(actual_position2,overlap,kmer_dictionary2)

            actual_position1 -= overlap
            actual_position2 -= overlap

            if next_kmer1[0] == next_kmer2[0]  and  (next_kmer1 != (None,None) or next_kmer2 != (None,None)) :
                
                into_dictionary(kmer_of_anchors1, next_kmer1[0], actual_position1)
                into_dictionary(kmer_of_anchors2, next_kmer2[0], actual_position2)
                break
            else : 
                break
        else :
            actual_position1 += kmer_length
            actual_position2 += kmer_length
            into_dictionary(kmer_of_anchors1, next_kmer1, actual_position1)
            into_dictionary(kmer_of_anchors2, next_kmer2, actual_position2)

    return kmer_of_anchors1, kmer_of_anchors2
            
#print(extend_anchor_front(ancre[1][0],ancre[2][0],100,kmer_ass,kmer_ref,scaffold_assemblage[scaffold_ancre_assemblage][0],scaffold_reference[scaffold_ancre_reference][0],len(scaffold_ancre_assemblage),len(scaffold_ancre_reference)))
def extend_anchor_back (start_position1 : int,start_position2 : int,kmer_length : int,kmer_dictionary1 : dict, kmer_dictionary2 : dict ,absolute_position_of_scaffold1 : int,absolute_position_of_scaffold2 : int) -> tuple[dict,dict]:
    """
    function that extend the back part of the anchor until the end of the scaffold
    Input : int : The position of the first K-mer of the anchor, dict : The dictionary of the K-mer of the scaffold, str : The scaffold, int : The absolute position of the scaffold
    Output : dict : The dictionary of the K-mer of the scaffold extended (form : {"sequence" : [position(s)]})
    """
    
    actual_position1 = start_position1
    actual_position2 = start_position2
    kmer_of_anchors1 = {}
    kmer_of_anchors2 = {}

    while actual_position1 > absolute_position_of_scaffold1 or actual_position2 > absolute_position_of_scaffold2 : 
        previous_kmer1 = previous_kmer_no_overlap(actual_position1,kmer_length,kmer_dictionary1)
        previous_kmer2 = previous_kmer_no_overlap(actual_position2,kmer_length,kmer_dictionary2)

        if actual_position1 - kmer_length < 0 or previous_kmer1 == None or previous_kmer2 == None :
            overlap = kmer_length - 1
            previous_kmer1 = previous_kmer_overlap(actual_position1,kmer_length,overlap,kmer_dictionary1)
            previous_kmer2 = previous_kmer_overlap(actual_position2,kmer_length,overlap,kmer_dictionary2)
            while overlap > 0  and (previous_kmer1 == (None,None) or previous_kmer2 == (None,None)) :
                overlap -= 1
                previous_kmer1 = previous_kmer_overlap(actual_position1,overlap,kmer_dictionary1)
                previous_kmer2 = previous_kmer_overlap(actual_position2,overlap,kmer_dictionary2)
            actual_position1 -= overlap
            actual_position2 -= overlap

            if previous_kmer1[0] == previous_kmer2[0]  and  (previous_kmer1 != (None,None) or previous_kmer2 != (None,None)) :
                
                into_dictionary(kmer_of_anchors1, previous_kmer1[0], actual_position1)
                into_dictionary(kmer_of_anchors2, previous_kmer2[0], actual_position2)
                break
            else : 
                break

        else :
            actual_position1 -= kmer_length
            actual_position2 -= kmer_length
            into_dictionary(kmer_of_anchors1, previous_kmer1, actual_position1)
            into_dictionary(kmer_of_anchors2, previous_kmer2, actual_position2)

    return kmer_of_anchors1, kmer_of_anchors2

#print(extend_anchor_back(ancre[1][0],ancre[2][0],100,kmer_ass,kmer_ref,scaffold_assemblage[scaffold_ancre_assemblage][0],scaffold_reference[scaffold_ancre_reference][0]))

def ordered_dict_by_position (dictionary : dict) -> dict :
    """
    Function who return the dictionary ordered by values (smaller to higher) (form = {"sequence" : [value(s)],...})
    Input : dict : Non ordered dictionary
    Output : dict : The same dictionary but ordered by values 
    """
    return dict(sorted(dictionary.items(), key=lambda x: x[1]))

def find_next_anchor_from_position(actual_position,anchor_dictionary) :
    """
    Function that find the next anchor from the actual position you are 
    Input : int : actual position of your program, dict : A dictionary that contains all the anchors of the scaffold with their position (form = {"anchor" : [position(s)]})
    Output : int : The position of the next kmer
    """
    ordered_dict = ordered_dict_by_position(anchor_dictionary)
    for position in anchor_dictionary.values() :
        print(position)
        if position > actual_position :
            return position 

def multi_anchor_extension (anchor_dictionary : dict, absolute_position_scaffold1 : int, absolute_position_scaffold2 : int,scaffold_length1 :  int, scaffold_length2 : int, kmer_length: int, kmer_dictionary1 : dict,kmer_dictionary2 : dict ) -> tuple[dict,dict] :
    max_position1 = absolute_position_scaffold1 + scaffold_length1
    max_position2 = absolute_position_scaffold2 + scaffold_length2
    
    anchor_of_scaffold = {}
    for anchor in anchor_dictionary.keys() :
        if anchor_dictionary[anchor][0] > absolute_position_scaffold1 and anchor_dictionary[anchor][1] and (anchor_dictionary[anchor][0] + kmer_length) <= max_position1 and (anchor_dictionary[anchor][1] + kmer_length) <= max_position2  :
            into_dictionary(anchor_of_scaffold,anchor,anchor_dictionary[anchor][0])
            into_dictionary(anchor_of_scaffold,anchor,anchor_dictionary[anchor][1])

    front_anchor1 = {}
    front_anchor2 = {}
    actual_position1 = absolute_position_scaffold1
    actual_position2 = absolute_position_scaffold2
   

    while actual_position1 < max_position1 and actual_position2 < max_position2 :
        anchor = extend_anchor_front(actual_position1,actual_position2,kmer_length,kmer_dictionary1,kmer_dictionary2,absolute_position_scaffold1,absolute_position_scaffold2,scaffold_length1,scaffold_length2)
       
        dictionary_in_dictionary(front_anchor1,anchor[0])
        dictionary_in_dictionary(front_anchor2,anchor[1])
        
        positions1 = extract_position_from_dictionary(anchor[0])
        positions2 = extract_position_from_dictionary(anchor[1])
        print(positions1[1])
        position_next_anchor1 = find_next_anchor_from_position(positions1[1],anchor_of_scaffold)
        position_next_anchor2 = find_next_anchor_from_position(positions2[1],anchor_of_scaffold)
        print(position_next_anchor1)
        print(position_next_anchor2)
        
        actual_position1 = position_next_anchor1[1]
        actual_position2 = position_next_anchor2[1]
    # print(front_anchor1)
    # print(front_anchor2)


        
    back_anchor = {}


multi_anchor_extension(all_anchors,scaffold_assemblage[scaffold_ancre_assemblage][0],scaffold_reference[scaffold_ancre_reference][0],len(scaffold_ancre_assemblage),len(scaffold_ancre_reference),100,kmer_ass,kmer_ref)


    