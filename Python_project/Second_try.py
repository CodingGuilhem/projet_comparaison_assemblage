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

def position_in_scaffold (scaffold : str,kmer : str) -> int :
    """
    Function that returns the position of the kmer in the given scaffold if it is present 
    Input: string : scaffold that contains the kmer , string  : kmer searched 
    Output: int position of the start of the kmer in the scaffold, or None
    """
    iterator =0
    
    while iterator + len(kmer) < len(scaffold) :
        if scaffold[iterator : iterator+len(kmer)] == kmer :
            return iterator
        iterator+=1
    print("Position not found",file = sys.stderr)
    return None


position_assemblage = position_in_scaffold(scaffold_ancre_assemblage,ancre[0])
position_reference = position_in_scaffold(scaffold_ancre_reference,ancre[0])

def next_kmer_overlap (position : int ,overlap_value : int, kmer_dictionary : dict ) -> tuple[str,int]:
    """
    Function to find the next overlapping K-mer of the sequence using his absolute position
    Input : int : The absolute position of the K-mer, int : The number of nucleotide that overlap the K-mer, dict : A K-mer dictionary of the form : {"sequence" : [Absolute position], ...}
    Output : str : The K-mer overlapping the inputted K-mer, int : The absolute position of the K-mer overlapping
    """

    for kmer in kmer_dictionary.keys() : 
        if position + overlap_value == kmer_dictionary[kmer][0] :
            return  kmer,kmer_dictionary[kmer][0]

    print("Next K-mer not found", file = sys.stderr)
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

    print("Next K-mer not found", file = sys.stderr)
    return None,None

def previous_kmer_no_overlap (position : int, kmer_length : int , kmer_dictionary : dict) -> str :
    """
    Function to find the previous K-mer without overlap of the sequence using his absolute position
    Input : int : The absolute position of the K-mer, int : The number of nucleotide that overlap the K-mer, dict : A K-mer dictionary of the form : {"sequence" : [Absolute position], ...}
    Output : str : The K-mer overlapping the inputted K-mer, int : The absolute position of the K-mer overlapping

    """
    for kmer in kmer_dictionary.keys() :
        if position - kmer_length == kmer_dictionary[kmer][0] :
            return kmer



def extend_anchor_front (kmer_length : int ,position1 : int , position2 : int , kmer_dict1: dict , kmer_dict2 : dict) :
    """
    Function that return the next kmer of the two scaffold if he is different return None
    Input  :  int : length of the K-mer, int : aABSOLUTE position of the first K-mer , int : ABSOLUTE position of the second K-mer, dict : dictionary of all the K-mer of the first scaffold, dict : dictionary of all the K-mers of the second scaffold ( form  = {"sequence" : [position(s)]})
    Output :  dict :  The next K-mer (form = {"sequence" : [Absolute position in dic1,  Absolute position in dic 2]})

    """

    kmer1 = next_kmer_no_overlap(position1,kmer_length,kmer_dict1)
    kmer2 = next_kmer_no_overlap(position2,kmer_length,kmer_dict2)
    if kmer1 == kmer2  :
        return kmer1,position1+kmer_length,position2+kmer_length
    else :
        overlap = kmer_length - 1
        while overlap != 0 : 
            kmer1 = next_kmer_overlap(position1,overlap,kmer_dict1)
            kmer2 = next_kmer_overlap(position2,overlap,kmer_dict2)
            if kmer1 == kmer2 :
                return kmer1,position1+kmer_length,position2+kmer_length
            overlap -= 1
            
        #print("No K-mer found" , file = sys.stderr)
        return None,None,None

#print(extend_anchor_front(100,ancre[2][0],ancre[1][0],kmer_ass,kmer_ref))

def extend_anchor_back (kmer_length : int ,position1 : int , position2 : int , kmer_dict1: dict , kmer_dict2 : dict) :
    """
    Function that return the previous kmer of the two scaffold if he is different return None
    Input  :  int : length of the K-mer, int : aABSOLUTE position of the first K-mer , int : ABSOLUTE position of the second K-mer, dict : dictionary of all the K-mer of the first scaffold, dict : dictionary of all the K-mers of the second scaffold ( form  = {"sequence" : [position(s)]})
    Output :  dict :  The previous K-mer (form = {"sequence" : [Absolute position in dic1,  Absolute position in dic 2]})

    """
    kmer1 = previous_kmer_no_overlap(position1,kmer_length,kmer_dict1)
    kmer2 = previous_kmer_no_overlap(position2,kmer_length,kmer_dict2)
    if kmer1 == kmer2  :
        return kmer1,position1-kmer_length,position2-kmer_length
    else :
        overlap = kmer_length - 1
        while overlap != 0 : 
            
            kmer1 = previous_kmer_overlap(position1,kmer_length,overlap,kmer_dict1)
            kmer2 = previous_kmer_overlap(position2,kmer_length,overlap,kmer_dict2)
            if kmer1 == kmer2 :
                return kmer1,position1-kmer_length,position2-kmer_length
            overlap -= 1
            
        #print("No K-mer found" , file = sys.stderr)
        return None,None,None


#print(extend_anchor_back(100,ancre[1][0],ancre[2][0],kmer_ref,kmer_ass))

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



def extend_anchor_full(anchor : str, scaffold1 : str, scaffold2 : str ,position1 : int, position2 : int, kmer_dict1 : dict, kmer_dict2 : dict) -> tuple[dict,dict]: 

    """
    Function to extend the anchor, return all the K-mers present in the two sequences linked directly with the anchor 
    Input : str : anchor , str : First scaffold that contains the anchor ,  str :  Second scaffold that contains the anchor, int : ABSOLUTE position of the anchor in the first scaffold , int : ABSOLUTE position of the anchor in the second scaffold , dict : dictionary of all the K-mers present in the second scaffold (dict form = {"sequence" : [position(s),...]})
    Output : Dictionary of all the scaffold of the extended anchor (form = {"K-mer" : [absolute position in scaffold1, absolute position in scaffold2]})

    """
    kmer_of_scaffold1 = {}
    kmer_of_scaffold2 = {}
    kmer_length = len(anchor)
    position_in_scaffold1 = position_in_scaffold(scaffold1, anchor)
    position_in_scaffold2 = position_in_scaffold(scaffold2, anchor)

    absolute_position1 = position1
    absolute_position2 = position2
    if position_in_scaffold1 == None or position_in_scaffold2 == None :
        return kmer_of_scaffold1,kmer_of_scaffold2
    
    while len(scaffold1) > position_in_scaffold1+kmer_length and len(scaffold2) > position_in_scaffold2 +kmer_length :
        next_kmer = extend_anchor_front(kmer_length,absolute_position1,absolute_position2,kmer_dict1,kmer_dict2)
        if next_kmer[0] != None :
            position_in_scaffold1 += next_kmer[1] - absolute_position1
            position_in_scaffold2 += next_kmer[2] - absolute_position2 
            
            absolute_position1 = next_kmer[1]
            absolute_position2 = next_kmer[2]
    
            into_dictionary(kmer_of_scaffold1,next_kmer[0],absolute_position1)
            into_dictionary(kmer_of_scaffold2,next_kmer[0],absolute_position2)
        else : #Getting out of the while loop because the anchor aren't the same for both scaffold 
            position_in_scaffold1 = len(scaffold1) 
            position_in_scaffold2 = len(scaffold2)

    position_in_scaffold1 = position_in_scaffold(scaffold1, anchor)
    position_in_scaffold2 = position_in_scaffold(scaffold2, anchor)
    
    absolute_position1 = position1
    absolute_position2 = position2
    
    while position_in_scaffold1 > 0 and position_in_scaffold2 > 0 :   

       
        previous_kmer = extend_anchor_back(kmer_length,absolute_position1,absolute_position2,kmer_dict1,kmer_dict2)

        if previous_kmer[0] != None :
            position_in_scaffold1  -= absolute_position1 -previous_kmer[1] 
            position_in_scaffold2 -= absolute_position2 - previous_kmer[2] 
            absolute_position1 = previous_kmer[1]

            absolute_position2 = previous_kmer[2]
            into_dictionary(kmer_of_scaffold1,previous_kmer[0],absolute_position1)
            into_dictionary(kmer_of_scaffold2,previous_kmer[0],absolute_position2)
        else :
            position_in_scaffold1 = 0
            position_in_scaffold2 = 0

    return kmer_of_scaffold1,kmer_of_scaffold2
 
anchors = extend_anchor_full(ancre[0],scaffold_ancre_assemblage,scaffold_ancre_reference,ancre[2][0],ancre[1][0],kmer_ass,kmer_ref)

def ordered_dict_by_position (dictionary : dict) -> dict :
    """
    Function who return the dictionary ordered by values (smaller to higher) (form = {"sequence" : [value(s)],...})
    Input : dict : Non ordered dictionary
    Output : dict : The same dictionary but ordered by values 

    """
    return dict(sorted(dictionary.items(), key=lambda x: x[1]))

def create_full_anchor(anchor_dictionary : dict, anchor : str , position : int ) -> tuple[str,int,int] :
    """
    Function that makes all the K-mer of the dictionary as one single string 
    Input :  dict : A dictionary of all the K-mer needed to be concatenated (form = {"sequence" : [position(s)]}) , str : The anchor K-mer : will be added into the dictionary, int : The anchor position
    Output : str : The concatenation of all the K-mer in the dictionary and the anchors in the order of their potions , int : the position of the first K-mer of the full anchor, int : the last position of the full anchor

    """
    into_dictionary(anchor_dictionary,anchor,position)
    ordered_dict = ordered_dict_by_position(anchor_dictionary)
    full_anchor = ""
    for kmer in ordered_dict.keys() :
        full_anchor += kmer
    
    return full_anchor, list(ordered_dict.values())[0][0],list(ordered_dict.values())[-1][0]
full_anchor1 = create_full_anchor(anchors[0],ancre[0],position_assemblage)

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
# Return the max and min value of a dictionary (form = {"sequence" : [value(s)],...})
def dictionary_in_dictionary(dictionary : dict, dictionary_to_add : dict) -> dict :

    """
    Function that add a dictionary into another dictionary
    Input : dict : The dictionary that will be added to the other dictionary, dict : The dictionary that will be added
    Output : dict : The first dictionary with the second dictionary added to it

    """
    for key in dictionary_to_add.keys() :
        into_dictionary(dictionary,key,dictionary_to_add[key][0])
    return dictionary

def multi_anchor_extension_one_scaffold(scaffold1 : str,scaffold2 : str ,scaffold_dictionary1 : dict, scaffold_dictionary2 : dict, anchor_dictionary : dict, kmer_dictionary1 : dict, kmer_dictionary2 : dict ) -> dict:
   
    anchor_of_scaffolds = {}
    for anchor in anchor_dictionary.keys() :
        if kmer_dictionary1[anchor][0]>= scaffold_dictionary1[scaffold1][0] and kmer_dictionary2[anchor][0] >= scaffold_dictionary2[scaffold2][0] and  kmer_dictionary1[anchor][0]< scaffold_dictionary1[scaffold1][0]+len(scaffold1) and kmer_dictionary2[anchor][0] < scaffold_dictionary2[scaffold2][0]+len(scaffold2) :
            into_dictionary(anchor_of_scaffolds,anchor, anchor_dictionary[anchor][0] )
    anchor_of_scaffolds = ordered_dict_by_position(anchor_of_scaffolds)

    kmer_number =0
    kmer_list = list(anchor_of_scaffolds.keys())
    
    full_anchor = extend_anchor_full(kmer_list[kmer_number],scaffold1,scaffold2,kmer_dictionary1[kmer_list[kmer_number]][0],kmer_dictionary2[kmer_list[kmer_number]][0],kmer_dictionary1,kmer_dictionary2)
    
    if full_anchor[0] == {} or full_anchor[1] == {} :
        kmer_number += 1
        full_anchor = extend_anchor_full(kmer_list[kmer_number],scaffold1,scaffold2,kmer_dictionary1[kmer_list[kmer_number]][0],kmer_dictionary2[kmer_list[kmer_number]][0],kmer_dictionary1,kmer_dictionary2)
    
    anchor_scaffold1 = full_anchor[0]
    anchor_scaffold2 = full_anchor[1]
    position1 =  extract_position_from_dictionary(full_anchor[0])
    position2 = extract_position_from_dictionary(full_anchor[1])
    

    while kmer_number < len(kmer_list) and position1[1]< scaffold_dictionary1[scaffold1][0]+len(scaffold1) and position2[1] < scaffold_dictionary1[scaffold1][0]+ len(scaffold2) :
        kmer_number += 1
        if kmer_number < len(kmer_list) and (anchor_dictionary[kmer_list[kmer_number]][0] < position1[0] or anchor_dictionary[kmer_list[kmer_number]][0] > position1 [1]) and (anchor_dictionary[kmer_list[kmer_number]][0] < position2[0] or anchor_dictionary[kmer_list[kmer_number]][0] > position2 [1]):
            plus_kmer = extend_anchor_full(kmer_list[kmer_number],scaffold1,scaffold2,kmer_dictionary1[kmer_list[kmer_number]][0],kmer_dictionary2[kmer_list[kmer_number]][0],kmer_dictionary1,kmer_dictionary2)
            dictionary_in_dictionary(anchor_scaffold1,plus_kmer[0])
            dictionary_in_dictionary(anchor_scaffold2,plus_kmer[1])
            position1 =  extract_position_from_dictionary(plus_kmer[0])
            position2 = extract_position_from_dictionary(plus_kmer[1])
        
        display_charging_bar((kmer_number/len(kmer_list)*100))
    
    with open("resultat2.txt","w") as file :
        file.write(str(full_anchor))
        file.write("\n")
        file.write(str(position1)+ str(position2)+"\n")
        file.write(str(scaffold_dictionary1[scaffold1][0]+len(scaffold1))+str(scaffold_dictionary2[scaffold2][0]+len(scaffold2)))

    return full_anchor


print(multi_anchor_extension_one_scaffold(scaffold_ancre_assemblage,scaffold_ancre_reference,scaffold_assemblage,scaffold_reference,all_anchors,kmer_ass,kmer_ref))

def multi_scaffold_extension(scaffold_dictionary1 : dict, scaffold_dictionary2 : dict, anchor_dictionary : dict,kmer_dictionary1 : dict,kmer_dictionary2 : dict) :
    extreme_position1 = extract_position_from_dictionary(scaffold_dictionary1)
    extreme_position2 = extract_position_from_dictionary(scaffold_dictionary2)
    scaffold_number1 = 0
    scaffold_number2 = 0
    scaffolds1 = list(scaffold_dictionary1.keys())
    scaffolds2 = list(scaffold_dictionary2.keys())

    current_scaffold  = multi_anchor_extension_one_scaffold(scaffolds1[scaffold_number1],scaffolds2[scaffold_number2],scaffold_dictionary1,scaffold_dictionary2,anchor_dictionary,kmer_dictionary1,kmer_dictionary2)

    absolute_position1 = extract_position_from_dictionary(current_scaffold[0])
    absolute_position2 = extract_position_from_dictionary(current_scaffold[1])

    while absolute_position1 < extreme_position1[1] or absolute_position2 < extreme_position2[1] :
        print("ok")
#multi_scaffold_extension(scaffold_assemblage,scaffold_reference,all_anchors,kmer_ass,kmer_ref)


def supprimer_clef(dictionnaire, clef):
    if clef in dictionnaire:
        del dictionnaire[clef]
    