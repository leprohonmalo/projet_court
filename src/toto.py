#! /usr/bin/env python3


##################################Script Projet court#################################

# author : Malo Leprohon

# For now, it is just a draft and a testing script

# To implement:
#   -function to get arguments, helps function and interaction with user (only pdb or pdb + molecular
#   dynamics files?) <---- Doing this 
#   -pbxplore part, waiting for test dataset <---- Doing this
#   -create matrixs from pbxplore output to compute MI <---- Done
#   -analysis of MI results, weigthed network building (maybe try visualization on PyMol) <--- heatmap needs axis labels, pymol thing
#   -other stuff if time allows (probably not)

#To do: 
# verify that each argument is correct : --output dir is an existing directory, topo file in an existing, topo file, traj file in an existing, traj file
# Write a true help
# change value for MI between same pos

#Change done : removed some checking print in differents functions, wrote documentation, help, and added commentary 

import os
import sys
import getopt
import pbxplore as pbx
import pandas as pd
from numpy import NaN
import math
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from copy import deepcopy

def args_check(argv):
    """This function collect parameters from the input command line, and check that every needed
        parameter is provided and correct. Then it returns a list of parameter or call use().

        Parameter:
            argv : a list of the different argument in the input commande line

        Output:
            path_list : a list of parameter values (paths)

    """
    # Check that every needed argument are provided.
    try:
        opts, args = getopt.getopt(argv, "ho:", ["help", "output=", "traj=", "topo="])
    except getopt.GetoptError:
        print("You did not enter correctly each argument : output directory, trajectory file, topology file.\n")
        use()
        sys.exit(2)

    # Initiation of parameters and then assignation.
    output_dir = None
    traj_file = None
    topo_file = None

    for opt, arg in opts:
        if opt in ("-h", "--help"):      
            use()                     
        elif opt in ("-o", "--output"):
            output_dir = arg
        elif opt == "--traj":
            traj_file = arg
        elif opt == "--topo":
            topo_file = arg

    # Expand ~/ paths
    path_list = [os.path.expanduser(output_dir), os.path.expanduser(traj_file), os.path.expanduser(topo_file)]
    if None in path_list or "" in path_list:
        print("You did not enter correctly each argument : output directory, trajectory file, topology file.\n")
        use()

    # Check that every path is correct    
    check_path(output_dir, "directory")
    check_path(traj_file, "file")
    check_path(topo_file, "file")
    return path_list


def use():
    """This function prints the usage of this program and ends it.
    """
    print("""\
    This program compute mutual information during a molecular dynamic simulation between positions 
    in protein sequences of Protein Blocks (PB) from the structural alphabet developped by Barnoud et al. 

    Input :

        It requires a trajectory file and topology file in order to computes PB sequences. 

    Output :
    
        This program output four .csv file and one heatmap (pdf format):

        -PB_seq_table.csv is a table containing all PB sequences extracted from the molecular dynamic.

        -PB_pos_freq_table.csv is a table containing frequences of each PB at each position of a protein.

        -PB_dbpos_freq_table.csv is a table containing frequence of each couple compination of PB for a 
        couple of position, and for each couple of position.

        -MI_table.csv is a table containing mutual information between each position of a protein sequence.

        -mi_heatmap.pdf is an heatmap of MI_table.csv.

    Parameters :
        
        -o, --output : Path to the output directory for output .csv files and figures .

        --traj path/to/trajectory : Path to trajectory file.

        --topo path/to/topology : Path to topology file.

    Option :

        -h, --help : Print this message describing usage of this program. Any other parameters will be ignored.

    References :
        Jonathan Barnoud, Hubert Santuz,Pierrick Craveur, Agnel Praveen Joseph, Vincent
        Jallu, Alexandre G. de Brevern, Pierre Poulain, PBxplore: A Tool To Analyze Local Protein
        Structure And Deformability With Protein Blocks
        (bioRxiv 136408; doi: https://doi.org/10.1101/136408).
    """)
    sys.exit()

def check_path(path, path_type):
    """This function check that a path is correct and is corresponding to the expected path_type (file or directory).
        If not it prints a message and call use().

        Parameters:

            path: a string representing a path

            path_type: either "file" or "directory", represents the expected type of path.

    """
    #Boolean for : is the path a file ? is the path a directory?
    path_to_file = os.path.isfile(path)
    path_to_dir = os.path.isdir(path)
    flag = True
    if path_to_file and path_type == "file":
            flag = False
    elif path_to_dir and path_type == "directory":
            flag = False
    if flag:
        print("The path: {} is not a {} or does not exists.\n".format(path, path_type))
        use()

def mk_table_seq(traj_file, topo_file):
    """This function extracts Protein Block sequences from a trajectory file using pbxplore 
    librairy https://pbxplore.readthedocs.io/en/latest/ and stock the sequences in table_seq 
    DataFrame.

    Parameters: 

        traj_file: string, path to a trajectory file.

        topo_file: string, path to a topology file.

    Output: 

        table_seq: a pandas.DataFrame containing each sequence of PB extracted.

    """
    table_seq = pd.DataFrame()
    # Counter for frame for row name
    i = 0
    for chain_name, chain in pbx.chains_from_trajectory(traj_file, topo_file):
        i += 1
        # Get dihedrals angles to assign a PB to each position
        dihedrals = chain.get_phi_psi_angles()
        pb_seq = pbx.assign(dihedrals)
        # pbxplore need a tresholds of 5 positions to a assign a PB to a positions. Consequently the two first
        # and the two last PB of a sequence are Z which undertermined, we remove them.
        table_seq = pd.concat(
            [table_seq, pd.DataFrame(list(pb_seq)[2:-2], columns=[i])], axis=1
        )
    # For some reasons, the sequences are assigned by columns, so the table is transposed at the end
    # to put them in rows
    table_seq = table_seq.transpose()
    table_seq.columns = list(range(2,len(table_seq.columns) + 2))

    return(table_seq)


def mk_table_pos_freq(table_seq):
    """This function create a table of Protein Blocks frequencies at each position in a sequence
    by looking at all sequences in table_seq. 

    Parameter:

        table_seq: a pandas.DataFrame containing each sequence of PB extracted.

    Output:

        table_pos_freq: a pandas.Dataframe containing frequencies of each Protein Blocks at each position in a sequence.        

    """
    table_pos_freq = pd.DataFrame(index=list(pbx.PB.NAMES))
    for i in range(len(table_seq.columns)):
        #Compute each frequence of values in a column
        table_pos_freq = pd.concat(
            [table_pos_freq, table_seq.iloc[:, i].value_counts(normalize=True)],
            sort=True,
            axis=1,
        )
    table_pos_freq = table_pos_freq.fillna(0)

    return table_pos_freq


def mk_table_dbpos_freq(table_seq):
    """This function create a table of couple of Protein Blocks frequencies at each couple of position in a sequence
    by looking at all sequences in table_seq. 

    Parameter:

        table_seq: a pandas.DataFrame containing each sequence of PB extracted.

    Output:

        table_dbpos_freq: a pandas.Dataframe containing frequencies of each couple of Protein Blocks at each couple of position in a sequence. 
    """

    #Define all possible combination of PB to later use them as keys
    dic_keys = combine_element(
        list(pbx.PB.NAMES),replacement=True, order=True
    ) 
    #Define a table with combination of position by pair in columns and combination of PB by pair rows
    table_dbpos_freq = pd.DataFrame(
        index=list(combine_element(list(pbx.PB.NAMES), replacement=True, order=True)),
        columns=combine_element(list(table_seq.columns), replacement=False, order=False),
    )

    for i in list(table_dbpos_freq.columns):
        #Create a dict that will count occurences of combinations of PB for two positions
        pb_couple_count = {key: 0 for key in dic_keys}
        #Extract position from key
        pos1 = i.split(":")[0]
        pos2 = i.split(":")[1]
        for j in range(len(table_seq.index)):
            # Extract PB combination, -2 is because whe use iloc for j and the positions extracted
            #previously start from 2
            comb = [table_seq.iloc[j, int(pos1)-2], table_seq.iloc[j, int(pos2)-2]]
            comb = ":".join(comb)
            pb_couple_count[comb] = pb_couple_count[comb] + 1
        for key in pb_couple_count.keys(): 
            #conversion of counts (absolute frequencies) to relativ frequencies
            pb_couple_count[key] = pb_couple_count[key] / len(table_seq.index)
        table_dbpos_freq[i].update(pd.Series(pb_couple_count))

    return table_dbpos_freq

def combine_element(
    listi, replacement=True, order=True
):  
    """ This function returns a list of string representing all possible couple combinations of every element of listi.
    Elements of a couple are separated by a ":".

    Parameters:

        listi : a list of elements

        replacement : a boolean, define if combinations of the same item are possible

        order : a boolean, define if items oder in a combination matters (for instance ab != ba)
    
    Output:

        list_fini : a list of string representing all possible couple combinations of every element of listi.
    """
    # Flag that indicate the state of the replacement parameter for one part of the function
    if replacement:
        lever = 0
    else:
        lever = 1
    list_ini = deepcopy(listi)
    list_fini = []

    # Part where orders matters (this is the case for PB combination)
    if order :
        for i in range(len(list_ini)):
            for j in range(len(list_ini)):
                #Avoid repetition in a combination, but not needed for PB
                if (not replacement
                    and list_ini[i] == list_ini[j]):
                        continue
                list_fini.append(str(list_ini[i]) + ":" + str(list_ini[j]))

    #The first item of the list here is deleted at each iteration to avoid repetion of 
    # combination with different order (a:c happens, c:a never happens because a is deleted 
    # afterwards)
    else:
        while len(list_ini) > 0:
            for i in range(len(list_ini) - lever):
                list_fini.append(str(list_ini[0]) + ":" + str(list_ini[i + lever]))
            del list_ini[0]

    return list_fini

def pos_mutual_information(pos1, pos2, freq_df, comb_freq_df):
    """ Compute mutual information between two positions in a PB sequence.

    Parameters

        pos1 : int or str, a position in a sequence of PB to compute MI with pos2

        pos2 : int or str, a position in a sequence of PB to compute MI with pos1

        freq_df : a pandas.Dataframe, a table of PB frequences at a given position

        comb_freq_df: a pandas.Dataframe, a table of PB combination frequences at two given position
    
    Output

        mutual_information: a float value representing mutual information.
    """
    comb_pos = str(pos1) + ":" + str(pos2)
    mutual_information = 0
    for i in freq_df.index:
        #Extraction of value from the differents table
        for j in freq_df.index:
            comb_pb = str(i) + ":" + str(j)
            probx = freq_df.loc[i, pos1]
            proby = freq_df.loc[j, pos2]
            probxy = comb_freq_df.loc[comb_pb, comb_pos]
            # No need to compute if a value is equal to 0 (log(0) undefined)
            if probx != 0 and proby != 0 and probxy !=0:
                        mutual_information += probxy * math.log(probxy / (probx * proby), 16)
    return mutual_information

#Matrix of MI
def mk_mi_table(table_pos_freq, table_dbpos_freq):
    """This function computes mutual information between each position of a set of PB sequences 
    and the results in a Dataframe

    Parameters:
    
        table_pos_freq: a pandas.Dataframe containing frequencies of each Protein Blocks at each position in a sequence. 

        table_dbpos_freq: a pandas.Dataframe containing frequencies of each couple of Protein Blocks at each couple of position in a sequence.

    Output:

        table_mi: a pandas.Dataframe containing each mutual information between each postion in a Protein Block sequence.

    """
    table_mi = pd.DataFrame(columns=table_pos_freq.columns, index=table_pos_freq.columns, dtype="float")
    for i in list(table_mi.index):
        for j in list(table_mi.columns):
            # MI(X,X) makes no sense so NaN is assigned
            if i == j:
                table_mi.loc[i,j] = NaN
            elif i > j:
                table_mi.loc[i,j] = pos_mutual_information(j,i,table_pos_freq, table_dbpos_freq)
            else:
                table_mi.loc[i,j] = pos_mutual_information(i,j,table_pos_freq, table_dbpos_freq)
    return table_mi

def mi_heatmap(table_mi, output_dir):
    """ This function plots an heatmap of table_mi and saves it in output_dir.
        
        Parameters:

            table_mi : pandas.Dataframe containing each mutual information between each postion in a Protein Block sequence.

            output_dir : a string representing a path towards an ouput directory where the heatmap will be saved.
    """
    fig = plt.figure(figsize=(50,25))
    sns.set(font_scale = 0.7)
    mi_heatmap = sns.heatmap(table_mi, mask=table_mi.isnull(), annot=False, xticklabels=2, square = True, linewidths=0.5, cmap="coolwarm", cbar_kws={"label" : "Mutual Information"})
    mi_heatmap.set_yticklabels(mi_heatmap.get_yticklabels(), rotation=0)
    mi_heatmap.set_xticklabels(mi_heatmap.get_xticklabels(), rotation=0)
    plt.xlabel=("position")
    plt.ylabel=("position")
    plt.tittle=("Mutal Information between positions")
    fig.savefig(output_dir+"mutual_information_heatmap.pdf")


def main():
    """ The main function of this script, that uses all previous function to extract Protein Block sequences from a
    trajectory and a topology file and computes mutual information between positions of these sequences.
    """
    # Check parameters
    parameters = args_check(sys.argv[1:])
    print("\nExtracting PB sequences from trajectory.\n")
    # Extraction of PB sequences with pbxplore
    table_seq = mk_table_seq(parameters[1], parameters[2])
    print("\nComputing frequences of PB (This may take a few minutes)\n")
    # Computation of frequencies
    table_pos_freq = mk_table_pos_freq(table_seq)
    table_dbpos_freq = mk_table_dbpos_freq(table_seq)
    print("Computing Mutual Information (This may take a few minutes)\n")
    # Computation of Mutual Information
    table_mi = mk_mi_table(table_pos_freq, table_dbpos_freq)
    mi_heatmap(table_mi, parameters[0])
    # Write results in output directory
    table_seq.to_csv(parameters[0]+"PB_seq_table.csv")
    table_pos_freq.to_csv(parameters[0]+"PB_pos_freq_table.csv")
    table_dbpos_freq.to_csv(parameters[0]+"PB_dbpos_freq_table.csv")
    table_mi.to_csv(parameters[0]+"MI_table.csv")
    print("Done!\n")
    sys.exit()

if __name__ == "__main__":
    main()
