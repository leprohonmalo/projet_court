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

#Change done : replace pb alphabet string by pb.PB.NAMES, added pbxplore processing

import pbxplore as pbx
import pandas as pd
import math
import copy
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import getopt

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "ho:", ["help", "output=", "traj=", "topo="])
    except getopt.GetoptError:
        print("You did not enter correctly each argument : output directory, trajectory file, topology file.\n")
        use()
        sys.exit(2)

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

    path_list = [output_dir, traj_file, topo_file]
    if None in path_list or "" in path_list:
        print("You did not enter correctly each argument : output directory, trajectory file, topology file.\n")
        use()
    print(path_list)
    
    check_path(output_dir, "directory")
    check_path(traj_file, "file")
    check_path(topo_file, "file")
    return path_list


def use():
    print("Pour le moment y'a pas d'aide parce que y'a des trucs plus important à faire.\n")
    sys.exit()

def check_path(path, path_type):
    path_exist = os.path.exists(path)
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
        

def combine_element(
    listi, replacement=True, order=True
):  
    """ Return a list of string representing all possible couple combinations of every element of listi
    parameters:
        -listi : list, a list of elements
        -replacemnt : boolean, define if combinations of the same item are possible
        -order : boolean, define if items oder in a combination matters (for instance ab != ba)
    """
    if replacement:
        lever = 0
    else:
        lever = 1
    list_ini = copy.deepcopy(listi)
    list_fini = []

    if order :
        for i in range(len(list_ini)):
            for j in range(len(list_ini)):
                if (not replacement
                    and list_ini[i] == list_ini[j]):
                        continue
                list_fini.append(str(list_ini[i]) + ":" + str(list_ini[j]))

    else:
        while len(list_ini) > 0:
            for i in range(len(list_ini) - lever):
                list_fini.append(str(list_ini[0]) + ":" + str(list_ini[i + lever]))
            del list_ini[0]

    return list_fini

def pos_mutual_information(pos1, pos2, freq_df, comb_freq_df):
    """ Compute mutual information between two positions in a PB sequence 
    parameters 
        -pos1 : int or str, a position in a sequence of PB to compute MI with pos2
        -pos2 : int or str, a position in a sequence of PB to compute MI with pos1
        -freq_df : pd.Dataframe, a table of PB frequences at a given position
        -comb_freq_df: pd.Dataframe, a table of PB combination frequences at two given position
    """
    comb_pos = str(pos1) + ":" + str(pos2)
    mutual_information = 0
    print(pos1, pos2)
    for i in freq_df.index:
        for j in freq_df.index:
            comb_pb = str(i) + ":" + str(j)
            probx = freq_df.loc[i, pos1]
            proby = freq_df.loc[j, pos2]
            probxy = comb_freq_df.loc[comb_pb, comb_pos]
            print("comb: {}, probx: {}, proby: {}, probxy: {}".format(comb_pb, probx, proby, probxy))
            if probx != 0 and proby != 0 and probxy !=0:
                        mutual_information += probxy * math.log(probxy / (probx * proby), 16)
    return mutual_information

def mk_table_seq(traj_file, topo_file):
    table_seq = pd.DataFrame()
    i = 0
    for chain_name, chain in pbx.chains_from_trajectory(traj_file, topo_file):
        i += 1
        dihedrals = chain.get_phi_psi_angles()
        pb_seq = pbx.assign(dihedrals)
        table_seq = pd.concat(
            [table_seq, pd.DataFrame(list(pb_seq)[2:-2], columns=[i])], axis=1
        )

    table_seq = table_seq.transpose()
    table_seq.columns = list(range(2,len(table_seq.columns) + 2))
    print(table_seq)
    return(table_seq)




# This create a dataframe with a sequence in each row.
# Matrix of sequences
#table_seq = pd.DataFrame()
#mock_pb_seq = "ZZdddfklonbfklmmmmmmmmnopafklnoiaklmmmmmnoopacddddddehkllmmmmngoilmmmmmmmmmmmmnopacdcddZZ"
#n = 20
#for i in range(n):
#    table_seq = pd.concat(
#        [table_seq, pd.DataFrame(list(mock_pb_seq)[2:-2], columns=[i])], axis=1
#    )

#Use pbxplore to extract pb sequences
# This create a dataframe with a sequence in each row.
# Matrix of sequences

def mk_table_seq(traj_file, topo_file):
    table_seq = pd.DataFrame()
    for chain_name, chain in pbx.chains_from_trajectory(traj_file, topo_file):
        dihedrals = chain.get_phi_psi_angles()
        pb_seq = pbx.assign(dihedrals)
        table_seq = pd.concat(
            [table_seq, pd.DataFrame(list(pb_seq)[2:-2], columns=[i])], axis=1
        )

    table_seq = table_seq.transpose()
    table_seq.columns = list(range(2,len(table_seq.columns) + 2))
    print(table_seq)
    return(table_seq)

# Matrix of frequence of PB at a position
def mk_table_pos_freq(table_seq):
    table_pos_freq = pd.DataFrame(index=list("abc"))#pb cons
    for i in range(len(table_seq.columns)):
        table_pos_freq = pd.concat(
            [table_pos_freq, table_seq.iloc[:, i].value_counts(normalize=True)],
            sort=True,
            axis=1,
        )
    table_pos_freq = table_pos_freq.fillna(0)

    print(table_pos_freq)
    return table_pos_freq


# Matrix of frequece of combination of PB at two positions
def mk_table_dbpos_freq(table_seq):

    dic_keys = combine_element(
        list("abc"),replacement=True, order=True #pb cons
    ) 
    table_dbpos_freq = pd.DataFrame(
        index=list(combine_element(list("abc"), replacement=True, order=True)), #pb cons
        columns=combine_element(list(table_seq.columns), replacement=False, order=False),
    )
    print(dic_keys)

    for i in list(table_dbpos_freq.columns):
        pb_couple_count = {key: 0 for key in dic_keys}
        pos1 = i.split(":")[0]
        pos2 = i.split(":")[1]
        for j in range(len(table_seq.index)):
            comb = [table_seq.iloc[j, int(pos1)-2], table_seq.iloc[j, int(pos2)-2]]
            #comb.sort(key=str.lower)
            comb = ":".join(comb)
            pb_couple_count[comb] = pb_couple_count[comb] + 1
        for key in pb_couple_count.keys(): 
            pb_couple_count[key] = pb_couple_count[key] / len(table_seq.index)
        table_dbpos_freq[i].update(pd.Series(pb_couple_count))
    print(pb_couple_count)
    print(pos1, pos2)

    print(table_dbpos_freq)
    return table_dbpos_freq

print(pbx.PB.NAMES)
#parameters = main(sys.argv[1:])
#seq_table = mk_table_seq(parameters[1], parameters[2])
a = ["zzaaazz" for i in range(10)] + ["zzcabzz" for i in range(10)]
table_seq = pd.DataFrame()
i = 0
for mock_seq in a:
    i+=1
    table_seq = pd.concat(
       [table_seq, pd.DataFrame(list(mock_seq)[2:-2], columns=[i])], axis=1
    )
    print(mock_seq)
table_seq = table_seq.transpose()
table_seq.columns = list(range(2,len(table_seq.columns) + 2))
table_pos_freq = mk_table_pos_freq(table_seq)
table_dbpos_freq = mk_table_dbpos_freq(table_seq)

print(table_seq)
#print(table_pos_freq)
#print(table_dbpos_freq)





#print(pos_mutual_information(2,3,table_pos_freq, table_dbpos_freq))

#Matrix of MI
table_mi = pd.DataFrame(columns=table_seq.columns, index=table_seq.columns, dtype="float")
for i in list(table_mi.index):
    for j in list(table_mi.columns):
        if i == j:
            table_mi.loc[i,j] = 1
        elif i > j:
            table_mi.loc[i,j] = pos_mutual_information(j,i,table_pos_freq, table_dbpos_freq)
        else:
            table_mi.loc[i,j] = pos_mutual_information(i,j,table_pos_freq, table_dbpos_freq)
print(table_mi)

sys.exit()

mk_table_dbpos_freq

#fig = plt.figure()
#ax = fig.adsubplot()
#plt.imshow(table_mi, cmap="hot", interpolation="nearest")
#ax.set_xticks(range(len(table_mi.index)))
#ax.set_xticks(range(len(table_mi.columns)))
#ax.set_xticklabels(table_mi.index)
#ax.ser_xticklabels(table_mi.columns)
#plt.show()
plt.figure(figsize=(50,25))
#ax.set_xlabel = ("position")
sns.set(font_scale = 0.7)
mi_heatmap = sns.heatmap(table_mi, annot=False, xticklabels=2, square = True, linewidths=0.5, cmap="coolwarm", cbar_kws={"label" : "Mutual Information"})
mi_heatmap.set_yticklabels(mi_heatmap.get_yticklabels(), rotation=0)
mi_heatmap.set_xticklabels(mi_heatmap.get_xticklabels(), rotation=0)
plt.xlabel=("position")
plt.ylabel=("position")
plt.tittle=("Mutal Information between positions")
plt.show()