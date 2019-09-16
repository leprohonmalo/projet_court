##################################Script Projet court#################################

# author : Malo Leprohon

# For now, it is just a draft and a testing script

# To implement:
#   -function to get arguments, helps function and interaction with user (only pdb or pdb + molecular
#   dynamics files?) <---- Doing this
#   -pbxplore part, waiting for test dataset
#   -create matrixs from pbxplore output to compute MI <---- Done
#   -analysis of MI results, weigthed network building (maybe try visualization on PyMol) <--- heatmap needs axis labels, pymol thing
#   -other stuff if time allows (probably not)


import pbxplore as pb
import pandas as pd
import math
import copy
import matplotlib.pyplot as plt
import seaborn as sns

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
                if not replacement:
                    if list_ini[i] == list_ini[j]:
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
    for i in freq_df.index:
        for j in freq_df.index:
            comb_pb = str(i) + ":" + str(j)
            probx = freq_df.loc[i, pos1]
            proby = freq_df.loc[j, pos2]
            probxy = comb_freq_df.loc[comb_pb, comb_pos]
            if probx != 0:
                if proby != 0:
                    if probxy != 0:
                        mutual_information += probxy * math.log(probxy / (probx * proby), 2)
    return mutual_information


# This create a dataframe with a sequence in each row.
# Matrix of sequences
table_seq = pd.DataFrame()
mock_pb_seq = "ZZdddfklonbfklmmmmmmmmnopafklnoiaklmmmmmnoopacddddddehkllmmmmngoilmmmmmmmmmmmmnopacdcddZZ"
n = 20
for i in range(n):
    table_seq = pd.concat(
        [table_seq, pd.DataFrame(list(mock_pb_seq)[2:-2], columns=[i])], axis=1
    )

table_seq = table_seq.transpose()
table_seq.columns = list(range(2,len(table_seq.columns) + 2))
print(table_seq)


# Matrix of frequence of PB at a position
table_pos_freq = pd.DataFrame(index=list("abcdefghijklmnop"))
for i in range(len(table_seq.columns)):
    table_pos_freq = pd.concat(
        [table_pos_freq, table_seq.iloc[:, i].value_counts(normalize=True)],
        sort=True,
        axis=1,
    )
table_pos_freq = table_pos_freq.fillna(0)
print(table_pos_freq)


# Matrix of frequece of combination of PB at two positions
dic_keys = combine_element(
    list("abcdefghijklmnop"),replacement=True, order=True
)  # See if this list can be extracted from dataframe
table_dbpos_freq = pd.DataFrame(
    index=list(combine_element(list("abcdefghijklmnop"), replacement=True, order=True)),
    columns=combine_element(list(table_seq.columns), replacement=False, order=False),
)
print(dic_keys)
for i in list(table_dbpos_freq.columns):
    pb_couple_count = {key: 0 for key in dic_keys}
    pos1 = i.split(":")[0]
    pos2 = i.split(":")[1]
    for j in range(len(table_seq.index)):
        comb = [table_seq.iloc[j, int(pos1)-2], table_seq.iloc[j, int(pos2)-2]]
        comb.sort(key=str.lower)
        comb = ":".join(comb)
        pb_couple_count[comb] = pb_couple_count[comb] + 1/len(table_seq.index)
    table_dbpos_freq[i].update(pd.Series(pb_couple_count))

print(pb_couple_count)
print(table_dbpos_freq)

print(pos_mutual_information(2,3,table_pos_freq, table_dbpos_freq))

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

#fig = plt.figure()
#ax = fig.adsubplot()
#plt.imshow(table_mi, cmap="hot", interpolation="nearest")
#ax.set_xticks(range(len(table_mi.index)))
#ax.set_xticks(range(len(table_mi.columns)))
#ax.set_xticklabels(table_mi.index)
#ax.ser_xticklabels(table_mi.columns)
#plt.show()
plt.figure(figsize=(50,25))
ax.set_xlabel = ("position")
sns.set(font_scale = 0.7)
mi_heatmap = sns.heatmap(table_mi, annot=False, xticklabels=2, square = True, linewidths=0.5, cmap="coolwarm", cbar_kws={"label" : "Mutual Information"})
mi_heatmap.set_yticklabels(mi_heatmap.get_yticklabels(), rotation=0)
mi_heatmap.set_xticklabels(mi_heatmap.get_xticklabels(), rotation=0)
plt.xlabel=("position")
plt.ylabel=("position")
plt.tittle=("Mutal Information between positions")
plt.show()