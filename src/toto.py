##################################Script Projet court#################################

# author : Malo Leprohon

# For now, it is just a draft and a testing script

# To implement:
#   -function to get arguments, helps function and interaction with user (only pdb or pdb + molecular
#   dynamics files?)
#   -pbxplore part, waiting for test dataset
#   -create matrixs from pbxplore output to compute MI <---- Doing this
#   -analysis of MI results, weigthed network building (maybe try visualization on PyMol)
#   -other stuff if time allows (probably not)

# Problem to fix:
#   -remove Z from pb sequence (and maintain position data)

import pbxplore as pb
import pandas as pd
import math
import copy

#to modify
def combine_element(
    listi, replacement=True, order=True
):  
    """ Return a list of string representing all possible couple combinations of every element of listi """
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
    """ Compute mutual information between two positions in a PB sequence """
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
                        print(probx, proby, probxy)
                        mutual_information += probxy * math.log(probxy / (probx * proby), 2)
    return mutual_information



# MUST REMOVE Z AT SOME POINT

# This create a dataframe with a sequence in each row.
# Matrix of sequences
table_seq = pd.DataFrame()
mock_pb_seq = "ZZdddfklonbfklmmmmmmmmnopafklnoiaklmmmmmnoopacddddddehkllmmmmngoilmmmmmmmmmmmmnopacdcddZZ"
n = 20
for i in range(n):
    table_seq = pd.concat(
        [table_seq, pd.DataFrame(list(mock_pb_seq), columns=[i])], axis=1
    )

table_seq = table_seq.transpose()
print(table_seq)


# Matrix of frequence of PB at a position
table_pos_freq = pd.DataFrame(index=list("abcdefghijklmnopZ"))
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
    list("abcdefghijklmnopZ"),replacement=True, order=True
)  # See if this list can be extracted from dataframe
table_dbpos_freq = pd.DataFrame(
    index=list(combine_element(list("abcdefghijklmnopZ"), replacement=True, order=True)),
    columns=combine_element(list(range(len(table_seq.columns))), replacement=False, order=False),
)
print(dic_keys)
for i in list(table_dbpos_freq.columns):
    pb_couple_count = {key: 0 for key in dic_keys}
    pos1 = i.split(":")[0]
    pos2 = i.split(":")[1]
    for j in range(len(table_seq.index)):
        comb = [table_seq.iloc[j, int(pos1)], table_seq.loc[j, int(pos2)]]
        comb.sort(key=str.lower)
        comb = ":".join(comb)
        pb_couple_count[comb] = pb_couple_count[comb] + 1/len(table_seq.index)
    table_dbpos_freq[i].update(pd.Series(pb_couple_count))

print(pb_couple_count)
print(table_dbpos_freq)

print(pos_mutual_information(2,3,table_pos_freq, table_dbpos_freq))
