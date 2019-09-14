##################################Script Projet court#################################

# author : Malo Leprohon

# For now, it is just a draft and a testing script

#To implement:
#   -function to get arguments, helps function and interaction with user (only pdb or pdb + molecular 
#   dynamics files?)
#   -pbxplore part, waiting for test dataset 
#   -create matrixs from pbxplore output to compute MI <---- Doing this 
#   -analysis of MI results, weigthed network building (maybe try visualization on PyMol)
#   -other stuff if time allows (probably not)

#Problem to fix:
#   -remove Z from pb sequence (and maintain position data)
#   -modify combine element to exclude doublet or not

import pbxplore as pb
import pandas as pd 
import copy

def combine_element(listi, replacement=True): #A changer car bien que aa soit possible, 1,1 doit ne pas l'être
    # si on passe une liste en argument ca va la modifier, donc les passer avec list
    ''' Return a list of string representing all possible couple combinations of every element of listi '''
    if replacement :
        lever = 0
    else:
        lever = 1
    list_ini = copy.deepcopy(listi)
    list_fini = []
    while len(list_ini) > 0:
        for i in range(len(list_ini)-lever):
            list_fini.append(str(list_ini[0]) + ":" + str(list_ini[i+lever]))
        del list_ini[0]
    return list_fini

#MUST REMOVE Z AT SOME POINT

# This create a dataframe with a sequence in each row.
# Matrix of sequences
table_seq = pd.DataFrame()
mock_pb_seq = "ZZdddfklonbfklmmmmmmmmnopafklnoiaklmmmmmnoopacddddddehkllmmmmngoilmmmmmmmmmmmmnopacdcddZZ"
n = 20
for i in range(n):
    table_seq = pd.concat([table_seq, pd.DataFrame(list(mock_pb_seq), columns=[i])], axis=1)

table_seq = table_seq.transpose()
print(table_seq)



# Matrix of frequence of PB at a position 
table_pos_freq = pd.DataFrame(index=list("abcdefghijklmnopZ"))
for i in range(len(table_seq.columns)):
    table_pos_freq = pd.concat([table_pos_freq, table_seq.iloc[:,i].value_counts(normalize = True)], sort=True, axis =1)
table_pos_freq = table_pos_freq.fillna(0)
print(table_pos_freq)


# Matrix of frequece of combination of PB at two positions
dic_keys = combine_element(list("abcdefghijklmnopZ")) # See if this list can be extracted from dataframe
table_dbpos_freq = pd.DataFrame(index=list(combine_element(list("abcdefghijklmnopZ"))), columns=combine_element(list(range(len(table_seq.columns))), replacement=False))
print(table_dbpos_freq)
for i in list(table_dbpos_freq.columns):
    pb_couple_count = {key : 0 for key in dic_keys}
    pos1 = i.split(":")[0]
    pos2 = i.split(":")[1]
    for j in range(len(table_seq.index)):
        comb = [table_seq.iloc[j,int(pos1)], table_seq.loc[j,int(pos2)]]
        comb.sort(key=str.lower)
        comb = ":".join(comb)
        pb_couple_count[comb] = pb_couple_count[comb] + 1
    table_dbpos_freq[i].update(pd.Series(pb_couple_count))

print(pb_couple_count)
print(table_dbpos_freq)
