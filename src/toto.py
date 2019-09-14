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

import pbxplore as pb
import pandas as pd 
import copy

def combine_element(listi):
    # si on passe une liste en argument ca va la modifier, donc les passer avec list
    ''' Return a list of string representing all possible couple combinations of every element of listi '''
    list_ini = copy.deepcopy(listi)
    list_fini = []
    while len(list_ini) > 0:
        for i in range(len(list_ini)):
            list_fini.append(str(list_ini[0]) + ":" + str(list_ini[i]))
        del list_ini[0]
    return list_fini


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
table_pos_freq = pd.DataFrame(index=list("Zabcdefghijklmnop"))
for i in range(len(table_seq.columns)):
    table_pos_freq = pd.concat([table_pos_freq, table_seq.iloc[:,i].value_counts(normalize = True)], sort=True, axis =1)
table_pos_freq = table_pos_freq.fillna(0)
print(table_pos_freq)

test1 = combine_element(list(range(89)))
print(test1)

