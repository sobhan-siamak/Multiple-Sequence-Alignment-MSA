


##### @copy by sobhan siamak

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from Bio.Seq import Seq
# from Bio.Alphabet import Gapped
from collections import defaultdict
from Bio import Seq
from random import randrange, random
import random
from datetime import datetime



start_time = datetime.now()

######Pam250 example
# se = ['S','W','R','C','Q']
# seq1 = se[1]
# seq2 = se[3]
#
pmatrix = matlist.pam250
# print(pmatrix[(seq1, seq2)])
# print("pam is", pmatrix)
# print(pmatrix.values())
# print(pmatrix.keys())




def sequences(name, form):
    bio = []
    se = SeqIO.parse(name,form)
    for i in se:
        bio.append(i.seq)
    return bio


def alignindex(name, form):
    se = SeqIO.parse(name, form)
    dictse = SeqIO.to_dict(se)
    index = []
    for x in dictse.keys():
      index.append(x)
    numbers = []
    for x in index:
        numbers.append(x.split('seq')[1])
    numbers = [int(i) for i in numbers]
    return numbers

def seq2dict(bio):
    #https://stackoverflow.com/questions/51539434/filling-empty-dictionary-in-python
    biodict = defaultdict(list)
    for i in range(20):
        for j in bio[i]:
            biodict[int(i)].append(j)
    return biodict

def seq2list(bio):
    biolist = []
    for i in range(len(bio)):
        biolist.append(list(bio[i]))
    return biolist

def pairwise_align(bio):
    m = np.size(bio)
    poplevel = []
    ####Select random sequence from pairwise alignment
    for i in range(m):
        j = random.choice([x for x in range(0, m) if x != i])
        align = pairwise2.align.globaldx(bio[i], bio[j], pmatrix)
        poplevel.append(align[0][0])



        # pairwise = []
        # for j in range(m):
        #     if (i != j):
        #         align = pairwise2.align.globaldx(bio[i], bio[j], pmatrix)
        #         pairwise.append(align[0][0])  # return first alignment sequence
        # poplevel.append(pairwise[np.random.randint(0, len(pairwise))])
    return poplevel

def random_gap(poplevel, extension):
    mx = max(len(i) for i in poplevel)  ##max size of sequence in all sequences
    mxrate = int(mx * extension)  ##max rate of extension in every sequence
    # print("mxrate is :", mxrate)
    for i in range(len(poplevel)):
        # print("len poplevel is:", len(poplevel))
        # print(len(poplevel[i]))
        while (len(poplevel[i]) < mxrate):
        # num = mxrate-len(poplevel[i])
        # for j in range(mxrate,num, -1):
            poplevel[i] = np.insert(poplevel[i], randrange(0, np.size(poplevel[i])), '-')
        # print(len(poplevel[i]))
        # print(mxrate)
    # print("finish")
    return poplevel




####### Reading BioData by biopython
name1 = "BBA0017_01.tfa"
name2 = "BBA0019_02.tfa"
name3 = "BBA0023_03.tfa"
name4 = "BBA0029_04.tfa"
name5 = "BBA0031_05.tfa"

form1 = "fasta"
form2 = "clustal"

alname1 = "1-align.clustal_num"
alname2 = "2-align.clustal_num"
alname3 = "3-align.clustal_num"
alname4 = "4-align.clustal_num"
alname5 = "5-align.clustal_num"

bio1 = sequences(name1, form1)
bio2 = sequences(name2, form1)
bio3 = sequences(name3, form1)
bio4 = sequences(name4, form1)
bio5 = sequences(name5, form1)
print("bio is:", bio1[0][0])
print(bio1[0])
print(np.size(bio1[0]))
print(len(bio1))
biodict1 = seq2dict(bio1)
print(str(biodict1[0]))
print(biodict1[1])
print(type(bio1[0]))
print(type(biodict1[1]))
a = pairwise2.align.globalxx(bio1[0], bio1[1])
print(list(a[0][0]))




mx = max(len(i) for i in bio3)
print("max len is:", mx)
albio1 = sequences(alname1, form2)
albio2 = sequences(alname2, form2)
albio3 = sequences(alname3, form2)
albio4 = sequences(alname4, form2)
albio5 = sequences(alname5, form2)
print("albio is:", list(albio1[0]))
print(albio2[0])
print(np.shape(albio1))

########Align index for finding and sorting
aindex1 = alignindex(alname1, form2)
aindex2 = alignindex(alname2, form2)
aindex3 = alignindex(alname3, form2)
aindex4 = alignindex(alname4, form2)
aindex5 = alignindex(alname5, form2)
print(aindex1)
print(aindex2)

albiodict1 = seq2dict(albio1)
albiodict2 = seq2dict(albio2)
albiodict3 = seq2dict(albio3)
albiodict4 = seq2dict(albio4)
albiodict5 = seq2dict(albio5)

print("dict albio1 is", pmatrix[(albiodict1[0][94], albiodict1[0][76])])
print("dict albio is:", albiodict1)
print(np.shape(albiodict1[0]))
seq = list(albiodict1.values())
print("seq is:", seq)
alignments = pairwise2.align.globaldx(bio1[0], bio1[17], pmatrix)
print("alignment is:", alignments[0][0])

# print(bio1)
# ins = bio1[0]
#
# # biodict1 ={}  https://stackoverflow.com/questions/51539434/filling-empty-dictionary-in-python
# biodict1 = defaultdict(list)
#
# for i in range(20):
#     for j in bio1[i]:
#         biodict1[int(i)].append(j)
# print(biodict1.get(16))
# biodict1_17 = biodict1.get(17)
# print("pam example is:", pmatrix[(biodict1_17[0],biodict1_17[1])])
#
# ins = bio1[0]
# print(type(ins))
# for i in range(10):
#   ins = np.insert(ins, np.random.randint(0,len(ins)), '-')
# # ins = Seq(str(ins))
# print(ins)
# print(len(ins))
# print("pam example is:", pmatrix[(ins[0],ins[1])])
# print(pmatrix['M', 'H'])



############ Genetic Algorithm #############
extension = 1
psize = 100#population size
pc=0.85#probability of crossover
pm=0.15#probability of mutation
seqnumber = 20# the size of sequnences
iteration = 2500
epsilon = 0.0001
# selectnumber = 70 #### for offspring select from parent selection
# selectnumber2 = 20#### from parent selection for elitism
# selectrandom = 10#### for random parent selection for chance to not best parents
selectnumber = round(0.7*psize)
selectnumber2 = round(0.2*psize)
selectrandom = round(0.1*psize)
##### init population = 100
def initpopulation(bio, extension, psize):
    population = []
    m = np.size(bio)##m=20
    for k in range(psize):
        print("population is", k)
        poplevel = pairwise_align(bio)
        poplevel = seq2list(poplevel)
        ###Add random gap
        # poplevel = list(seq2dict(poplevel).values())
        # poplevel = bio
        poplevel2 = random_gap(poplevel, extension)
        # poplevel2 = poplevel
        population.append(np.array(poplevel2))

    return population

def initpopulation2(bio, extension, psize):
    population = []
    m = np.size(bio)##m=20
    for k in range(psize):
        print("population is", k)
        # poplevel = pairwise_align(bio)
        bio = seq2list(bio)
        poplevel2 = random_gap(bio, extension)
        population.append(np.array(poplevel2))

    return population



####insert random gap
def insert_gap(poplevel):
    mx = max(len(i) for i in poplevel)  ##max size of sequence in all sequences
    # mx = mx + 10
    for i in range(len(poplevel)):
        # print("len poplevel is:", len(poplevel))
        # print(len(poplevel[i]))
        while (len(poplevel[i]) < mx):
             poplevel[i] = np.insert(poplevel[i], randrange(0, np.size(poplevel[i])), '-')
    return poplevel


####### Fitness Function and calculate fitness
######first fitness is sum of pairs
def fitness(poplevel, pmatrix):
    score=0
    m1,n1 = np.shape(poplevel)
    # print("m1 , n1 in fitness are:", m1, n1)
    for i in range(n1):
        for j in range(m1):
            for k in range(j+1, m1):
                  if(poplevel[j,i]=='-' or poplevel[k,i]=='-'):
                      score = score - 10
                  elif(poplevel[j,i]=='-' and poplevel[k,i]=='-'):
                      score = score + 0

                  else:
                    try:
                        score = score + pmatrix[(poplevel[j,i],poplevel[k,i])]
                    except:
                        score = score + pmatrix[(poplevel[k,i], poplevel[j,i])]

    return score
#####second fitness is Percentage of non gaps
def fitness2():
    pass
######## Iteration = 2000

def selection_simple(population, selectnumber, pmatrix):
    scores = []
    parents = []
    for i in range(psize):
        sc = fitness(population[i], pmatrix)
        scores.append(sc)
    sindex = np.argsort(scores)[::-1]
    for i in range(selectnumber):
        parents.append(np.array(population[sindex[i]]))

    return parents


#####Parent Selection Rw
def selectionRW(population, selectnumber, pmatrix):
    scores=[]
    parents=[]
    for i in range(psize):
        sc = fitness(population[i], pmatrix)
        scores.append(sc)
    mn = min(scores)
    newscores = [scores[i]-mn for i in range(len(scores))]
    s = sum(newscores)
    newscores = [(newscores[i]/s) for i in range(len(newscores))]
    newscores = np.cumsum(newscores)
    for i in range(selectnumber):
         q1 = np.random.rand()
         index1 = tuple(np.where(q1 <= newscores))[0]
         ind1 = min(index1)
         q2 = np.random.rand()
         index2 = tuple(np.where(q2 <= newscores))[0]
         ind2 = min(index2)
         # while(ind1==ind2):
         #     q2 = np.random.rand()
         #     index2 = tuple(np.where(q2 <= newscores))[0]
         #     ind2 = min(index2)
         parents.append(population[ind1])
         parents.append(population[ind2])
         #######remove selected parents from population
         # population.remove(population[ind1])
         # population.remove(population[ind2])
         # newscores = list(newscores)
         # newscores.remove(newscores[ind1])
         # newscores.remove(newscores[ind2])
         # scores.remove(scores[ind1])
         # scores.remove(scores[ind2])
    return parents










#####Parent Selection Tournoment
def TSelection(population, psize, pmatrix):
    pass
#### CrossOver
def crossover(generation, selectnumber):
    childs = []
    childs2 = []
    for i in range(0, selectnumber, 2):
        parent1 = generation[i]
        parent2 = generation[i + 1]
        child1 = []
        child2 = []
        for j in range(seqnumber):  ####20 = len(parent1)= seqnumber
            if (np.random.rand() >= 0.5):
                child1.append(parent2[j])
                child2.append(parent1[j])
            else:
                child1.append(parent1[j])
                child2.append(parent2[j])
        childs.append(child1)
        childs.append(child2)
    # for k in range(len(childs)):
    #     childs2.append(np.array(list(childs[k])))




    return childs














#### Mutation
def mutation(offspring):
        pm2 = 0.4
        newpopulation = []
        pop = []
        for i in range(len(offspring)):
            # try:
               newpopulation.append(insert_gap(offspring[i]))
            # except:
            #   newpopulation = newpopulation.copy()
        for j in range(len(offspring)):
             poplevel = seq2list(offspring[j])
             ###Add random gap
             # poplevel = list(seq2dict(poplevel).values())
             # poplevel = bio
             poplevel2 = random_gap(poplevel, extension)
             # poplevel2 = poplevel
             pop.append(np.array(poplevel2))

        return pop








#### Survival Selectionoffspring, pm
def survival_selection(population, offspring):
       newpopulation = []
       sc1 = []
       newscore = []
       # offspring = np.array(seq2list(offspring))
       # newpop = all offspring + 30% population=20%best for elitism + 10% random for chance to not best
       #### select all offspring
       for i in range(psize):
           sctemp1 = fitness(population[i], pmatrix)
           sc1.append(sctemp1)
       for j in range(len(offspring)):
           # print("run is:",j)
           sctemp2 = fitness(offspring[j], pmatrix)
           newpopulation.append(offspring[j])
           newscore.append(sctemp2)
       #### select 20 best old population for elitism
       sindex = np.argsort(sc1)[::-1]
       for i in range(selectnumber2):
           newpopulation.append(population[sindex[i]])
           newscore.append(sc1[sindex[i]])
       #### select 10 worst old population
       sindex2 = np.argsort(sc1)
       for i in range(selectrandom):
           newpopulation.append(population[sindex2[i]])
           newscore.append(sc1[sindex2[i]])


       return newpopulation, newscore


#########Genetic Algorithm Call

population = initpopulation(bio1, extension, psize)
# population = initpopulation(bio2, extension, psize)
# population = initpopulation(bio3, extension, psize)
# population = initpopulation(bio4, extension, psize)
# population = initpopulation(bio5, extension, psize)

# score = []
# for i in range(psize):
#     sc = fitness(population[i], pmatrix)
#     score.append(sc)
# print("scores are:", score)
# amx = np.argsort(score)[::-1]
# mxsore = max(score)
# print("sorted index scores are:", amx)
# print("max score is:", mxsore)
# mxm = max([np.size(population[i][1]) for i in range(psize)])
# print(mxm)


#########main loop of GA

for i in range(iteration):

    scr = []
    print("iteration is:", i)
    generation = selectionRW(population, selectnumber, pmatrix)
    # generation = selection_simple(population, selectnumber, pmatrix)
    # generation = TSelection(population, selectnumber, pmatrix)
    # if (np.random.rand()<pc)
    offspring = crossover(generation, selectnumber)
    offspring = mutation(offspring)
    population, scrs = survival_selection(population, offspring)####the size of population must be psize=100
    print("max score in this iteration is:", max(scrs))
    scr.append(max(scrs))
print(scrs)





###########Score of Multiple Alignment in Website Tool by SOP method

alscore = []
sc1= fitness(np.array(seq2list(albio1)), pmatrix)
alscore.append(sc1)
sc2= fitness(np.array(seq2list(albio2)), pmatrix)
alscore.append(sc2)
sc3= fitness(np.array(seq2list(albio3)), pmatrix)
alscore.append(sc3)
sc4= fitness(np.array(seq2list(albio4)), pmatrix)
alscore.append(sc4)
sc5= fitness(np.array(seq2list(albio5)), pmatrix)
alscore.append(sc5)
print("align scores for albios are:", alscore)









end_time = datetime.now()
print('Time for running this Code is: {}'.format(end_time - start_time))

