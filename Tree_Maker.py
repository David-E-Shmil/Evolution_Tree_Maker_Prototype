
#Initial stages of a phylogenetic tree creater. The end goal is to create a phyologenetic tree using minimal packages.
#This code will be updated as progress is made. As of now it can take in multiple sequences and tell you which two are closest.
import Bio
import string
import random
from urllib.request import urlopen
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Seq import Seq
from Bio import Align
import math

#Takes in a DNA string and converts it into an RNA string
def DNA_to_RNA_Converter(x):
    q=""
    for i in range(len(x)):
        if x[i]=="T":
            q+="U"
        else:
            q+=x[i]
    return q

#Takes in a fasta file and converts it to a string
def fasta_to_String(x):
    stringer = open(x, 'r')
    for line in range(x):
        if line.startswith(">"):
            i=0

    next(stringer)
    better_stringer = stringer.read()
    final_string = "".join(line.strip() for line in better_stringer.splitlines())
    final_string2 = DNA_to_RNA_Converter(final_string)
    return final_string2

#This function takes in a bunch of sequences from a fasta file and spits out a a dictionary just containing the sequences and a reference number
def Seq_dict_maker(string):
    seq_list = {}
    temp=""
    seq_temp = ""
    i=0
    for x in string:
        if len(x)!=0:
            if x[0] ==">":
                if temp != "":
                    seq_list[i] = [seq_temp]
                    i+=1
                temp = x
                seq_temp = ""
            else:
                seq_temp += x
    seq_list[i] = seq_temp
    return seq_list

#This function takes in a String and makes a list that just contains the names of each sequence
def Name_list_maker(string):
    temp=""
    name_list = []

    for i in string:
        if len(i)!= 0:
            if i[0]==">":
                temp=i
                name_list.append((temp))
    return name_list



#This function calculates the score between two RNA sequences
def Score_Calculator(RNA1,RNA2,match,mismatch,gap_s,gap_e):

    aligning = Align.PairwiseAligner()
    aligning.match_score = match

    aligning.mismatch_score = mismatch
    aligning.open_gap_score = gap_s

    aligning.extend_gap_score = gap_e

    alignment = aligning.align(RNA1, RNA2)
    score = aligning.score(RNA1, RNA2)
    print(score)
    return score

#This function runs through every single combination of of sequence comparison and gives out the highest number
def The_Math(list,list2):
    temp1=[]
    temp2=[]
    x=float('-inf')
    score=0
    temp3 = []
    vari = 0
    for i in range(len(list)-1):
        for a in range(len(list[i])):
            temp1 = []
            vari = 0
            for z in range(i+1,len(list)):
                for y in range(len(list[z])):
                    score = Score_Calculator(list[i][a], list[z][y], 1, -1, -3, -2)
                    if score > x:
                        x=score
                        temp3.append((i,z))
    temp1.append((i,z))

#this function asks the user for the data needed then calls on all the necessary functions.
#Currently just uses a url fasta file from a previous project
def Tree_Maker():
    #fasta = input("Pleast input a fasta file ")

    #url1 = fasta
    url1 = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/SARS_MERS_coronavirus.raw_sequence.fasta"
    response1 = urlopen(url1)
    covid_gen = response1.read().decode("utf-8", "ignore")
    final_string = (line.strip() for line in covid_gen.splitlines())

    seq_dict = Seq_dict_maker(final_string)

    response2 = urlopen(url1)
    covid_gen2 = response2.read().decode("utf-8", "ignore")
    final_string2 = (line.strip() for line in covid_gen2.splitlines())
    Name_List = Name_list_maker(final_string2)

    The_Math(seq_dict,Name_List)
    #match_score = float(input("Input the match score "))
    #mismatch_score = float(input("Input the mismatch score "))
    #gap_opening_score = float(input("Input the gap opening score "))
    #gap_extending_score = float(input("Input the gap extending score "))

if __name__ == "__main__":
    Tree_Maker()