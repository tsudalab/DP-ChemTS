from subprocess import Popen, PIPE
from math import *
import random
import numpy as np
import random as pr
from copy import deepcopy
from types import IntType, ListType, TupleType, StringTypes
import itertools
import time
import math


import tensorflow as tf

import argparse
import subprocess
from load_model import loaded_model
from keras.preprocessing import sequence
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import MolFromSmiles, MolToSmiles
import sys
from make_smile import zinc_data_with_bracket_original, zinc_processed_with_bracket
from threading import Thread, Lock, RLock
import threading
from Queue import Queue
#from mpi4py import MPI
from RDKitText import tansfersdf
from SDF2GauInput import GauTDDFT_ForDFT


import sascorer
import pickle
import gzip
import networkx as nx
from rdkit.Chem import rdmolops





class chemical:

    def __init__(self):

        self.position=['&']
    def Clone(self):

        st = chemical()
        st.position= self.position[:]
        return st

    def SelectPosition(self,m):

        self.position.append(m)

    def Getatom(self):
        return [i for i in range(self.num_atom)]

class Node:

    def __init__(self, position = None, parent = None, state = None, nodelock=threading.Lock()):
        self.position = position
        self.parentNode = parent
        self.childNodes = []
        self.child=None
        self.wins = 0
        self.visits = 0
        self.depth=0
        self.expanded=[]
        self.nodeadded=[]
        self.all_posible=[]
        self.generate_smile=[]
        self.node_index=[]
        self.valid_smile=[]
        self.new_compound=[]
        self.nodelock=nodelock
        self.ucb=[]
        self.core_id=[]



    def Selectnode(self):
        self.nodelock.acquire()

        ucb=[]
        for i in range(len(self.childNodes)):
            ucb.append(self.childNodes[i].wins/(self.childNodes[i].visits)+1.0*sqrt(2*log(self.visits)/(self.childNodes[i].visits)))
        m = np.amax(ucb)
        indices = np.nonzero(ucb == m)[0]
        ind=pr.choice(indices)
        s=self.childNodes[ind]
        #print "which thread's ucb:",threading.currentThread().getName()
        #print ucb
        self.nodelock.release()
        return s

    def Addnode(self, m):

        #n = Node(position = m, parent = self, state = s)
        self.nodeadded.remove(m)
        n = Node(position = m, parent = self)
        self.childNodes.append(n)
        return n


    def Update(self, result):
        #self.nodelock.acquire()
        self.visits += 1
        self.wins += result
        #self.nodelock.release()

    def expanded_node(self,model, state, val):
        all_nodes=[]

        end="\n"
        position=[]
        position.extend(state)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))

        get_int=get_int_old
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad= sequence.pad_sequences(x, maxlen=82, dtype='int32',
            padding='post', truncating='pre', value=0.)

        for i in range(1):
            global graph
            with graph.as_default():
                predictions=model.predict(x_pad)
                #print "shape of RNN",predictions.shape
                preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
                preds = np.log(preds) / 1.0
                preds = np.exp(preds) / np.sum(np.exp(preds))
                next_probas=np.argsort(preds)[-5:]
		#print next_probas
		#next_probas = np.random.multinomial(1, preds, 1)
                #next_int=np.argmax(next_probas)
                #get_int.append(next_int)
                #all_nodes.append(next_int)

        #all_nodes=list(set(all_nodes))
        all_nodes=next_probas
	print all_nodes

        self.expanded=all_nodes
        #print self.expanded




    def expanded_node1(self, model,state,val):
        all_nodes=[]

        end="\n"
        position=[]
        position.extend(state)
        total_generated=[]
        new_compound=[]
        get_int_old=[]
        for j in range(len(position)):
            get_int_old.append(val.index(position[j]))

        get_int=get_int_old
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad= sequence.pad_sequences(x, maxlen=82, dtype='int32',
            padding='post', truncating='pre', value=0.)

        for i in range(60):
            global graph
            with graph.as_default():
                predictions=model.predict(x_pad)
                #print "shape of RNN",predictions.shape
                preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
                preds = np.log(preds) / 1.0
                preds = np.exp(preds) / np.sum(np.exp(preds))
                next_probas = np.random.multinomial(1, preds, 1)
                next_int=np.argmax(next_probas)
                #get_int.append(next_int)
                all_nodes.append(next_int)

        all_nodes=list(set(all_nodes))
        #print all_nodes

        self.expanded=all_nodes
        #print self.expanded

    def node_to_add(self, all_nodes,val):
        added_nodes=[]
        for i in range(len(all_nodes)):
            added_nodes.append(val[all_nodes[i]])

        self.nodeadded=added_nodes
        #print "childNodes:",self.nodeadded





"""Define some functions used for RNN"""



def chem_kn_simulation(model,state,val,added_nodes):
    all_posible=[]

    end="\n"

    position=[]
    position.extend(state)
    position.append(added_nodes)
    total_generated=[]
    new_compound=[]
    get_int_old=[]
    for j in range(len(position)):
        get_int_old.append(val.index(position[j]))

    get_int=get_int_old

    x=np.reshape(get_int,(1,len(get_int)))
    x_pad= sequence.pad_sequences(x, maxlen=82, dtype='int32',
        padding='post', truncating='pre', value=0.)
    while not get_int[-1] == val.index(end):
        predictions=model.predict(x_pad)
        #print "shape of RNN",predictions.shape
        preds=np.asarray(predictions[0][len(get_int)-1]).astype('float64')
        preds = np.log(preds) / 1.0
        preds = np.exp(preds) / np.sum(np.exp(preds))
        next_probas = np.random.multinomial(1, preds, 1)
        next_int=np.argmax(next_probas)
        a=predictions[0][len(get_int)-1]
        next_int_test=sorted(range(len(a)), key=lambda i: a[i])[-10:]
        get_int.append(next_int)
        x=np.reshape(get_int,(1,len(get_int)))
        x_pad = sequence.pad_sequences(x, maxlen=82, dtype='int32',
            padding='post', truncating='pre', value=0.)
        if len(get_int)>82:
            break
    total_generated.append(get_int)
    all_posible.extend(total_generated)


    return all_posible






def predict_smile(all_posible,val):
    new_compound=[]
    for i in range(len(all_posible)):
        total_generated=all_posible[i]

        generate_smile=[]

        for j in range(len(total_generated)-1):
            generate_smile.append(val[total_generated[j]])
        generate_smile.remove("&")
        new_compound.append(generate_smile)

    return new_compound


def make_input_smile(generate_smile):
    new_compound=[]
    for i in range(len(generate_smile)):
        middle=[]
        for j in range(len(generate_smile[i])):
            middle.append(generate_smile[i][j])
        com=''.join(middle)
        new_compound.append(com)
    return new_compound


def simulator(new_compound):
    score=[]
    for i in range(len(new_compound)):
        try:
            m = Chem.MolFromSmiles(str(new_compound[i]))
        except:
            print None
        if m!=None and len(new_compound[i])<=81:
            logp=Descriptors.MolLogP(m)
            t=random.randint(100,300)
            #time.sleep(t)
            for i in range(t):
                i=i+1
        else:
            logp=-10
            #t=random.randint(10,50)
            t=random.randint(100,300)
            #time.sleep(t)
            for i in range(t):
                i=i+1
        score.append(logp)

    return score



def ChemTS_run(rootnode,chem_model,num_simulations):
    """----------------------------------------------------------------------"""
    """----------------------------------------------------------------------"""
    maxnum=0
    #global gau_file_index
    all_compounds=[]
    wave_compounds=[]
    gau_file_index=0
    wave=[]
    depth=[]
    start_time=time.time()
    #while maxnum<1000:
    while time.time()-start_time<3600:
        node = rootnode
        state=['&']
        """selection step"""


        while node.expanded!=[] and node.nodeadded==[]:
            node = node.Selectnode()
            state.append(node.position)

        print "state position:",state
        depth.append(len(state))
        """this if condition makes sure the tree not exceed the maximum depth"""
        if len(state)>=81:
            re=-2.0
            #lock.acquire()
            while node != None:
                node.Update(re)
                node = node.parentNode
            #lock.release()
        else:
            """expansion step"""
            if node.expanded==[]:
	        maxnum+=1
                node.expanded_node(chem_model,state,val)
                node.node_to_add(node.expanded,val)
                if node.nodeadded!=[]:
                    m=random.choice(node.nodeadded)
                    node=node.Addnode(m)

            else:
                if node.nodeadded!=[]:
                    m=random.choice(node.nodeadded)
                    node=node.Addnode(m)


            #maxnum+=1
            data=gaussion_workers(chem_model,val,state,m)
            """simulation step"""

            if data[0]!=-1000:
                re=(0.8*data[0])/(1.0+abs(0.8*data[0]))
                wave_compounds.append(data[1])
                wave.append(data[0])
            else:
                re=-1

            all_compounds.append(data[1])

            """backpropation step"""

            while node != None:
                node.Update(re)
                node = node.parentNode

    print "valid_com=",wave_compounds
    print "depth=",depth
    print "score=",wave
    print len(all_compounds)
    print len(wave_compounds)
    print max(wave)
    print len(wave)
    print maxnum


def gaussion_workers(chem_model,val,state,m):

    all_posible=chem_kn_simulation(chem_model,state,val,m)
    generate_smile=predict_smile(all_posible,val)
    new_compound=make_input_smile(generate_smile)
    score=[]
    kao=[]

    try:
        m = Chem.MolFromSmiles(str(new_compound[0]))
    except:
        m=None
    #if m!=None and len(task[i])<=81:
    if m!=None:
        try:
           logp=Descriptors.MolLogP(m)
        except:
           logp=-1000
        SA_score = -sascorer.calculateScore(MolFromSmiles(new_compound[0]))
        cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(MolFromSmiles(new_compound[0]))))
        if len(cycle_list) == 0:
            cycle_length = 0
        else:
            cycle_length = max([ len(j) for j in cycle_list ])
        if cycle_length <= 6:
            cycle_length = 0
        else:
            cycle_length = cycle_length - 6
        cycle_score = -cycle_length
            #print cycle_score
            #print SA_score
            #print logp
        SA_score_norm=(SA_score-SA_mean)/SA_std
        logp_norm=(logp-logP_mean)/logP_std
        cycle_score_norm=(cycle_score-cycle_mean)/cycle_std
        score_one = SA_score_norm + logp_norm + cycle_score_norm
        score.append(score_one)

    else:
        score.append(-1000)
    score.append(new_compound[0])

    return score




if __name__ == "__main__":
    smile_old=zinc_data_with_bracket_original()
    val,smile=zinc_processed_with_bracket(smile_old)
    print val
    val2=['C', '(', ')', 'c', '1', '2', 'o', '=', 'O', 'N', '3', 'F', '[C@@H]', 'n', '-', '#', 'S', 'Cl', '[O-]', '[C@H]', '[NH+]', '[C@]', 's', 'Br', '/', '[nH]', '[NH3+]', '4', '[NH2+]', '[C@@]', '[N+]', '[nH+]', '\\', '[S@]', '5', '[N-]', '[n+]', '[S@@]', '[S-]', '6', '7', 'I', '[n-]', 'P', '[OH+]', '[NH-]', '[P@@H]', '[P@@]', '[PH2]', '[P@]', '[P+]', '[S+]', '[o+]', '[CH2-]', '[CH-]', '[SH+]', '[O+]', '[s+]', '[PH+]', '[PH]', '8', '[S@@+]']
    logP_values = np.loadtxt('logP_values.txt')
    SA_scores = np.loadtxt('SA_scores.txt')
    cycle_scores = np.loadtxt('cycle_scores.txt')
    SA_mean =  np.mean(SA_scores)
    #print len(SA_scores)

    SA_std=np.std(SA_scores)
    logP_mean = np.mean(logP_values)
    logP_std= np.std(logP_values)
    cycle_mean = np.mean(cycle_scores)
    cycle_std=np.std(cycle_scores)
    chem_model=loaded_model()
    graph = tf.get_default_graph()
    chemical_state = chemical()

    num_simulations=3
    thread_pool=[]
    lock=Lock()
    gau_file_index=0

    """initialization of the chemical trees and grammar trees"""
    root=['&']
    rootnode = Node(position= root)
    maxnum=0
    reward_dis=[]
    compounds=[]

    ChemTS_run(rootnode,chem_model,num_simulations)
