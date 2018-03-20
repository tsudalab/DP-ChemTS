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
from mpi4py import MPI
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
        self.random_node=[]
        self.all_posible=[]
        self.generate_smile=[]
        self.node_index=[]
        self.valid_smile=[]
        self.new_compound=[]
        self.nodelock=nodelock
        self.ucb=[]
        self.core_id=[]



    def Selectnode(self):
        #self.nodelock.acquire()

        ucb=[]
        for i in range(len(self.childNodes)):
            ucb.append(self.childNodes[i].wins/(self.childNodes[i].visits)+1.0*sqrt(2*log(self.visits)/(self.childNodes[i].visits)))
        m = np.amax(ucb)
        indices = np.nonzero(ucb == m)[0]
        ind=pr.choice(indices)
        s=self.childNodes[ind]
        #print "which thread's ucb:",threading.currentThread().getName()
        #print ucb
        #self.nodelock.release()
        return s

    def Addnode(self, m):

        #n = Node(position = m, parent = self, state = s)
        self.nodeadded.remove(m)
        n = Node(position = m, parent = self)
        self.childNodes.append(n)
        return n



    def Update(self, result):
        #self.nodelock.acquire()
        #print "update visits:",self.visits
        self.visits += 1
        self.wins += result
        #self.nodelock.release()

    def expanded_node(self, model, state, val):
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
		next_probs=np.argsort(preds)[-5:]
		next_probs=list(next_probs)
                #next_probas = np.random.multinomial(1, preds, 1)
                #next_int=np.argmax(next_probas)
                #get_int.append(next_int)
                #all_nodes.append(next_int)
		

        #all_nodes=list(set(all_nodes))
        #print all_nodes
	if 0 in next_probs:
	   next_probs.remove(0)
	all_nodes=next_probs
        self.expanded=all_nodes
        print self.expanded

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

        #print "childNodes of current node:", self.nodeadded

    def random_node_to_add(self, all_nodes,val):
        added_nodes=[]
        for i in range(len(all_nodes)):
            added_nodes.append(val[all_nodes[i]])

        self.random_node=added_nodes



        #print "node.nodeadded:",self.nodeadded





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



def ChemTS_run(rootnode,result_queue,lock,chem_model):
    """----------------------------------------------------------------------"""
    """----------------------------------------------------------------------"""
    global maxnum
    global gau_file_index
    start_time=time.time()
    #while maxnum<1000:
    while time.time()-start_time<1800:

        node = rootnode
        state=['&']
        """selection step"""
        node_pool=[]
        lock.acquire()
        while node.expanded!=[] and node.nodeadded==[] and len(node.childNodes)==len(node.expanded):
            s=[]
            for i in range(len(node.childNodes)):
                if node.childNodes[i].visits!=0:
                    s.append(i)
                else:
                    pass
            #print len(s),len(node.childNodes)
            if len(s)==len(node.childNodes):
                node = node.Selectnode()
                state.append(node.position)
            else:
                break

        #print "current tree branch:",state
        depth.append(len(state))


        """this if condition makes sure the tree not exceed the maximum depth"""
        if len(state)>=81:
            re=-2.0
            lock.acquire()    
            while node != None:
                node.Update(re)
                node = node.parentNode
            lock.release()
        else:
            """expansion step"""
            if node.expanded==[]:
	        maxnum+=1
                node.expanded_node(chem_model,state,val)
                node.node_to_add(node.expanded,val)
                node.random_node_to_add(node.expanded,val)
                if node.nodeadded!=[]:
                    m=random.choice(node.nodeadded)
                    node=node.Addnode(m)
                else:
                    node=random.choice(node.childNodes)
                    m=node.position
            else:
                if node.nodeadded!=[]:
                    m=random.choice(node.nodeadded)
                    node=node.Addnode(m)

                else:
                    node=random.choice(node.childNodes)
                    m=node.position

            #maxnum+=1
            lock.release()
            """simulation step"""
            lock.acquire()
            dest_core=random.choice(free_core_id)
            free_core_id.remove(dest_core)
            #message=np.array([state,m])
            comm.send([state,m], dest=dest_core, tag=START)
            lock.release()
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            lock.acquire()
	    #tag=status.Get_tag()
	    #if tag==DONE:
	     #  lock.acquire()
	    free_core_id.append(data[2])
	    lock.release()
            #print data
            #print "free_core_id:",free_core_id
            #print "mpi core:",MPI.ANY_SOURCE
            #print "status:",status

            tag = status.Get_tag()
            if tag == DONE:
                lock.acquire()
                all_compounds.append(data[1])
                lock.release()
                if data[0]!=-1000:
                    re=(0.8*data[0])/(1.0+abs(0.8*data[0]))
                    lock.acquire()
                    wave_compounds.append(data[1])
                    wave.append(data[0])
                    lock.release()
                else:
                    re=-1



            lock.acquire()
            """backpropation step"""
            while node!= None:
                #print "node.parentNode:",node.parentNode
                node.Update(re)
                node = node.parentNode
            lock.release()
    result_queue.put([all_compounds,wave_compounds,depth,wave,maxnum])
  
    
def gaussion_workers(chem_model,val):
    while True:
        #comm.send(None,dest=0,tag=READY)
        task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()
        if tag==START:
	    si_time=time.time()
            state=task[0]
            m=task[1]
            all_posible=chem_kn_simulation(chem_model,state,val,m)
            generate_smile=predict_smile(all_posible,val)
            new_compound=make_input_smile(generate_smile)
            score=[]
            kao=[]

            try:
                m = Chem.MolFromSmiles(str(new_compound[0]))
            except:
                m= None
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
            score.append(rank)
	    fi_time=time.time()-si_time
	    print "si_time:",fi_time

            comm.send(score, dest=0, tag=DONE)

        if tag==EXIT:
            #MPI.Abort(MPI.COMM_WORLD)
            break
    comm.send(None, dest=0, tag=EXIT)




if __name__ == "__main__":
    comm=MPI.COMM_WORLD
    size=comm.size
    rank=comm.rank
    status=MPI.Status()
    READY, START, DONE, EXIT = 0, 1, 2, 3

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

    num_simulations=15
    thread_pool=[]
    lock=Lock()
    gau_file_index=0

    """initialization of the chemical trees and grammar trees"""
    root=['&']
    rootnode = Node(position= root)
    maxnum=0
    reward_dis=[]
    all_compounds=[]
    wave_compounds=[]
    depth=[]
    wave=[]
    result=[]
    result_queue=Queue()
    free_core_id=range(1,num_simulations+1)
    if rank==0:
        #comm.Abort()
        for thread_id in range(num_simulations):
            thread_best = Thread(target=ChemTS_run,args=(rootnode,result_queue,lock,chem_model))
            thread_pool.append(thread_best)

        for i in range(num_simulations):
            thread_pool[i].start()

        for i in range(num_simulations):
            thread_pool[i].join()
        for i in range(num_simulations):
            result.append(result_queue.get())
        print result[0][1]
	print result[0][2]
	print result[0][3]
	print len(result[0][0])
	print len(result[0][1])
        #print result[0][1]
        print max(result[0][3])
        print len(result[0][3])
	print result[0][4]
        comm.Abort()   
        for i in range(num_simulations):
            comm.send(None, dest=i+1, tag=EXIT)

    else:
        gaussion_workers(chem_model,val)
