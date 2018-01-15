

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import AtomInfo
from SDF2xyzV2 import Read_sdf

def tansfersdf(com,index):


    #m2 = Chem.MolFromSmiles('C1CNC1')
    #m2 = Chem.MolFromSmiles('CCNC(CC)C')
    #m2 = Chem.MolFromSmiles('N[C@@H](C)C(=O)O')
    #m2 = Chem.MolFromSmiles('O=C(C1=CC=CC2=C1C3=CC(N)=C2)N(CCN(C)C)C3=O')
    #m2 = Chem.MolFromSmiles('O=C(C1=CC=CC2=C1C3=CC([N+]([O-])=O)=C2)N(CCN(C)C)C3=O')
    m2 = Chem.MolFromSmiles(com)
    AllChem.EmbedMolecule(m2)
    m3 = Chem.AddHs(m2)
    AllChem.EmbedMolecule(m3)
    Chem.MolToMolFile(m3,'CheckMol'+str(index)+'.sdf')
    try:
        opt = AllChem.UFFOptimizeMolecule(m3,maxIters=200)
    except:
        opt=None
    if opt!=None:
        Chem.MolToMolFile(m3,'CheckMolopt'+str(index)+'.sdf')
        SpinMulti=Read_sdf('CheckMolopt'+str(index)+'.sdf')
    else:
        SpinMulti=0

    print SpinMulti




    return SpinMulti




#hi=tansfersdf('O[C@H]1[C@]2(C)COC[C@@H]1C2')
