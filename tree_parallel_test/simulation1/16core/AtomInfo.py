
import sys, math

def AtomicWeight(SMILES):


     AtomicWeight = {'H': 1.008,'He': 4.003, \
     'Li': 6.941, 'Be': 9.012, 'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'F': 19.00, 'Ne': 20.18, \
     'Na': 22.99, 'Mg': 24.81, 'Al':26.98, 'Si':28.09, 'P': 30.97, 'S': 32.07, 'Cl': 35.45, 'Ar': 39.95, \
     'K': 39.10, 'Ca': 40.08, 'Sc': 44.96, 'Ti': 47.87, 'V': 50.94, 'Cr': 52.00, 'Mn': 54.94,'Fe': 55.58, 'Co': 58.93,\
     'Ni': 58.69, 'Cu': 63.55, 'Zn': 65.38,'Ga': 69.72, 'Ge': 72.63, 'As': 74.92,'Se': 78.96, 'Br': 79.90, 'Kr': 83.80,\
     'Rb': 85.47}

     if(Element in AtomicWeight):
      return AtomicWeight[Element]
     else:
      print ("We don't have the information about %-s!" % (Element))
      exit()


def AtomicNumElec(Element):

 AtomicNumElec = {'H': 1,'He': 2, \
 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, \
 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, \
 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,'Fe': 26, 'Co': 27,\
 'Ni': 28, 'Cu': 29, 'Zn': 30,'Ga': 31, 'Ge': 32, 'As': 33,'Se': 34, 'Br': 35, 'Kr': 36,\
 'Rb': 37}

 if(Element in AtomicNumElec):
  return AtomicNumElec[Element]
 else:
  print ("We don't have the information about %-s!" % (Element))
  exit()


#var = input()
#print (AtomicWeight(var))
