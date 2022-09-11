import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

l = 20
seeds =[[None]*20]*20

for ii in range(len(seeds)):
	seeds[ii] = [random.uniform(0,10),random.uniform(0,10)]
	
tri = Delaunay(np.asarray(seeds))

Fox = np.asarray([5,5])
Hunter = np.asarray([3,6])


def Distance_Compute(Fox,Hunter,tri):
	idf = tri.find_simplex(Fox)
	idh = tri.find_simplex(Hunter)
	if idh==-1 or idf == -1 :
		return 0
	else:
		xmin,xmax = np.sort([Fox[0],Hunter[0]])
		if xmin==Fox[0]:
			ymin,ymax = Fox[1],Hunter[1]
		else:
			ymin,ymax = Hunter[1],Fox[1]
		Line = [np.arange(xmin,xmax,0.1),
			np.arange(ymin,ymax,0.1*np.sign(ymax-ymin))]	
		#print(Line)
		idl = [tri.find_simplex([Line[0][ii],Line[1][ii]]) for ii in range(min(len(Line[0]),len(Line[1])))]
		idl_clean=np.asarray([x for x in idl if x!=-1])
		return len(np.unique(idl_clean))

def Fox_Move(Fox,dr):
	Fox_dr = [Fox]*8
	for ii in range(len(Fox_dr)):
		Fox_dr[ii] = [Fox[0]+dr*np.cos(np.pi/4.*ii),
				Fox[1]+dr*np.sin(np.pi/4.*ii)]
	return Fox_dr
	
def Dist_Euclid(A,B):
	return np.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2)

def Position_Iteration(Hunter,Fox,tri):
	dr = 0.1
	Fox_dr = Fox_Move(Fox,dr)
	dist = Distance_Compute(Fox,Hunter,tri)
	dist_dr = [Distance_Compute(Fox_dr[ii],Hunter,tri) for ii in range(8)] 
	while np.max(dist_dr)<=dist and dr<10:
		dr=dr+0.1
		Fox_dr = Fox_Move(Fox,dr)
		dist_dr = [Distance_Compute(Fox_dr[ii],Hunter,tri) for ii in range(8)] 	
	Hunter = [Hunter[0] + 0.1*np.arccos((Fox[0]-Hunter[0])/Dist_Euclid(Hunter,Fox)),
	Hunter[1] + 0.1* np.arcsin((Fox[1]-Hunter[1])/Dist_Euclid(Hunter,Fox))]
	Fox = [Fox[0]+0.1*np.cos(np.pi/4.*ii)*np.argmax(dist_dr),
	Fox[1]+0.1*np.sin(np.pi/4.*ii)*np.argmax(dist_dr)]
	return Hunter,Fox

Hunter,Fox = Position_Iteration(Hunter,Fox,tri)
print(Hunter,Fox)



