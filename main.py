import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

l = 500
#seeds =[[None]*2]*l

seeds= np.random.normal(4,6,size=(l,2))
#for ii in range(len(seeds)):
#	seeds[ii] = [random.uniform(0,10),random.uniform(0,10)]
	
tri = Delaunay(np.asarray(seeds))

Fox = np.asarray([3.,3.])
Hunter = np.asarray([2.,2.])


def Topological_Distance(Fox,Hunter,tri): #compute the topological distance
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

def Fox_Move(Fox,dr):#potential position after a move of dr
	Fox_dr = [Fox]*8
	for ii in range(len(Fox_dr)):
		Fox_dr[ii] = [Fox[0]+dr*np.cos(np.pi/4.*ii),
				Fox[1]+dr*np.sin(np.pi/4.*ii)]
	return Fox_dr

def Hunter_Move(Hunter,Fox,dr):#makes the Hunter moves from dr in direction of the Fox
	v = [0,0]
	v[0] = Fox[0]-Hunter[0] #distance vector
	v[1] = Fox[1]-Hunter[1]
	norm = (v[0]**2+v[1]**2)**0.5
	Hunter[0] += dr*v[0]/norm
	Hunter[1] += dr*v[1]/norm
	return Hunter
	
	
def Dist_Euclid(A,B): #euclidean distance
	return np.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2)

def Position_Iteration(Hunter,Fox,tri):#move fox (topological distance) and hunter (past position tracking)
	dr = 5.
	f = 0.1
	Fox_dr = Fox_Move(Fox,dr)
	dist = Topological_Distance(Fox,Hunter,tri)
	dist_dr = [Topological_Distance(Fox_dr[ii],Hunter,tri) for ii in range(8)] 
	while np.max(dist_dr)<=dist and dr<10:
		dr=dr+0.1
		if dr>=10:
			return Hunter,Fox
		Fox_dr = Fox_Move(Fox,dr)

		dist_dr = [Topological_Distance(Fox_dr[ii],Hunter,tri) for ii in range(8)]

	#print("=======\n Fox_dr:\n")
	#print(Fox_dr[0])
	#print(dist_dr)			
	Hunter = Hunter_Move(Hunter,Fox,f*0.75)
	theta = np.pi*np.argmax(dist_dr)/4.
	#print(np.argmax(dist_dr),theta)
	#print("------------------") 
	#print(theta)
	#Fox[0] = Fox[0]+dr*np.cos(theta)
	#Fox[1] = Fox[1]+dr*np.sin(theta)
	Fox = [Fox[0]+f*np.cos(theta),
		Fox[1]+f*np.sin(theta)]
	return Hunter,Fox

plt.figure()
plt.plot(Fox[0],Fox[1],'or',mfc='white')
#plt.plot(Hunter[0],Hunter[1],'ok',mfc='white')

#for i in range(20):
#	Hunter = Hunter_Move(Hunter,Fox,0.1)
for i in range(200):
	print(i)
	Hunter,Fox = Position_Iteration(Hunter,Fox,tri)
	E =Dist_Euclid(Hunter,Fox)
	T = Topological_Distance(Fox,Hunter,tri)
	plt.plot(Fox[0],Fox[1],'or')
	plt.plot(Hunter[0],Hunter[1],'ok')
	print(E,T)
	if E<0.2 : 
		break

#Hunter,Fox = Position_Iteration(Hunter,Fox,tri)
plt.plot(Fox[0],Fox[1],'xr')
#plt.plot(Hunter[0],Hunter[1],'xk')
plt.triplot(np.asarray(seeds)[:,0],
	np.asarray(seeds)[:,1],tri.simplices)
plt.xlim([0,15]),plt.ylim([0,15])
plt.show()
