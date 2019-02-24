from random import uniform
import numpy as np
from random import sample, random
from math import floor, exp, log
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

N = 256 #Number of particules
T = 90 # Temperature in Kelvin
Rc = 10# Cutoff radius for energy calculation
sigma = 3.7 #Sigma for methane
epsilon=1.38064852e-23 * 120     #  
Nav = 6.02214076e23  #Avogadro number
density = 400 #methane densitty
molar_mass = 16 #methane molar mass
Nb_cycle = 10000 #Number of monteCarlo cycles
k =1.38064852e-23

def Boltzmann_factor(delta_U, T = 90):
	#return Boltzmann factor for energy acceptation 
	k =1.38064852e-23  #   1.38064852e-23 J/K   
	return exp(-delta_U/(k*T))

def Calc_box_length(N, density, molar_mass):
	#Return montecarlo box size 
    return ((N*molar_mass*1e4)/(6.022*density))**(1./3)


def Place_particule(box_length):
	#Random particule positionning for initial state
	pos=[uniform(0,box_length),uniform(0,box_length),uniform(0,box_length)]
	return pos

def Move_particule(position, box_length, dampening_factor = 0.75):
	""" Move the particule from previous position, if the new position
		is outside the box, it is replaced according periodics conditions
		Dampening factor is used to improve configuration acceptation """
	x = position[0] + (uniform(0,1)- 0.5)*dampening_factor
	y = position[1] + (uniform(0,1)- 0.5)*dampening_factor
	z = position[2] + (uniform(0,1)- 0.5)*dampening_factor

	#Conditions de periodicité
	x = x- floor(x/box_length)*box_length
	y = y- floor(y/box_length)*box_length
	z = z- floor(z/box_length)*box_length

	return [x,y,z]

def Calc_distance(Particule, Positions, box_length):
	#Distance calculator inbetween particules
	distance = np.zeros(N)
	i = Particule
	#print(Positions)
    #Conditions de periodicité
	#abs_place = Positions -Positions[i,:]
	#print(abs_place)

	perio = ((Positions - Positions[i,:]) /(L/2)).astype(int) #entre -L et L
	Positions = Positions - perio*L
	#print(Positions)
	abs_place = Positions -Positions[i,:]

	distance = np.sqrt(abs_place[:,0]**2 + abs_place[:,1]**2 + abs_place[:,2]**2)
	return distance

def Liste_positions_initiales(N, box_length, tol=3.5):
	#Initial positions generator
	positions = np.zeros([N,3])

	for i in range(N):
		is_ill_placed = True
		while(is_ill_placed):
			positions[i,:] = Place_particule(box_length)
			distance = Calc_distance(i, positions,box_length)
			is_ill_placed = (distance[:i]<tol).any() or (distance[i+1:]<tol).any()
	return positions

L = Calc_box_length(N, density, molar_mass)
Positions = Liste_positions_initiales(N, L)


print(L)
print("ICIIIIIIIIIII")
#ax.scatter(Positions[:,0], Positions[:,1], Positions[:,2])
#plt.show(False)


def Calc_energy(Particule,distance_vector,Rc, sigma, epsilon):
	#Energy calculation, return result in kJ
	i = Particule
	#Prise en compte du rayon de coupure et de la distance de la particule avec elle même
	distance_vector = distance_vector[np.where(distance_vector<=Rc)]
	distance_vector = distance_vector[np.where(distance_vector!= 0)]

	epot = 4*epsilon*(np.sum(((sigma/distance_vector)**12 -(sigma/distance_vector)**6)))
	return epot



def Change_configuration(N, box_length, positions,old_energy, Rc, sigma, epsilon, acceptation_rate, dampening_factor):
	particule_energy = old_energy
	moving_order = sample(range(N),N)
	new_positions = positions.copy()
	counter1 = 0.; counter2 =0.

	#work needs to be done around this
	if acceptation_rate <0.4:
		dampening_factor = dampening_factor*0.75

	for i in moving_order:
		new_positions[i,:] = Move_particule(positions[i,:], box_length, dampening_factor)
		distance = Calc_distance(i, new_positions,box_length)
		new_energy = Calc_energy(i, distance,Rc,sigma,epsilon)
		if (new_energy < particule_energy[i]):
			particule_energy[i] = new_energy
			counter1 += 1
		else:
			B_fac = Boltzmann_factor(new_energy-old_energy[i])
			if (B_fac>random()):    #Acceptation rate is too high
				particule_energy[i] = new_energy
				counter2 += 1
			else:	
				new_positions[i,:]=positions[i,:]

	acceptation_rate = (counter1+counter2)/N
	Cv = ((old_energy - old_energy.mean())**2).std()/ (k*T**2)

	return old_energy.mean(),Cv,  new_positions, acceptation_rate, dampening_factor



#Energie de départ
epot=np.zeros(N) 
for i in range(N):
	distance = Calc_distance(i, Positions, L)
	epot[i]= Calc_energy(i, distance, Rc, sigma,epsilon)
epot_total = np.sum(epot)

epot_evolution = np.zeros(Nb_cycle)
Cv_evolution = np.zeros(Nb_cycle)
epot_evolution[0] = epot_total #initial state energy, useless
print(epot_total)
acceptation = np.zeros(Nb_cycle)


old_positions = Positions
dampening_factor = 0.50 ; acceptation_rate = 1
for i in range(1,Nb_cycle):
	epot_total, Cv, old_positions, acceptation_rate, dampening_factor = Change_configuration(N, L, old_positions,epot, Rc, sigma, epsilon, acceptation_rate, dampening_factor)
	epot_evolution[i]=epot_total
	Cv_evolution[i] = Cv
	acceptation[i] = acceptation_rate
	#ax.scatter(new_positions[:,0], new_positions[:,1], new_positions[:,2])
	#fig.canvas.draw()
	#ax.cla()

plt.plot(Cv_evolution[20:]*Nav)
plt.xlabel("Montecarlo Cycles")
plt.ylabel("Heat capacity in J/mol.K")
plt.show()

plt.plot(epot_evolution[20:]*1e-3*Nav)
plt.xlabel("Montecarlo Cycles")
plt.ylabel("Energy in kJ/mol")
plt.show()

plt.plot(acceptation[20:])
plt.show()

