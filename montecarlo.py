from random import uniform
import numpy as np
from random import sample
from math import floor
from matplotlib import pyplot as plt
#TODO kJ/mol

dim = 3
N = 256
T = 90
Rc = 12
sigma = 3.7
epsilon=5.953e-24           #J 1.4e-5 kJ/mol
density = 400
molar_mass = 16
tol = 0.5
Nb_cycle = 1000

def Calc_box_length(N, density, molar_mass):
    return ((N*molar_mass*1e4)/(6.022*density))**(1./3)


def Place_particule(box_length):
	pos = np.random.rand(3)*box_length
	return pos

def Move_particule(position, box_length):
	x = position[0] + (np.random.rand()- 0.5) 
	y = position[1] + (np.random.rand()- 0.5)
	z = position[2] + (np.random.rand()- 0.5)

	#Conditions de periodicité
	x = x- floor(x/box_length)*box_length
	y = y- floor(y/box_length)*box_length
	z = z- floor(z/box_length)*box_length
	return [x,y,z]

def Calc_distance(Particule, Positions, box_length):
    distance = np.zeros(256)
    i = Particule

    #Conditions de periodicité
    abs_place = abs(Positions - Positions[i,:]) 
    abs_place = abs_place - (abs_place > box_length/2)*box_length

    distance = np.sqrt(abs_place[:,0]**2 + abs_place[:,1]**2 + abs_place[:,2]**2)

    return distance

def Liste_positions_initiales(N, box_length, tol=0.5):
	positions = np.zeros([N,3])

	for i in range(N):
		distance = Calc_distance(i, positions, box_length)
		is_ill_placed = True
		while(is_ill_placed):
			positions[i,:] = Place_particule(box_length)
			distance = Calc_distance(i, positions,box_length)
			is_ill_placed = (distance[:i]>tol).all() & (distance[i+1:]>tol).all()
		#print("Particule num:", i)
	return positions

L = Calc_box_length(N, density, molar_mass)
Positions = Liste_positions_initiales(N, L)



def Calc_energy(Particule,distance_vector,Rc, sigma, epsilon):
	i = Particule
	#Prise en compte du rayon de coupure et de la distance de la particule avec elle même
	distance_vector = distance_vector[np.where(distance_vector<Rc)]
	distance_vector = distance_vector[np.where(distance_vector!= 0)]

	epot = 4*epsilon*(np.sum(((sigma/distance_vector)**12 -(sigma/distance_vector)**6)))
	return epot



def Change_configuration(N, box_length, positions,old_energy, tol, Rc, sigma, epsilon):
	particule_energy = old_energy
	moving_order = sample(range(N),N)
	new_positions = positions.copy()

	for i in moving_order:
		is_ill_placed = True
		while(is_ill_placed):
			new_positions[i,:] = Move_particule(positions[i,:], box_length)
			distance = Calc_distance(i, new_positions,box_length)
			is_ill_placed = (distance[:i]>tol).all() & (distance[i+1:]>tol).all()
			if is_ill_placed:
				new_positions[i,:] = positions[i,:]
		new_energy = Calc_energy(i, distance,Rc,sigma,epsilon)
		if (new_energy < particule_energy[i]):
			particule_energy[i] = new_energy
		else:
			new_positions[i,:]=positions[i,:]

	return [np.sum(particule_energy), new_positions]


#Energie de départ
epot=np.zeros(N) 
for i in range(N):
	distance = Calc_distance(i, Positions, L)
	epot[i]= Calc_energy(i, distance, Rc, sigma,epsilon)
epot_total = np.sum(epot)

epot_evolution = np.zeros(Nb_cycle)
epot_evolution[0] = epot_total
print(epot_total)

old_positions = Positions

for i in range(1,Nb_cycle):
	epot_total, new_positions = Change_configuration(N, L, old_positions,epot, tol, Rc, sigma, epsilon)
	old_positions = new_positions
	epot_evolution[i]=epot_total
	print (epot_total)


plt.plot(epot_evolution[50:])
plt.show()

#for i in range(Nb_cycle):
#	epot_total = Change_configuration(N, L, Positions, tol, Rc, sigma, epsilon)
