# This file is part of ccheminfolib, a library for the generation of 
# chemical descriptors and modeling of these descriptors (QSAR/QSSR/QSPR)
#
# Copyright (C) 2011-2015 Jeremy Henle and the Denmark Lab
#
# This library is confidential and should not be distributed until
# the time deemed appropriate (Jeremy H., 2015)
"""This file contains the modeling classes for the Molecule libraries"""
##standard python libs
import sys
import os
import math
import random
import time
from operator import itemgetter
from collections import Counter
from copy import deepcopy
##special imports
#ccheminfolib
from ccheminfolib.cchemlib import datatypes as dt
from ccheminfolib.cchemlib import atomtypes as at
from ccheminfolib.cqsarlib import library as lib
#sklearn
from sklearn.cross_decomposition import PLSRegression as PLS
from sklearn import preprocessing
from sklearn import linear_model
from sklearn.cross_validation import KFold
#pyearth (for MARS)
from pyearth import Earth
#numpy/scipy
import numpy as np
from scipy import linalg as L
#plotting
import matplotlib.pyplot as plt

##modes
PLS_RAW 			= 0
PLS_GRIND 			= 1
MARS_RAW 			= 2
MARS_GRIND			= 3 

##pls modes
COMB = 0
IND = 1

##CV modes
RAN_R2 = 0
Q2_LOO = 1

class Modeler(object):
	"""Superclass for all modeler objects"""
	def __init__(self, library, n_components = 6, scaled=False):
		##placeholder variables
		
		self.best_solution = []
		self.best_Q2 = np.float64(0)
		self.scaled = scaled
		self.n_components = n_components
		##member vars
		self.library = library
		#self.mode = mode
		self.descriptor_ids = self.library.getDescriptorIDs()
		##scalers
		self.esp_scaler = preprocessing.StandardScaler()
		self.vdw_scaler = preprocessing.StandardScaler()
		self.GRIND_scaler = preprocessing.StandardScaler()
		self.observable_scaler = preprocessing.StandardScaler()
		##from the library, we can generate a smaller library of just
		##the training set -- so we don't have to later
		self.training_set = lib.DescriptorLibrary(library.label + '_TS', self.library.n_descriptors)
		self.solutions = []
		#we add the molecules that have observable values
		for molecule in self.library.training_set:
			self.training_set.addMolDescriptors(molecule, self.library.getMolDescriptors(molecule))
			self.training_set.addObservable(molecule, self.library.getObservable(molecule))
		super(Modeler, self).__init__()
	def acceptance_probability(self, old_cost, new_cost, T):
		"""returns the acceptance probability [0, 1] """
		a_p = None
		try:
			a_p = math.exp((new_cost-old_cost)/T)
		except OverflowError:
			a_p = float('inf')
		if a_p > 1:
			return np.float64(1.0)
		else:
			return a_p
	def cost(self, solution, k):
		"""determines the cost of a particular solution. 
		For this method, this is related to the cross-validation of the model
		This method must be implemented on a per-model basis
		"""
		return 0.0
	def neighbor(self, solution, descriptor_id=-1):
		"""determines a neighbor solution to the current solution
		This can be based on descriptor ID, descriptor location, etc.
		##This method must be implemented on a per-model basis
		"""
		return []
	def ran_r2(self, solution):
		pass
	def anneal(self,initial_solution,T=np.float64(1.0), T_min = np.float64(0.00001), alpha = np.float64(0.9999), iter=100, max_T_iter=1000, threshold=0.95, k=2, output_file=None):
		"""actual annealing process"""
		self.old_cost = self.cost(initial_solution, k)
		self.best_solution = deepcopy(initial_solution)
		self.best_cost =self.old_cost
		self.iteration_best_solution = deepcopy(initial_solution)
		self.T = T
		self.T_min = T_min
		self.alpha = alpha
		self.max_T_iter = max_T_iter
		t_iter = 0
		while self.T > self.T_min and t_iter < self.max_T_iter:
			#we iterate over a certain number of solutions per temperature point
			#determine which descriptor to play with
			descriptor_num = random.randint(0,self.n_components-1)
			for x in range(int(iter)):
				#first, get a new solution that is a neighbor to the best current solution
				new_solution = self.neighbor(self.iteration_best_solution, descriptor_num)
				#self.solutions.append(deepcopy(sorted(new_solution)))
				#determine the cost of the new solution
				self.new_cost = self.cost(new_solution, k)
				#now we get the acceptance probability; if new solution better than old solution, ap = 1, always
				#NOTE THAT IN THIS CASE, INCREASE VALUE IN COST IS GOOD. I.E., COST FUNCTION EVALUATES
				#GOODNESS OF SOLUTION, NOT "ENERGY" WHERE A SMALLER VALUE WOULD BE BETTER!
				ap = self.acceptance_probability(self.old_cost, self.new_cost, self.T)
				#should we move to the new solution?
				if ap > random.random():
					#record the new best solution
					self.iteration_best_solution = new_solution
					self.old_cost = self.new_cost
			#after iteration, decrease temperature
			iteration_R2_data = self.ran_r2(self.iteration_best_solution)
			overall_R2_data = self.ran_r2(self.best_solution)
			print "*******************************************"
			print "Iteration best configuration: " + str(self.iteration_best_solution)
			print "CV RMSD: " + str(self.old_cost)
			print "R2 train: " + str(iteration_R2_data[0])
			print "R2 test: " + str(iteration_R2_data[1])
			print "Overall best configuration: " + str(self.best_solution)
			print "CV RMSD: " + str(self.best_cost)
			print "R2 train: " + str(overall_R2_data[0])
			print "R2 test: " + str(overall_R2_data[1])
			print "T: " + str(self.T)
			print "*******************************************"
			print ''
			if output_file != None:
				output_file.write("*******************************************\n")
				output_file.write("Iteration best configuration: " + str(self.iteration_best_solution)+'\n')
				output_file.write("CV RMSD: " + str(self.old_cost)+'\n')
				output_file.write("R2 train: " + str(iteration_R2_data[0])+'\n')
				output_file.write("R2 test: " + str(iteration_R2_data[1])+'\n')
				output_file.write("Overall best configuration: " + str(self.best_solution)+'\n')
				output_file.write("CV RMSD: " + str(self.best_cost)+'\n')
				output_file.write("R2 train: " + str(overall_R2_data[0])+'\n')
				output_file.write("R2 test: " + str(overall_R2_data[1])+'\n')
				output_file.write("T: " + str(self.T)+'\n')
				output_file.write("*******************************************\n")
				output_file.write('\n')
			if self.old_cost > self.best_cost:
				self.best_solution = deepcopy(self.iteration_best_solution)
				self.best_cost = self.old_cost
				if  self.best_cost > threshold:
					
					print "ACCEPTABLE SOLUTION FOUND: " + str(self.best_solution)
					print self.best_cost
					if output_file != None:
						output_file.write("ACCEPTABLE SOLUTION FOUND: " + str(self.best_solution) + '\n')
					break
			self.T = self.T*alpha
			t_iter += 1
		##determine cross-validation value
		self.cross_validation_value = self.cost(self.best_solution, k=k)
		return dt.SUCCESS
	def generate_model(self, solution):
		"""using the training_set, generate a model based on the solution"""
		pass
	def cross_validation(self, solution):
		"""determines Q2 using random_split as the method of choice"""
		return 0.0
		
	def predict_for_single_molecule(self, molecule_ID):
		"""predicts a single molecule based on the current best solution
		Note: If run before anneal(), this uses a random solution
		Must be implemented for each model system
		"""
		return 0.0
	def predict(self, molecule_IDs):
		"""returns predicted observable for a list of molecules"""
		predictions = {}
		#model = self.generate_model(self.current_best_prediction)
		for molecule in self.library.molecules:
			#prediction = []
			predicted_value = self.predict_for_single_molecule(self.library.molecules[molecule].label)
			predictions[molecule] = predicted_value
			#prediction.append(self.library.getObservable()
		return predictions

class GRINDModeler(Modeler):
	"""This modeler uses GRIND descriptors as inputs to a PLS model
		Model generation uses sklearn.linear_model.PLS 
		as alternative to recoding PLS
	"""
	def __init__(self, library, n_components=6,model_type=PLS_GRIND, mode=RAN_R2, threshold=0.6, scaled=False):
		##how many components should PLS keep?
		#self.n_components=n_components
		self.model_type=model_type
		self.mode = mode
		self.library = deepcopy(library)
		self.n_grind = self.library.n_descriptors
		self.split_training_sets = None
		self.split_test_sets = None
		
		
		super(GRINDModeler, self).__init__(library, n_components = n_components, scaled=scaled)
		self.fit_scalers()
	def fit_scalers(self):
		#grind
		GRIND = []
		for molecule in self.library.mol_descriptors:
			for x in range(len(self.library.mol_descriptors[molecule])):
				GRIND.append(self.library.mol_descriptors[molecule][x])
		GRIND = np.array(GRIND).reshape(-1,1)
		self.GRIND_scaler.fit(GRIND)
		observables = []
		for x in range(-99, 100):
			observables.append(math.log((100.0+float(x))/(100.0-float(x))))
		#for molecule in self.library.training_set:
		#	observables.append(self.library.getObservable(molecule))
		observables = np.array(observables).reshape(-1,1)
		self.observable_scaler.fit(observables)
		return dt.SUCCESS
	def cost(self, solution, k):
		"""determines the cost (cross-validation step) of a given solution"""
		return self.k_fold_cv(solution,k)#, k)
	def neighbor(self, solution, descriptor_id=-1):
		neighbor_solution = []
		#determine which descriptor of the solution to change
		descriptor_num = None
		if descriptor_id == -1:
			descriptor_num = random.randint(0,len(solution)-1)
		else:
			descriptor_num = descriptor_id
		acceptable_neighbor_found = False
		potential_neighbor = None
		while not acceptable_neighbor_found:
			neighbor_solution = deepcopy(solution)
			#get new descriptor id for the descriptor to add to the solution
			potential_neighbor = random.randint(0, self.n_grind-1)
			#make sure we don't already have that one
			if potential_neighbor in solution:
				continue
			else:
				pass
			##need to see if this is all zeros
			allZeros = True
			for molecule in self.library.training_set:
				if self.library.getMolDescriptor(molecule, potential_neighbor) != 0.0:
					allZeros = False
				else:
					pass
			if allZeros:
				continue
			else:
				pass
			neighbor_solution[descriptor_num] = potential_neighbor
			if sorted(neighbor_solution) in self.solutions:
				continue
			else:
				acceptable_neighbor_found = True
				return neighbor_solution
		#return the solution
		return neighbor_solution
	def generate_initial_solution(self):
		initial_solution = []
		#rand_int = None
		while len(initial_solution) < self.n_components:
			rand_int = random.randint(0,self.n_grind-1)
			if rand_int in initial_solution:
				continue
			else:
				initial_solution.append(rand_int)
		return initial_solution
	def get_config_training_set_input(self, solution, split_list_id=0): ##can be overloaded if scaling is 100% needed
		"""returns a list [X,Y]"""
		X_TS = [] ##descriptor array
		Y_TS = [] ## observable list(1D)
		X_ts = [] ##descriptor array, test set
		Y_ts = [] ##observable list, test set
		
		if self.split_training_sets == None:
			print "ERROR -- Split datasets not defined!"
			return dt.FAIL
		if self.scaled:
			for catalyst in self.library.training_set:
				if catalyst in self.split_training_sets[split_list_id]:
					Y_TS.append(self.observable_scaler.transform(np.array([self.library.getObservable(catalyst)]).reshape(-1,1)).ravel()[0])
				else:
					Y_ts.append(self.observable_scaler.transform(np.array([self.library.getObservable(catalyst)]).reshape(-1,1)).ravel()[0])
				descriptors = []
				for descriptor in solution:
					cat_descriptor = self.library.getMolDescriptor(catalyst, descriptor)
					descriptors.append(cat_descriptor)
					del cat_descriptor
				if catalyst in self.split_training_sets[split_list_id]:
					X_TS.append(self.GRIND_scaler.transform(np.array(descriptors).reshape(1,-1)).ravel())
				else:
					X_ts.append(self.GRIND_scaler.transform(np.array(descriptors).reshape(1,-1)).ravel())
		else:
			for catalyst in self.library.training_set:
				if catalyst in self.split_training_sets[split_list_id]:
					Y_TS.append(self.library.getObservable(catalyst))
				else:
					Y_ts.append(self.library.getObservable(catalyst))
				descriptors = []
				for descriptor in solution:
					cat_descriptor = self.library.getMolDescriptor(catalyst, descriptor)
					descriptors.append(cat_descriptor)
					del cat_descriptor
				if catalyst in self.split_training_sets[split_list_id]:
					X_TS.append(deepcopy(descriptors))
				else:
					X_ts.append(deepcopy(descriptors))
		return X_TS, Y_TS, X_ts, Y_ts
	def set_splits(self, percentage=0.7, num_splits=1):
		#determine which catalysts end up in training set, and which end up in test set
		self.n_split = num_splits
		self.split_training_sets = []
		self.split_test_sets = []
		for x in range(self.n_split):
			self.split_training_sets.append([])
			self.split_test_sets.append([])
			while len(self.split_training_sets[x]) < math.floor(len(self.library.training_set)*percentage):
				candidate = self.library.training_set[random.randint(0,len(self.library.training_set)-1)]
				if candidate not in self.split_training_sets[x]:
					self.split_training_sets[x].append(candidate)
				else:
					continue
			for molecule in self.library.training_set:
				if molecule not in self.split_training_sets:
					self.split_test_sets[x].append(molecule)
				else:
					continue
		return dt.SUCCESS
	def get_split_data_for_solution(self, solution, num_splits):
		"""returns data split into multiple parts"""
		split_data = []
		for x in range(num_splits): 
			split_data.append([[],[]])
		#print split_data
		num_per_list = len(self.library.training_set)/int(num_splits)
		current_list = 0
		if self.scaled:
			for molecule in self.library.training_set:
				descriptors = []
				if len(split_data[current_list][0]) > num_per_list:
					current_list += 1
				split_data[current_list][1].append(self.observable_scaler.transform(np.array([self.library.getObservable(molecule)]).reshape(-1,1)).ravel()[0])
				for descriptor in solution:
					descriptor = np.array([self.library.getMolDescriptor(molecule, descriptor)]).reshape(-1,1)
					descriptors.append(descriptor)
					#descriptors.append(self.GRIND_scaler.transform(descriptor).ravel()[0])
				split_data[current_list][0].append(self.GRIND_scaler.transform(np.array(descriptors).reshape(1,-1)).ravel())
				del descriptors
		else:
			for molecule in self.library.training_set:
				descriptors = []
				if len(split_data[current_list][0]) > num_per_list:
					current_list += 1
				split_data[current_list][1].append(self.library.getObservable(molecule))
				for descriptor in solution:
					descriptors.append(self.library.getMolDescriptor(molecule, descriptor))
				split_data[current_list][0].append(deepcopy(descriptors))
				del descriptors
		return split_data
	def generate_model(self, solution):

		#prepare the dater
		X_TS, Y_TS, X_ts, Y_ts = self.get_config_training_set_input(solution, 0)
		X = np.vstack((X_TS,X_ts))
		Y = np.hstack((Y_TS,Y_ts))
		#print "Y: " + str(Y)
		model = None
		if self.model_type == PLS_GRIND:
			model = PLS(n_components=len(solution), scale=False)
			try:
				model.fit(X,Y)
			except:
				print "ERROR -- Model generation failed!"
				return dt.FAIL
		elif self.model_type == MARS_GRIND:
			try:
				model = Earth()
				model.fit(X,Y)
			except:
				print "ERROR -- Model generation failed!"
				return dt.FAIL
		else:
			print "ERROR -- Unknown GRIND model mode selected!"
			return dt.FAIL
		return model
	def ran_r2(self, solution):
		R2 = []
		X_TS,Y_TS,X_ts,Y_ts = self.get_config_training_set_input(solution, 0)
		X_TS = np.array(X_TS)
		Y_TS = np.array(Y_TS)
		X_ts = np.array(X_ts)
		Y_ts = np.array(Y_ts)
		#new_training_set = [[],[]]
		#prediction_set = [[],[]]
		#random_index = []
		#num_new_ts = int(math.ceil(len(data[0])*0.7))
		model = None
		##generate the model
		if self.model_type == PLS_GRIND:
			model = PLS(n_components=len(solution), scale=False)
			model.fit(X_TS,Y_TS)
		elif self.model_type == MARS_GRIND:
			model = Earth()
			model.fit(np.array(X_TS), np.array(Y_TS))
		else:
			return dt.FAIL
		##get the predicted values for the prediction set
		y_hat = model.predict(X_ts)
		reg = linear_model.LinearRegression()
		#print prediction_set[1]
		#print y_hat
		##test
		Y = Y_ts.reshape(-1,1)
		#for x in range(len(Y_ts)):
		#	Y.append([Y_ts[x]])
		#Y = np.array(Y)	
		y_hat = np.array(y_hat)
		reg.fit(Y,y_hat)
		R2 = reg.score(Y, y_hat)
		##figure out R2 for the rest of the training set
		##determine R2 for the training set
		
		#print Y_TS
	#	print Y_ts
		Y_stack = np.hstack((Y_TS,Y_ts))
		#print Y_stack
	#	print Y_stack.reshape(-1,1)
		#sys.exit()
		#model.fit(X_TS,Y_TS)
		y_hat_tot = model.predict(X_TS)
		
		#y_hat_ts = model.predict(new_training_set[0])
		#Y_ts  = []
		#Y_fs = np.vstack(Y_TS,Y_ts)
		Y_vs = Y_TS.reshape(-1,1)
		#for x in range(len(Y_stack)):
			#Y_vs.append(Y_stack[x])
		#Y_vs = np.array(Y_vs)
		reg.fit(Y_vs, y_hat_tot)
		R2_ts = reg.score(Y_vs,y_hat_tot)
		return [R2_ts, R2]
		
	def cross_validation(self, solution):
		if self.mode == RAN_R2:
			R2_summ = 0.0
			for x in range(self.n_split):
			
				X_TS,Y_TS,X_ts,Y_ts = self.get_config_training_set_input(solution, x)
				X_TS = np.array(X_TS)
				Y_TS = np.array(Y_TS)
				X_ts = np.array(X_ts)
				Y_ts = np.array(Y_ts)
				#new_training_set = [[],[]]
				#prediction_set = [[],[]]
				#random_index = []
				#num_new_ts = int(math.ceil(len(data[0])*0.7))
				model = None
				##generate the model
				if self.model_type == PLS_GRIND:
					model = PLS(n_components=len(solution), scale=False)
					model.fit(X_TS,Y_TS)
				elif self.model_type == MARS_GRIND:
					model = Earth()
					model.fit(X_TS, Y_TS)
				else:
					return dt.FAIL
				##get the predicted values for the prediction set
				y_hat = model.predict(X_ts)
				reg = linear_model.LinearRegression()
				#print prediction_set[1]
				#print y_hat
				##test
				Y = Y_ts.reshape(-1,1)
				#for x in range(len(Y_ts)):
				#	Y.append([Y_ts[x]])
				#Y = np.array(Y)	
				y_hat = np.array(y_hat)
				reg.fit(Y,y_hat)
				R2 = reg.score(Y, y_hat)
				##figure out R2 for the rest of the training set
				##determine R2 for the full set
				X_stack = np.vstack((X_TS,X_ts))
				#print Y_TS
			#	print Y_ts
				Y_stack = np.hstack((Y_TS,Y_ts))
				#print Y_stack
			#	print Y_stack.reshape(-1,1)
				#sys.exit()
				model.fit(X_stack,Y_stack)
				y_hat_tot = model.predict(X_stack)
				
				#y_hat_ts = model.predict(new_training_set[0])
				#Y_ts  = []
				#Y_fs = np.vstack(Y_TS,Y_ts)
				Y_vs = Y_stack.reshape(-1,1)
				#for x in range(len(Y_stack)):
					#Y_vs.append(Y_stack[x])
				#Y_vs = np.array(Y_vs)
				reg.fit(Y_vs, y_hat_tot)
				R2_ts = reg.score(Y_vs,y_hat_tot)
				
				if R2_ts < R2:
					R2_summ += R2_ts
				else:
					R2_summ += R2
			return R2_summ / self.n_split
		else:
			return 0.0
	def k_fold_cv(self, solution, k):
		sum_observables = 0.0
		sum_error_squared = 0.0
		SD = 0.0
		average = 0.0
		data =	self.get_split_data_for_solution(solution, k)
		full_observables = []
		for molecule in self.library.training_set:
			observable = None
			if self.scaled:
				observable = self.observable_scaler.transform(np.array([self.library.getObservable(molecule)]).reshape(-1,1)).ravel()[0]
			else:
				observable = self.library.getObservable(molecule)
			sum_observables += observable
			full_observables.append(observable)
		#print "Sum of observable data: " + str(sum_observables)
		average = sum_observables / float(len(self.library.training_set))
		RMSD_summ = 0.0
		#X_TS,Y_TS,X_ts,Y_ts = self.get_config_training_set_input(solution, 0)
		#print Y_ts
		N = []
		for ts_index in range(k):
			#cv numbers
			sum_error_squared = 0.0
			sum_null_model = 0.0
			X_test_set = np.array(data[ts_index][0])
			Y_test_set = np.array(data[ts_index][1])
			#print X_test_set
			X_training_set = None
			Y_training_set = None
			#X_training_set = np.array([])
			#Y_training_set = np.array([])
			for tr_index in range(k):
				if tr_index != ts_index:
					if X_training_set == None:
						X_training_set = np.array(data[tr_index][0])
						Y_training_set = np.array(data[tr_index][1])
					else:
						#print training_set_data.shape
						#print np.array(data[tr_index]).shape
						#print training_set_data
						#training_set_data = np.vstack((training_set_data, np.array(data[tr_index])))
						X_training_set = np.vstack((X_training_set, np.array(data[tr_index][0])))
						#print Y_training_set
						#print data[tr_index][1]
						Y_training_set = np.hstack((Y_training_set, np.array(data[tr_index][1])))
			X_training_set = X_training_set.tolist()
			Y_training_set = Y_training_set.tolist()
			model = None
			# Y_training_set
			##generate the model
			if self.model_type == PLS_GRIND:
				model = PLS(n_components=len(solution), scale=False)
				model.fit(X_training_set,Y_training_set)
			elif self.model_type == MARS_GRIND:
				X_training_set = np.array(X_training_set)
				Y_training_set = np.array(Y_training_set)
				model = Earth()
				model.fit(X_training_set, Y_training_set)
			else:
				return dt.FAIL
			#####################
			
			##now we need to get the data for if the full library was being used
			y_hat = model.predict(X_test_set)
			#find the error in the one we left out
			#print y_hat
			if self.model_type == PLS_GRIND:
				for x in range(len(y_hat)):
					sum_error_squared += (y_hat[x][0]-Y_test_set[x])**2
					sum_null_model += (Y_test_set[x]-average)**2
			else:
				for x in range(len(y_hat)):
					sum_error_squared += (y_hat[x]-Y_test_set[x])**2
			
		#X_TS,Y_TS,X_ts,Y_ts = self.get_config_training_set_input(solution, 0)
		#calculate the things we can calculate right now
		
			###now change the variable name cuz why not
			#PRESS = sum_error_squared
			##calculate SD
			#for x in range(len(self.library.training_set)):
			#	SD += (full_observables[x]-average)**2
			#RMSD = math.sqrt(PRESS/len(self.library.training_set))
			#RMSD_summ += RMSD
			N.append(max(1-(sum_error_squared/sum_null_model), 0))
		return np.mean(N)
	def predict(self, molecule_ID):
		#first generate the model
		model = self.generate_model(self.best_solution)
		if model == dt.FAIL:
			print "ERROR -- Prediction model generation failed!"
			return dt.FAIL
		#get the data for one catalyst
		data = []
		for descriptor in self.best_solution:
			data.append(self.library.getMolDescriptor(molecule_ID, descriptor))
		#print model.predict(np.array(data).reshape(1,-1))[0]
		#sys.exit()
		if self.scaled:
			return model.predict(self.GRIND_scaler.transform(np.array(data).reshape(1,-1)))[0]
		else:
			return model.predict(np.array(data).reshape(1,-1))[0]	
	def save_model_plot(self, filename='test.png'):
		##set filenameif self.is_scaled == True:
		self.plot_filename = filename
		##get the y_hats
		y_hat = {}
		for cat in self.library.training_set:
			pred_val = self.predict(cat)
			#print pred_val
			#sys.exit()
			if self.scaled:
				y_hat[cat] = [self.observable_scaler.inverse_transform(self.predict(cat)), self.library.getObservable(cat)]
			else: #assumes data is prescaled, observable AND input descriptors!
				y_hat[cat] = [self.observable_scaler.inverse_transform(self.predict(cat)), self.observable_scaler.inverse_transform([self.library.getObservable(cat)]).ravel()[0]]
		#y_hat = self.predict(model_config, scaled=scaled)
		##separate the dater
		X = []
		Y = []
		for cat in y_hat:
			X.append(y_hat[cat][0][0])
			Y.append(y_hat[cat][1])
		print X
		print Y
		f = open('model_logs/predicted_vs_observed.csv','w')
		for x in range(len(X)):
			f .write(str(X[x])+','+str(Y[x])+'\n')
		f.close()
		X = np.array(X).reshape(-1,1)
		Y = np.array(Y)
		reg = linear_model.LinearRegression()
		reg.fit(X,Y)
		print reg.coef_
		print reg.intercept_
		print reg.get_params()
		plt.figure(1)
		ax = plt.subplot(111, xlabel='Predicted Selectivity', ylabel='Predicted_Selectivity', title='title')
		plt.plot(X,Y, 'bs', [min(X)[0]-0.1, max(X)[0]+0.1], [reg.predict(min(X)[0]-0.1),reg.predict(max(X)[0]+0.1)], 'b-')
		
		for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
			item.set_fontsize(20)
		#plt.title(sys.argv[2].split('/')[-1].split('.')[0])
		#plt.xlabel('Predicted Selectivity', fontsize=20)
		#plt.ylabel('Observed Selectivity', fontsize=20)
		plt.axis([min(X)[0]-0.25,max(X)[0]+0.25,min(Y)-0.25, max(Y)+0.25])
		##get a regression model on X,Y
		
		plt.text(min(X)+0.5, max(Y)-0.2, 'R2 = ' + str(reg.score(X,Y)))
		#print reg.score(X,Y)
		#plt.show()
		plt.savefig(self.plot_filename)
		plt.clf()	
		
	
		
	
		