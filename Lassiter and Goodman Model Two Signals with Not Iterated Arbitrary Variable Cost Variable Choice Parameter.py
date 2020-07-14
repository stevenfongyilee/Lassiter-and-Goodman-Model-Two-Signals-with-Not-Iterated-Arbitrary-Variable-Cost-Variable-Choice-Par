# Here we allow arbitrary iteration up the hierarchy.

import numpy
numpy.set_printoptions(linewidth = 120)
import matplotlib
# matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from scipy.stats import beta
from scipy.stats import uniform
from scipy.stats import truncnorm
from scipy import integrate
import itertools

###################

###################

def receiver_0_signal_0(h):
	return state_distribution.pdf(h)

def receiver_0_signal_1(h, theta):
	if h < theta:
		return 0.
	else:
		return state_distribution.pdf(h) / (1. - state_distribution.cdf(theta))

def receiver_0_signal_not_1(h, theta):
	if h > theta:
		return 0.
	else:
		return state_distribution.pdf(h) / (state_distribution.cdf(theta))

def pragmatic_sender(n):
	if n not in pragmatic_sender_memo:
		print 'pragmatic_receiver(%s - 1)[0] = %s' % (n, pragmatic_receiver(n - 1)[0])
		pragmatic_receiver_signal_0_h_array = numpy.sum(pragmatic_receiver(n - 1)[0], axis = (0))
		print 'level %s pragmatic_receiver_signal_0_h_array = %s' % (n - 1, pragmatic_receiver_signal_0_h_array)
		pragmatic_receiver_signal_1_theta_1_by_h_array = pragmatic_receiver(n - 1)[1]
		pragmatic_receiver_signal_not_1_theta_1_by_h_array = pragmatic_receiver(n - 1)[2]
		
# note that we have two options here: first we could make the pragmatic sender sensitive
# to just h, second we could make it sensitive to h and theta. of course in Lassiter and
# Goodman's model the first level sender is only sensitive to h. if we want to keep a
# single model all the way up the hierarchy it seems like making the sender only h
# sensitive is the way to go. there is also a question about theta values: do we get them
# from the same level or the previous level? What way to think about it is consistent with
# the initial model? does the initial model (the first level sender) prescribe only one
# way of doing the calculations? Perhaps one way is to think of the receiver at any level
# as taking the background knowledge of h and then applying a literal interpretation which
# simply renormalizes above any given theta. at the bottom level 0 the background
# information is an absolute prior, and at level 2 the background for any given signal is
# the h posterior. Then the level 3 sender thinks of the level 2 receiver as renormalizing
# for any given theta above (or below) the h posterior for the given signal. That can give
# rise to a level 4 receiver who in turn derives a new posterior. Is this any different
# than what we have below for h sensitive?
		
		if pragmatic_sender_type == 'h_sensitive_version_3':
		
			pragmatic_receiver_signal_0_h_array = pragmatic_receiver_signal_0_h_array
			pragmatic_receiver_signal_1_h_array = numpy.sum(pragmatic_receiver(n - 1)[1], axis = (0))
			pragmatic_receiver_signal_not_1_h_array = numpy.sum(pragmatic_receiver(n - 1)[2], axis = (0))

			if choice_parameter != 0.:
				pragmatic_sender_signal_0_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.tile(pragmatic_receiver_signal_0_h_array / (num_states), (num_states, 1))) - cost_of_null_signal))
				pragmatic_sender_signal_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.tile(pragmatic_receiver_signal_1_h_array / (num_states), (num_states, 1))) - cost))
				pragmatic_sender_signal_not_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.tile(pragmatic_receiver_signal_not_1_h_array / (num_states), (num_states, 1))) - (cost + cost_of_not)))
			else:
				pragmatic_sender_signal_0_non_normalized = numpy.ones(num_states)
				pragmatic_sender_signal_1_non_normalized = numpy.ones((num_states, num_states))
				pragmatic_sender_signal_not_1_non_normalized = numpy.ones((num_states, num_states))
			
			print 'pragmatic_sender_signal_0_non_normalized = \n%s' % pragmatic_sender_signal_0_non_normalized
			print 'pragmatic_sender_signal_1_non_normalized = \n%s' % pragmatic_sender_signal_1_non_normalized
			print 'pragmatic_sender_signal_not_1_non_normalized = \n%s' % pragmatic_sender_signal_not_1_non_normalized
			
			denominator_array = pragmatic_sender_signal_0_non_normalized + pragmatic_sender_signal_1_non_normalized + pragmatic_sender_signal_not_1_non_normalized
		
			pragmatic_sender_signal_0_normalized = pragmatic_sender_signal_0_non_normalized / denominator_array
			pragmatic_sender_signal_1_normalized = pragmatic_sender_signal_1_non_normalized / denominator_array
			pragmatic_sender_signal_not_1_normalized = pragmatic_sender_signal_not_1_non_normalized / denominator_array
		
			print 'level %s pragmatic_sender_signal_0_normalized = \n%s' % (n, pragmatic_sender_signal_0_normalized)
			print 'level %s pragmatic_sender_signal_1_normalized = \n%s' % (n, pragmatic_sender_signal_1_normalized)
			print 'level %s pragmatic_sender_signal_not_1_normalized = \n%s' % (n, pragmatic_sender_signal_not_1_normalized)

		elif pragmatic_sender_type == 'h_sensitive':
		
			pragmatic_receiver_signal_0_pre_normalization_factor_array = numpy.ones(num_states)
			pragmatic_receiver_signal_1_pre_normalization_factor_array = numpy.tile(numpy.sum(pragmatic_receiver(n - 1)[1], axis = (1)), (num_states, 1)).T
			pragmatic_receiver_signal_not_1_pre_normalization_factor_array = numpy.tile(numpy.sum(pragmatic_receiver(n - 1)[2], axis = (1)), (num_states, 1)).T
			
			print 'pragmatic_receiver_signal_0_h_array = \n%s' % pragmatic_receiver_signal_0_h_array
			print 'pragmatic_receiver_signal_1_theta_1_by_h_array = \n%s' % pragmatic_receiver_signal_1_theta_1_by_h_array
			print 'pragmatic_receiver_signal_not_1_theta_1_by_h_array = \n%s' % pragmatic_receiver_signal_not_1_theta_1_by_h_array
						
			print 'pragmatic_receiver_signal_0_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_0_pre_normalization_factor_array
			print 'pragmatic_receiver_signal_1_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_1_pre_normalization_factor_array
			print 'pragmatic_receiver_signal_not_1_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_not_1_pre_normalization_factor_array
						
			if choice_parameter != 0.:
				pragmatic_sender_signal_0_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_0_h_array / pragmatic_receiver_signal_0_pre_normalization_factor_array) - cost_of_null_signal))
				print 'pragmatic_receiver_signal_1_theta_1_by_h_array / pragmatic_receiver_signal_1_pre_normalization_factor_array = %s' % (pragmatic_receiver_signal_1_theta_1_by_h_array/pragmatic_receiver_signal_1_pre_normalization_factor_array)
				pragmatic_sender_signal_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_1_theta_1_by_h_array / pragmatic_receiver_signal_1_pre_normalization_factor_array) - cost))
				pragmatic_sender_signal_not_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_not_1_theta_1_by_h_array / pragmatic_receiver_signal_not_1_pre_normalization_factor_array) - (cost + cost_of_not)))
			else:
				pragmatic_sender_signal_0_non_normalized = numpy.ones(num_states)
				pragmatic_sender_signal_1_non_normalized = numpy.ones((num_states, num_states))
				pragmatic_sender_signal_not_1_non_normalized = numpy.ones((num_states, num_states))
			
			print 'pragmatic_sender_signal_0_non_normalized = \n%s' % pragmatic_sender_signal_0_non_normalized
			print 'pragmatic_sender_signal_1_non_normalized = \n%s' % pragmatic_sender_signal_1_non_normalized
			print 'pragmatic_sender_signal_not_1_non_normalized = \n%s' % pragmatic_sender_signal_not_1_non_normalized
		
			denominator_array = (numpy.tile(pragmatic_sender_signal_0_non_normalized, (num_states, 1)) + pragmatic_sender_signal_1_non_normalized + pragmatic_sender_signal_not_1_non_normalized)
		
			pragmatic_sender_signal_0_normalized = numpy.tile(pragmatic_sender_signal_0_non_normalized, (num_states, 1)) / denominator_array
			pragmatic_sender_signal_1_normalized = pragmatic_sender_signal_1_non_normalized / denominator_array
			pragmatic_sender_signal_not_1_normalized = pragmatic_sender_signal_not_1_non_normalized / denominator_array
		
			print 'level %s pragmatic_sender_signal_0_normalized = \n%s' % (n, pragmatic_sender_signal_0_normalized)
			print 'level %s pragmatic_sender_signal_1_normalized = \n%s' % (n, pragmatic_sender_signal_1_normalized)
			print 'level %s pragmatic_sender_signal_not_1_normalized = \n%s' % (n, pragmatic_sender_signal_not_1_normalized)

		elif pragmatic_sender_type == 'h_and_theta_sensitive':
		
			pragmatic_receiver_signal_0_pre_normalization_factor_array = numpy.ones(num_states)	* num_states		
			pragmatic_receiver_signal_1_pre_normalization_factor_array = numpy.ones((num_states, num_states))
			pragmatic_receiver_signal_not_1_pre_normalization_factor_array = numpy.ones((num_states, num_states))
			
			if choice_parameter != 0.:
				pragmatic_sender_signal_0_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_0_h_array / pragmatic_receiver_signal_0_pre_normalization_factor_array) - cost_of_null_signal))
				pragmatic_sender_signal_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_1_theta_1_by_h_array / pragmatic_receiver_signal_1_pre_normalization_factor_array) - cost))
				pragmatic_sender_signal_not_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(pragmatic_receiver_signal_not_1_theta_1_by_h_array / pragmatic_receiver_signal_not_1_pre_normalization_factor_array) - (cost + cost_of_not)))
			else:
				pragmatic_sender_signal_0_non_normalized = numpy.ones(num_states)
				pragmatic_sender_signal_1_non_normalized = numpy.ones((num_states, num_states))
				pragmatic_sender_signal_not_1_non_normalized = numpy.ones((num_states, num_states))

			print 'pragmatic_sender_signal_0_non_normalized = \n%s' % pragmatic_sender_signal_0_non_normalized
			print 'pragmatic_sender_signal_1_non_normalized = \n%s' % pragmatic_sender_signal_1_non_normalized
			print 'pragmatic_sender_signal_not_1_non_normalized = \n%s' % pragmatic_sender_signal_not_1_non_normalized
		
			denominator_array = (numpy.tile(pragmatic_sender_signal_0_non_normalized, (num_states, 1)) + pragmatic_sender_signal_1_non_normalized + pragmatic_sender_signal_not_1_non_normalized)
		
			pragmatic_sender_signal_0_normalized = numpy.tile(pragmatic_sender_signal_0_non_normalized, (num_states, 1)) / denominator_array
			pragmatic_sender_signal_1_normalized = pragmatic_sender_signal_1_non_normalized / denominator_array
			pragmatic_sender_signal_not_1_normalized = pragmatic_sender_signal_not_1_non_normalized / denominator_array
		
			print 'level %s pragmatic_sender_signal_0_normalized = \n%s' % (n, pragmatic_sender_signal_0_normalized)
			print 'level %s pragmatic_sender_signal_1_normalized = \n%s' % (n, pragmatic_sender_signal_1_normalized)
			print 'level %s pragmatic_sender_signal_not_1_normalized = \n%s' % (n, pragmatic_sender_signal_not_1_normalized)
		
		elif pragmatic_sender_type == 'modified h sensitive':
		
			pragmatic_receiver_signal_0_pre_normalization_factor_array = numpy.ones((num_states, num_states))
			pragmatic_receiver_signal_1_pre_normalization_factor_array = numpy.tile(numpy.sum(pragmatic_receiver(n - 1)[1], axis = (1)), (num_states, 1)).T
			pragmatic_receiver_signal_not_1_pre_normalization_factor_array = numpy.tile(numpy.sum(pragmatic_receiver(n - 1)[2], axis = (1)), (num_states, 1)).T
			
			print 'pragmatic_receiver_signal_0_h_array = \n%s' % pragmatic_receiver_signal_0_h_array
			print 'pragmatic_receiver_signal_1_theta_1_by_h_array = \n%s' % pragmatic_receiver_signal_1_theta_1_by_h_array
			print 'pragmatic_receiver_signal_not_1_theta_1_by_h_array = \n%s' % pragmatic_receiver_signal_not_1_theta_1_by_h_array
						
			print 'pragmatic_receiver_signal_0_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_0_pre_normalization_factor_array
			print 'pragmatic_receiver_signal_1_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_1_pre_normalization_factor_array
			print 'pragmatic_receiver_signal_not_1_pre_normalization_factor_array = \n%s' % pragmatic_receiver_signal_not_1_pre_normalization_factor_array
						
			pragmatic_sender_signal_0_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.sum(weighting_array.reshape(1, num_states, num_states) * (pragmatic_receiver_signal_0_h_array / pragmatic_receiver_signal_0_pre_normalization_factor_array).reshape(num_states, 1, num_states), axis = 2)) - cost_of_null_signal))
			pragmatic_sender_signal_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.sum(weighting_array.reshape(1, num_states, num_states) * (pragmatic_receiver_signal_1_theta_1_by_h_array / pragmatic_receiver_signal_1_pre_normalization_factor_array).reshape(num_states, 1, num_states), axis = 2)) - cost))
			pragmatic_sender_signal_not_1_non_normalized = numpy.exp(choice_parameter * (numpy.log(numpy.sum(weighting_array.reshape(1, num_states, num_states) * (pragmatic_receiver_signal_not_1_theta_1_by_h_array / pragmatic_receiver_signal_not_1_pre_normalization_factor_array).reshape(num_states, 1, num_states), axis = 2)) - (cost + cost_of_not)))
			
			print 'pragmatic_sender_signal_0_non_normalized = \n%s' % pragmatic_sender_signal_0_non_normalized
			print 'pragmatic_sender_signal_1_non_normalized = \n%s' % pragmatic_sender_signal_1_non_normalized
			print 'pragmatic_sender_signal_not_1_non_normalized = \n%s' % pragmatic_sender_signal_not_1_non_normalized
		
			denominator_array = (pragmatic_sender_signal_0_non_normalized + pragmatic_sender_signal_1_non_normalized + pragmatic_sender_signal_not_1_non_normalized)
		
			pragmatic_sender_signal_0_normalized = pragmatic_sender_signal_0_non_normalized / denominator_array
			pragmatic_sender_signal_1_normalized = pragmatic_sender_signal_1_non_normalized / denominator_array
			pragmatic_sender_signal_not_1_normalized = pragmatic_sender_signal_not_1_non_normalized / denominator_array
		
			print 'level %s pragmatic_sender_signal_0_normalized = \n%s' % (n, pragmatic_sender_signal_0_normalized)
			print 'level %s pragmatic_sender_signal_1_normalized = \n%s' % (n, pragmatic_sender_signal_1_normalized)
			print 'level %s pragmatic_sender_signal_not_1_normalized = \n%s' % (n, pragmatic_sender_signal_not_1_normalized)

		
		pragmatic_sender_memo[n] = numpy.asarray((pragmatic_sender_signal_0_normalized, pragmatic_sender_signal_1_normalized, pragmatic_sender_signal_not_1_normalized))
	return pragmatic_sender_memo[n]

def pragmatic_receiver(n):
	if n not in pragmatic_receiver_memo:
		if n==2:
		
			if theta_posterior_source == 'signal 1':
				theta_1_distribution_array_memo[0] = numpy.tile(numpy.sum(pragmatic_receiver(0)[1:2], axis = (0, 2)), (3, 1)).reshape(3, num_states, 1)
			elif theta_posterior_source == 'signal 1, not 1':
				theta_1_distribution_array_memo[0] = numpy.tile(numpy.sum(pragmatic_receiver(0)[1:], axis = (0, 2)) / 2., (3, 1)).reshape(3, num_states, 1)
			elif theta_posterior_source == 'signal 0, 1, not 1':
				theta_1_distribution_array_memo[0] = numpy.tile(numpy.sum(pragmatic_receiver(0)[0:], axis = (0, 2)) / 3., (3, 1)).reshape(3, num_states, 1)
			elif theta_posterior_source == 'level 0 prior':
				theta_1_distribution_array_memo[0] = numpy.tile(theta_1_distribution.pdf(array_0)/numpy.sum(theta_1_distribution.pdf(array_0)), (3, 1)).reshape(3, num_states, 1)
			elif theta_posterior_source == 'signal_specific':
				theta_1_distribution_array_memo[0] = numpy.sum(pragmatic_receiver(0), axis = 2).reshape(3, num_states, 1)
			elif theta_posterior_source == 'signal_weighted':
# 				here assume the level 0 receiver assumes all signals are equally likely
				theta_1_distribution_array_memo[0] = numpy.tile(numpy.sum(pragmatic_receiver(0)[0:], axis = (0, 2)) / 3., (3, 1)).reshape(3, num_states, 1)

			if h_posterior_source == 'signal_weighted':
				h_distribution_array_memo[0] = numpy.tile(numpy.sum(pragmatic_receiver(0)[0:], axis = (0, 1)) / 3., (3, 1)).reshape(3, 1, num_states)
# 				Here we preserve an option to have the signal weighted prior h distribution for level 2
# 				interpretation be the level 0 prior h distribution, thus making the signal weighted
# 				hierarchy a consistent extension of the original model. 
# 				h_distribution_array_memo[n-2] = state_distribution.pdf(array_0)/numpy.sum(state_distribution.pdf(array_0))
			elif h_posterior_source == 'signal_specific level n':
				h_distribution_array_memo[0] = numpy.tile(numpy.sum(pragmatic_receiver(0)[0], axis = (0)), (3, 1)).reshape(3, 1, num_states)
			elif h_posterior_source == 'signal 0 level n':
				h_distribution_array_memo[0] = numpy.tile(numpy.sum(pragmatic_receiver(0)[0], axis = (0)), (3, 1)).reshape(3, 1, num_states)
			elif h_posterior_source == 'level 0 prior':
				h_distribution_array_memo[0] = numpy.tile(state_distribution.pdf(array_0)/numpy.sum(state_distribution.pdf(array_0)), (3, 1)).reshape(3, 1, num_states)

		elif n>2:
		
			if theta_posterior_source == 'signal 1':
				theta_1_distribution_array_memo[n-2] = numpy.tile(numpy.sum(pragmatic_receiver(n-2)[1:2], axis = (0, 2)), (3, 1)).reshape(3, num_states, 1)
			elif theta_posterior_source == 'signal 1, not 1':
				theta_1_distribution_array_memo[n-2] = numpy.tile(numpy.sum(pragmatic_receiver(n-2)[1:], axis = (0, 2)) / 2., (3, 1)).reshape(3, num_states, 1)
			elif theta_posterior_source == 'signal 0, 1, not 1':
				theta_1_distribution_array_memo[n-2] = numpy.tile(numpy.sum(pragmatic_receiver(n-2)[0:], axis = (0, 2)) / 3., (3, 1)).reshape(3, num_states, 1)
			elif theta_posterior_source == 'level 0 prior':
				theta_1_distribution_array_memo[n-2] = numpy.tile(theta_1_distribution.pdf(array_0)/numpy.sum(theta_1_distribution.pdf(array_0)), (3, 1)).reshape(3, num_states, 1)
			elif theta_posterior_source == 'signal_specific':
				theta_1_distribution_array_memo[n-2] = numpy.sum(pragmatic_receiver(n-2), axis = 2).reshape(3, num_states, 1)
			elif theta_posterior_source == 'signal_weighted':
#				note here we have to call pragmatic_receiver(n-2) first, otherwise
#				theta_1_distribution_array_memo[n-4] will be undefined for n=6. otherwise
#				if we set max_level = 10, then pragmatic_sender(7) will be called, then
#				pragmatic_receiver(6) will be called, and then
#				theta_1_distribution_array_memo[2] will be undefined.
				receiver_posteriors = pragmatic_receiver(n-2)
				utterance_probabilities = numpy.sum(pragmatic_sender(n-3) * theta_1_distribution_array_memo[n-4] * h_distribution_array_memo[n-4], axis = (1, 2))
				theta_1_utterance_filter = numpy.reshape(numpy.tile(utterance_probabilities, (num_states, 1)).T, (3, num_states, 1))
				theta_1_distribution_array_memo[n-2] = numpy.tile(numpy.sum(receiver_posteriors * theta_1_utterance_filter, axis = (0, 2)), (3, 1)).reshape(3, num_states, 1)

			print 'theta_1_distribution_array_memo = \n%s' % theta_1_distribution_array_memo
			for level in sorted(theta_1_distribution_array_memo):
				print 'numpy.sum(theta_1_distribution_array_memo[%s]) = %s' % (level, numpy.sum(theta_1_distribution_array_memo[level]))

			if h_posterior_source == 'signal_weighted':
				receiver_posteriors = pragmatic_receiver(n-2)
				utterance_probabilities = numpy.sum(pragmatic_sender(n-3) * theta_1_distribution_array_memo[n-4] * h_distribution_array_memo[n-4], axis = (1, 2))
				h_utterance_filter = numpy.reshape(numpy.tile(utterance_probabilities, (num_states, 1)).T, (3, 1, num_states))
				h_distribution_array_memo[n-2] = numpy.tile(numpy.sum(receiver_posteriors * h_utterance_filter, axis = (0, 1)), (3, 1)).reshape(3, 1, num_states)
			elif h_posterior_source == 'signal_specific level n':
				h_distribution_array_memo[n-2] = numpy.sum(pragmatic_receiver(n-2), axis = 1).reshape(3, 1, num_states)
			elif h_posterior_source == 'signal 0 level n':
				h_distribution_array_memo[n-2] = numpy.tile(numpy.sum(pragmatic_receiver(n-2)[0], axis = (0)), (3, 1)).reshape(3, 1, num_states)
			elif h_posterior_source == 'level 0 prior':
				h_distribution_array_memo[n-2] = numpy.tile(state_distribution.pdf(array_0)/numpy.sum(state_distribution.pdf(array_0)), (3, 1)).reshape(3, 1, num_states)

			print 'h_distribution_array_memo = \n%s' % h_distribution_array_memo
			for level in sorted(h_distribution_array_memo):
				print 'numpy.sum(h_distribution_array_memo[%s]) = %s' % (level, numpy.sum(h_distribution_array_memo[level]))

		pragmatic_receiver_signal_0_pre_normalized = pragmatic_sender(n-1)[0] * theta_1_distribution_array_memo[n-2][0,:,:] * h_distribution_array_memo[n-2][0,:,:]
		pragmatic_receiver_signal_1_pre_normalized = pragmatic_sender(n-1)[1] * theta_1_distribution_array_memo[n-2][1,:,:] * h_distribution_array_memo[n-2][1,:,:]
		pragmatic_receiver_signal_not_1_pre_normalized = pragmatic_sender(n-1)[2] * theta_1_distribution_array_memo[n-2][2,:,:] * h_distribution_array_memo[n-2][2,:,:]
		
		print 'level %s pragmatic_receiver_signal_0_pre_normalized = \n%s' % (n, pragmatic_receiver_signal_0_pre_normalized)
		print 'level %s pragmatic_receiver_signal_1_pre_normalized = \n%s' % (n, pragmatic_receiver_signal_1_pre_normalized)
		print 'level %s pragmatic_receiver_signal_not_1_pre_normalized = \n%s' % (n, pragmatic_receiver_signal_not_1_pre_normalized)
		
		pragmatic_receiver_signal_0_array = pragmatic_receiver_signal_0_pre_normalized / numpy.sum(pragmatic_receiver_signal_0_pre_normalized)
		pragmatic_receiver_signal_1_array = pragmatic_receiver_signal_1_pre_normalized / numpy.sum(pragmatic_receiver_signal_1_pre_normalized)
		pragmatic_receiver_signal_not_1_array = pragmatic_receiver_signal_not_1_pre_normalized / numpy.sum(pragmatic_receiver_signal_not_1_pre_normalized)
		
		pragmatic_receiver_memo[n] = numpy.asarray((pragmatic_receiver_signal_0_array, pragmatic_receiver_signal_1_array, pragmatic_receiver_signal_not_1_array))

		pragmatic_sender_signal_0_average_h_array = numpy.sum(pragmatic_sender(n-1)[0], axis = 0) / num_states
		pragmatic_sender_signal_1_average_h_array = numpy.sum(pragmatic_sender(n-1)[1], axis = 0) / num_states
		pragmatic_sender_signal_not_1_average_h_array = numpy.sum(pragmatic_sender(n-1)[2], axis = 0) / num_states

		pragmatic_sender_signal_0_weighted_h_array = numpy.sum(pragmatic_sender(n-1)[0] * theta_1_distribution_array_memo[n-2][0,:,:], axis = 0)
		pragmatic_sender_signal_1_weighted_h_array = numpy.sum(pragmatic_sender(n-1)[1] * theta_1_distribution_array_memo[n-2][1,:,:], axis = 0)
		pragmatic_sender_signal_not_1_weighted_h_array = numpy.sum(pragmatic_sender(n-1)[2] * theta_1_distribution_array_memo[n-2][2,:,:], axis = 0)

		print 'level %s-1 pragmatic_sender_signal_0_weighted_h_array = \n%s' % (n, pragmatic_sender_signal_0_weighted_h_array)
		print 'level %s-1 pragmatic_sender_signal_1_weighted_h_array = \n%s' % (n, pragmatic_sender_signal_1_weighted_h_array)
		print 'level %s-1 pragmatic_sender_signal_not_1_weighted_h_array = \n%s' % (n, pragmatic_sender_signal_not_1_weighted_h_array)
		
		print 'pragmatic_sender_signal_0_weighted_h_array + pragmatic_sender_signal_1_weighted_h_array + pragmatic_sender_signal_not_1_weighted_h_array = %s' % (pragmatic_sender_signal_0_weighted_h_array + pragmatic_sender_signal_1_weighted_h_array + pragmatic_sender_signal_not_1_weighted_h_array)
		
		# note that here we calculate a normalized sending probability: since we divide by the probability of the utterance we need a 'dummy' uniform probability for the height such that multiplying the n-1 sending probabilities times the n-2 theta probabilities times the dummy uniform h probabilities, comes out to the same volume as multiplying the n-1 sending probabilities times the n-2 theta probabilities times the n-2 h probabilities. Since that latter product just is the utterance probability, we can set the h probabilities to the utterance probability. 
		normalizing_uniform_signal_0_h_probability = numpy.sum(pragmatic_sender_signal_0_weighted_h_array * h_distribution_array_memo[n-2][0,:,:])
		normalizing_uniform_signal_1_h_probability = numpy.sum(pragmatic_sender_signal_1_weighted_h_array * h_distribution_array_memo[n-2][1,:,:])
		normalizing_uniform_signal_not_1_h_probability = numpy.sum(pragmatic_sender_signal_not_1_weighted_h_array * h_distribution_array_memo[n-2][2,:,:])

		print 'normalizing_uniform_signal_0_h_probability = %s' % normalizing_uniform_signal_0_h_probability
		print 'normalizing_uniform_signal_1_h_probability = %s' % normalizing_uniform_signal_1_h_probability
		print 'normalizing_uniform_signal_not_1_h_probability = %s' % normalizing_uniform_signal_not_1_h_probability

		pragmatic_sender_signal_0_normalized_weighted_h_array = (pragmatic_sender_signal_0_weighted_h_array * normalizing_uniform_signal_0_h_probability) / numpy.sum(pragmatic_sender_signal_0_weighted_h_array * h_distribution_array_memo[n-2][0,:,:])
		pragmatic_sender_signal_1_normalized_weighted_h_array = (pragmatic_sender_signal_1_weighted_h_array * normalizing_uniform_signal_1_h_probability) / numpy.sum(pragmatic_sender_signal_1_weighted_h_array * h_distribution_array_memo[n-2][1,:,:])
		pragmatic_sender_signal_not_1_normalized_weighted_h_array = (pragmatic_sender_signal_not_1_weighted_h_array * normalizing_uniform_signal_not_1_h_probability) / numpy.sum(pragmatic_sender_signal_not_1_weighted_h_array * h_distribution_array_memo[n-2][2,:,:])

		pragmatic_sender_signal_0_fixed_theta_1_h_array = pragmatic_sender(n-1)[0][fixed_theta_1_num]
		pragmatic_sender_signal_1_fixed_theta_1_h_array = pragmatic_sender(n-1)[1][fixed_theta_1_num]
		pragmatic_sender_signal_not_1_fixed_theta_1_h_array = pragmatic_sender(n-1)[2][fixed_theta_1_num]

		pragmatic_receiver_signal_0_h_array = numpy.sum(pragmatic_receiver_memo[n][0], axis = 0)
		pragmatic_receiver_signal_0_h_array_densities = pragmatic_receiver_signal_0_h_array / ((upper_bound - lower_bound)/num_states)
		print 'pragmatic_receiver_signal_0_h_array = \n%s' % pragmatic_receiver_signal_0_h_array
	
		pragmatic_receiver_signal_1_h_array = numpy.sum(pragmatic_receiver_memo[n][1], axis = 0)
		pragmatic_receiver_signal_1_h_array_densities = pragmatic_receiver_signal_1_h_array / ((upper_bound - lower_bound)/num_states)
		print 'pragmatic_receiver_signal_1_h_array = \n%s' % pragmatic_receiver_signal_1_h_array
		
		pragmatic_receiver_signal_not_1_h_array = numpy.sum(pragmatic_receiver_memo[n][2], axis = 0)
		pragmatic_receiver_signal_not_1_h_array_densities = pragmatic_receiver_signal_not_1_h_array / ((upper_bound - lower_bound)/num_states)
		print 'pragmatic_receiver_signal_not_1_h_array = \n%s' % pragmatic_receiver_signal_not_1_h_array
	
		pragmatic_receiver_signal_1_theta_1_array = numpy.sum(pragmatic_receiver_memo[n][1], axis = 1)
		pragmatic_receiver_signal_1_theta_1_array_densities = pragmatic_receiver_signal_1_theta_1_array / ((upper_bound - lower_bound)/num_states)
		print 'pragmatic_receiver_signal_1_theta_1_array = \n%s' % pragmatic_receiver_signal_1_theta_1_array
	
		pragmatic_receiver_signal_not_1_theta_1_array = numpy.sum(pragmatic_receiver_memo[n][2], axis = 1)
		pragmatic_receiver_signal_not_1_theta_1_array_densities = pragmatic_receiver_signal_not_1_theta_1_array / ((upper_bound - lower_bound)/num_states)
		print 'pragmatic_receiver_signal_not_1_theta_1_array = \n%s' % pragmatic_receiver_signal_not_1_theta_1_array
	
		for signal in numpy.arange(len(pragmatic_sender(n-1))):
# 			plot_theta_1_distribution_array = numpy.swapaxes(theta_1_distribution_array_memo[n-2].reshape(num_states, 3), 0, 1)[signal]
			plot_theta_1_distribution_array = theta_1_distribution_array_memo[n-2][signal, :, :]
			fig, ax = pyplot.subplots(1,1)
			pyplot.plot(array_0, plot_theta_1_distribution_array)
			fig.text(0, 0, 'theta_1_distribution_array_memo[%s-2][%s, :, :]' % (n, signal), fontsize = 10)
			pyplot.show()
			pyplot.close()
			
			fig = pyplot.figure()
			ax = fig.gca(projection = '3d')
			ax.set_xlim(-4., 4.)
			ax.set_ylim(-4., 4.)
			ax.set_xlabel(r'$h$')
			ax.set_ylabel(r'$\theta_{1}$')
			surface = ax.plot_surface(numpy.tile(array_0, (num_states, 1)), numpy.tile(array_0, (num_states, 1)).T, pragmatic_sender(n-1)[signal], cmap = cm.coolwarm, linewidth = 0.0, antialiased = True, rstride = 2, cstride = 2, shade = False)
			fig.text(0, 0, 'pragmatic_sender(%s-1)[%s]' % (n, signal), fontsize = 10)
			pyplot.show()
			pyplot.close()
			
			fig = pyplot.figure()
			ax = fig.gca(projection = '3d')
			ax.set_xlim(-4., 4.)
			ax.set_ylim(-4., 4.)
			ax.set_xlabel(r'$h$')
			ax.set_ylabel(r'$\theta_{1}$')
			surface = ax.plot_surface(numpy.tile(array_0, (num_states, 1)), numpy.tile(array_0, (num_states, 1)).T, pragmatic_receiver_memo[n][signal], cmap = cm.coolwarm, linewidth = 0.0, antialiased = True, rstride = 2, cstride = 2, shade = False)
			fig.text(0, 0, 'pragmatic_receiver_memo[%s][%s]' % (n, signal), fontsize = 10)
			pyplot.show()
			pyplot.close()

		fig, ax = pyplot.subplots(1, 1, figsize = (12,5))
	
		pyplot.subplot(1, 1, 1)
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_0_average_h_array, lw = 2, color = 'k')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_1_average_h_array, lw = 2, color = 'b')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_average_h_array, lw = 2, color = 'c')
	
		line = pyplot.plot(array_0, pragmatic_sender_signal_0_weighted_h_array, lw = 10, color = 'k', linestyle = '-.')
		line = pyplot.plot(array_0, pragmatic_sender_signal_1_weighted_h_array, lw = 10, color = 'b', linestyle = '-.')
		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_weighted_h_array, lw = 10, color = 'c', linestyle = '-.')

		line = pyplot.plot(array_0, pragmatic_sender_signal_0_normalized_weighted_h_array, lw = 6, color = 'k', linestyle = '-')
		line = pyplot.plot(array_0, pragmatic_sender_signal_1_normalized_weighted_h_array, lw = 6, color = 'b', linestyle = '-')
		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_normalized_weighted_h_array, lw = 6, color = 'c', linestyle = '-')

# 		line = pyplot.plot(array_0, pragmatic_sender_signal_0_fixed_theta_1_h_array, lw = 2, linestyle = '--', color = 'k')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_1_fixed_theta_1_h_array, lw = 2, linestyle = '--', color = 'b')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_fixed_theta_1_h_array, lw = 2, linestyle = '--', color = 'c')
# 	
# 		line = pyplot.plot(array_0, pragmatic_sender(n-1)[0][:,fixed_theta_1_num], lw = 5, linestyle = ':', color = 'k')
# 		line = pyplot.plot(array_0, pragmatic_sender(n-1)[1][:,fixed_theta_1_num], lw = 5, linestyle = ':', color = 'b')
# 		line = pyplot.plot(array_0, pragmatic_sender(n-1)[2][:,fixed_theta_1_num], lw = 5, linestyle = ':', color = 'c')
	
# 		pyplot.subplot(1, 1, 1)
# 		pyplot.legend([r'$\sigma_{%s}(u_{0}|h)$' % (n-1), r'$\sigma_{%s}(u_{1}|h)$' % (n-1), r'$\sigma_{%s}(\neg u_{1}|h)$' % (n-1), r'$\sigma_{%s}(u_{0}|h, \theta_{1} \approx %s)$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2)), r'$\sigma_{%s}(u_{1}|h, \theta_{1} \approx %s)$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2)), r'$\sigma_{%s}(\neg u_{1}|h, \theta_{1} \approx %s)$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2)), r'$\sigma_{%s}(u_{0}|h \approx %s, \theta_{1})$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2)), r'$\sigma_{%s}(u_{1}|h \approx %s, \theta_{1})$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2)), r'$\sigma_{%s}(\neg u_{1}|h \approx %s, \theta_{1})$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2))], loc = (0.15, 0.15), fontsize = 14)
	
		fig.text(0, 0, r'$Lassiter\ and\ Goodman\ Two\ Signals\ with\ Not\ Iterated\ Arbitrary\ Variable\ Cost\ Variable\ Choice\ Parameter$' + '\n', fontsize = 10)
	
		fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \mu = %s, \sigma = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), mu, sigma, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(.4, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \alpha = %s, \beta = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), alpha_parameter, beta_parameter, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(.4, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, Uniform distribution, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)

		fig, ax = pyplot.subplots(1, 2, figsize = (12,5))
	
		pyplot.subplot(1, 2, 1)
		line = pyplot.plot(array_0, pragmatic_sender_signal_0_average_h_array, lw = 2, color = 'k')
		line = pyplot.plot(array_0, pragmatic_sender_signal_1_average_h_array, lw = 2, color = 'b')
		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_average_h_array, lw = 2, color = 'c')
	
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_0_weighted_h_array, lw = 2, color = 'k', linestyle = '-.')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_1_weighted_h_array, lw = 2, color = 'b', linestyle = '-.')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_weighted_h_array, lw = 2, color = 'c', linestyle = '-.')

# 		line = pyplot.plot(array_0, pragmatic_sender_signal_0_weighted_h_array * h_distribution_array_memo[n-2][0,:,:], lw = 2, color = 'k', linestyle = '-.')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_1_weighted_h_array * h_distribution_array_memo[n-2][1,:,:], lw = 2, color = 'b', linestyle = '-.')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_weighted_h_array * h_distribution_array_memo[n-2][2,:,:], lw = 2, color = 'c', linestyle = '-.')

# 		line = pyplot.plot(array_0, pragmatic_sender_signal_0_normalized_weighted_h_array_densities, lw = 2, color = 'k', linestyle = '-.')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_1_normalized_weighted_h_array_densities, lw = 2, color = 'b', linestyle = '-.')
# 		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_normalized_weighted_h_array_densities, lw = 2, color = 'c', linestyle = '-.')

		line = pyplot.plot(array_0, pragmatic_sender_signal_0_fixed_theta_1_h_array, lw = 2, linestyle = '--', color = 'k')
		line = pyplot.plot(array_0, pragmatic_sender_signal_1_fixed_theta_1_h_array, lw = 2, linestyle = '--', color = 'b')
		line = pyplot.plot(array_0, pragmatic_sender_signal_not_1_fixed_theta_1_h_array, lw = 2, linestyle = '--', color = 'c')
	
		line = pyplot.plot(array_0, pragmatic_sender(n-1)[1][:,fixed_theta_1_num], lw = 5, linestyle = ':', color = 'b')
		line = pyplot.plot(array_0, pragmatic_sender(n-1)[2][:,fixed_theta_1_num], lw = 5, linestyle = ':', color = 'c')
	
		pyplot.subplot(1, 2, 2)
		line = pyplot.plot(array_0, pragmatic_receiver_signal_0_h_array_densities, lw = 2, color = 'k')
		line = pyplot.plot(array_0, pragmatic_receiver_signal_1_h_array_densities, lw = 2, color = 'b')
		line = pyplot.plot(array_0, pragmatic_receiver_signal_1_theta_1_array_densities, lw = 2, linestyle = '--', color = 'b')
		line = pyplot.plot(array_0, pragmatic_receiver_signal_not_1_h_array_densities, lw = 2, color = 'c')
		line = pyplot.plot(array_0, pragmatic_receiver_signal_not_1_theta_1_array_densities, lw = 2, linestyle = '--', color = 'c')
	
# 		pyplot.subplot(1, 2, 1)
# 		pyplot.legend([r'$\sigma_{%s}(u_{0}|h)$' % (n-1), r'$\sigma_{%s}(u_{1}|h)$' % (n-1), r'$\sigma_{%s}(\neg u_{1}|h)$' % (n-1), r'$\sigma_{%s}(u_{0}|h, \theta_{1} \approx %s)$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2)), r'$\sigma_{%s}(u_{1}|h, \theta_{1} \approx %s)$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2)), r'$\sigma_{%s}(\neg u_{1}|h, \theta_{1} \approx %s)$' % ((n-1), numpy.around(array_0[fixed_theta_1_num], decimals = 2))], loc = 0, fontsize = 14)
	
		pyplot.subplot(1, 2, 2)
		pyplot.legend([r'$\rho_{%s}(h|u_{0})$' % n, r'$\rho_{%s}(h|u_{1})$' % n, r'$\rho_{%s}(\theta_{1}|u_{1})$' % n, r'$\rho_{%s}(h|\neg u_{1})$' % n, r'$\rho_{%s}(\theta_{1}|\neg u_{1})$' % n, ], loc = 0, fontsize = 14)

		fig.text(0, 0, r'$Lassiter\ and\ Goodman\ Two\ Signals\ with\ Not\ Iterated\ Arbitrary\ Variable\ Cost\ Variable\ Choice\ Parameter$' + '\n', fontsize = 10)
	
		fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \mu = %s, \sigma = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), mu, sigma, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(.4, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \alpha = %s, \beta = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), alpha_parameter, beta_parameter, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(.4, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, Uniform distribution, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
	
		# pyplot.savefig('Lassiter and Goodman Model Two Signals with Not Iterated Arbitrary Normal Distribution.pdf')
		# pyplot.savefig('Lassiter and Goodman Model Two Signals with Not Iterated Arbitrary Beta Distribution.pdf')
		# pyplot.savefig('Lassiter and Goodman Model Two Signals with Not Iterated Arbitrary Uniform Distribution.pdf')
	
		pyplot.show()
		pyplot.close()

	return pragmatic_receiver_memo[n]

###################

# Here we have the settings for a level 0 receiver decoding probabilities, given a fixed
# theta. This forms the common basis for both Lassiter and Goodman's original model and
# our modified model.

cost_of_null_signal = 0.
cost_list = [1.]
cost_of_not_list = [1./3.]
choice_parameter_list = [4.]
lower_bound = -8.
upper_bound = 8.
num_states = 160

mu = 0.
sigma = 1.
state_distribution = norm(mu,sigma)

# alpha_parameter = 1.
# beta_parameter = 9.
# location_parameter = lower_bound
# scale_parameter = upper_bound - lower_bound
# state_distribution = beta(alpha_parameter, beta_parameter, loc = location_parameter, scale = scale_parameter)

# state_distribution = uniform(lower_bound, upper_bound - lower_bound)

fixed_theta_1_num = numpy.int(numpy.ceil(num_states*(6.75/12.)))
# fixed_theta_2_num = numpy.int(numpy.ceil(num_states*(4./12.)))

theta_distribution_type = 'uniform'
# theta_distribution_relation = 'True'
theta_posterior_source = 'signal_specific'
h_posterior_source = 'level 0 prior'

if theta_distribution_type == 'normal':
	theta_1_distribution = norm(mu+1.5, 1. * sigma)
# 	theta_2_distribution = norm(mu, sigma)
elif theta_distribution_type == 'Beta':
	theta_1_distribution = beta(1, 9, loc = lower_bound, scale = upper_bound - lower_bound)
# 	theta_2_distribution = beta(1, 9, loc = lower_bound, scale = upper_bound - lower_bound)
elif theta_distribution_type == 'uniform':
	theta_1_distribution = uniform(lower_bound, upper_bound - lower_bound)
# 	theta_2_distribution = uniform(lower_bound, upper_bound - lower_bound)

array_0 = numpy.flipud(numpy.linspace(upper_bound, lower_bound, num_states, endpoint = False)) - ((numpy.flipud(numpy.linspace(upper_bound, lower_bound, num_states, endpoint = False)) - numpy.linspace(lower_bound, upper_bound, num_states, endpoint = False))/2)
print 'array_0 = %s' % array_0

pragmatic_sender_type = 'h_sensitive'
# pragmatic_sender_type = 'h_and_theta_sensitive'
# pragmatic_sender_type = 'h_sensitive_version_3'
# pragmatic_sender_type = 'modified h sensitive'

if pragmatic_sender_type == 'modified h sensitive':
	weighting_sigma = 1.

if pragmatic_sender_type == 'modified h sensitive':
	weighting_array = numpy.empty((0, num_states))
	for h_num in numpy.arange(num_states):
		weighting_array = numpy.insert(weighting_array, h_num, truncnorm.pdf(array_0, lower_bound, upper_bound, loc = array_0[h_num], scale = weighting_sigma), axis = 0)
	print 'weighting_array = %s' % weighting_array
	
theta_1_distribution_array_memo = {}
h_distribution_array_memo = {}

max_level = 6

#########################

theta_1_distribution_array = theta_1_distribution.pdf(array_0)
theta_1_distribution_array = theta_1_distribution_array/numpy.sum(theta_1_distribution_array)
initial_theta_1_distribution_array = theta_1_distribution_array
initial_theta_1_distribution_array = numpy.reshape(initial_theta_1_distribution_array, [num_states, 1])

print 'initial_theta_1_distribution_array = \n%s' % initial_theta_1_distribution_array
print numpy.sum(initial_theta_1_distribution_array)

plot_theta_1_distribution_array = initial_theta_1_distribution_array
# plot_theta_2_distribution_array = numpy.sum(initial_theta_2_by_theta_1_distribution_array, axis = 1)

fig, ax = pyplot.subplots(1,1)
pyplot.plot(array_0, plot_theta_1_distribution_array)
# pyplot.plot(array_0, plot_theta_2_distribution_array)
pyplot.show()


#########################

receiver_0_signal_0_array = state_distribution.pdf(array_0)
receiver_0_signal_0_array = receiver_0_signal_0_array / numpy.sum(receiver_0_signal_0_array)
receiver_0_signal_0_array = numpy.tile(receiver_0_signal_0_array, (num_states, 1))
receiver_0_signal_0_array = receiver_0_signal_0_array * initial_theta_1_distribution_array

# receiver_0_signal_1_array = numpy.empty([0, num_states])
# for theta_num in range(num_states):
# 	temp_signal_1_fixed_theta_array = numpy.empty(0)
# 	for h_num in range(num_states):
# 		value = receiver_0_signal_1(array_0[h_num], array_0[theta_num])
# 		temp_signal_1_fixed_theta_array = numpy.append(temp_signal_1_fixed_theta_array, value)
# 	receiver_0_signal_1_array = numpy.insert(receiver_0_signal_1_array, theta_num, temp_signal_1_fixed_theta_array, axis = 0)
# receiver_0_signal_1_array = receiver_0_signal_1_array / numpy.tile(numpy.sum(receiver_0_signal_1_array, axis = 1), (num_states, 1)).T
# # receiver_0_signal_1_array = numpy.tile(receiver_0_signal_1_array, [num_states, 1, 1])
# receiver_0_signal_1_array = receiver_0_signal_1_array * initial_theta_1_distribution_array


state_distribution_pdf_array = state_distribution.pdf(array_0)
receiver_0_signal_1_array = numpy.empty([0, num_states])
for theta_num in range(num_states):
	temp_signal_1_fixed_theta_array = numpy.empty(0)
	temp_signal_1_renormalization_sum = numpy.sum(state_distribution_pdf_array[theta_num:])
	for h_num in range(num_states):
		if h_num < theta_num:
			temp_signal_1_fixed_theta_array = numpy.append(temp_signal_1_fixed_theta_array, 0.)	
		else:
			temp_signal_1_fixed_theta_array = numpy.append(temp_signal_1_fixed_theta_array, state_distribution_pdf_array[h_num]/temp_signal_1_renormalization_sum)	
	receiver_0_signal_1_array = numpy.insert(receiver_0_signal_1_array, theta_num, temp_signal_1_fixed_theta_array, axis = 0)
receiver_0_signal_1_array = receiver_0_signal_1_array / numpy.tile(numpy.sum(receiver_0_signal_1_array, axis = 1), (num_states, 1)).T
# receiver_0_signal_1_array = numpy.tile(receiver_0_signal_1_array, [num_states, 1, 1])
receiver_0_signal_1_array = receiver_0_signal_1_array * initial_theta_1_distribution_array



receiver_0_signal_not_1_array = numpy.empty([0, num_states])
for theta_num in range(num_states):
	temp_signal_not_1_fixed_theta_array = numpy.empty(0)
	for h_num in range(num_states):
		value = receiver_0_signal_not_1(array_0[h_num], array_0[theta_num])
		temp_signal_not_1_fixed_theta_array = numpy.append(temp_signal_not_1_fixed_theta_array, value)
	receiver_0_signal_not_1_array = numpy.insert(receiver_0_signal_not_1_array, theta_num, temp_signal_not_1_fixed_theta_array, axis = 0)
receiver_0_signal_not_1_array = receiver_0_signal_not_1_array / numpy.tile(numpy.sum(receiver_0_signal_not_1_array, axis = 1), (num_states, 1)).T
# receiver_0_signal_2_array = numpy.tile(numpy.reshape(receiver_0_signal_2_array, [num_states, 1, num_states]), [1, num_states, 1])
receiver_0_signal_not_1_array = receiver_0_signal_not_1_array * initial_theta_1_distribution_array

for choice_parameter_num in numpy.arange(len(choice_parameter_list)):
	choice_parameter = choice_parameter_list[choice_parameter_num]
	for cost_num in numpy.arange(len(cost_list)):
		cost = cost_list[cost_num]
		cost_of_not = cost_of_not_list[cost_num]
		pragmatic_receiver_memo = {}
		pragmatic_sender_memo = {}
		pragmatic_receiver_memo[0] = numpy.asarray((receiver_0_signal_0_array, receiver_0_signal_1_array, receiver_0_signal_not_1_array))
		for signal in numpy.arange(len(pragmatic_receiver(0))):
			fig = pyplot.figure()
			ax = fig.gca(projection = '3d')
			ax.set_xlim(-4., 4.)
			ax.set_ylim(-4., 4.)
			ax.set_xlabel(r'$h$')
			ax.set_ylabel(r'$\theta_{1}$')
			surface = ax.plot_surface(numpy.tile(array_0, (num_states, 1)), numpy.tile(array_0, (num_states, 1)).T, pragmatic_receiver(0)[signal]/initial_theta_1_distribution_array, cmap = cm.coolwarm, linewidth = 0.0, antialiased = True, rstride = 2, cstride = 2, shade = False)
			pyplot.show()
			pyplot.close()
		pragmatic_receiver(max_level)

choice_parameters_costs_pragmatic_receiver_h_array_densities = numpy.empty([0, len(cost_list), max_level/2 + 1, 3, num_states])
for choice_parameter_num in numpy.arange(len(choice_parameter_list)):
	costs_pragmatic_receiver_h_array_densities = numpy.empty([0, max_level/2 + 1, 3, num_states])
	for cost_num in numpy.arange(len(cost_list)):
		levels_pragmatic_receiver_h_array_densities = numpy.empty([0, 3, num_states])
		for pragmatic_receiver_level in sorted(pragmatic_receiver_memo):
			levels_pragmatic_receiver_h_array_densities = numpy.insert(levels_pragmatic_receiver_h_array_densities, pragmatic_receiver_level/2, numpy.sum(pragmatic_receiver_memo[pragmatic_receiver_level], axis = 1)/((upper_bound - lower_bound)/num_states), axis =0)
		costs_pragmatic_receiver_h_array_densities = numpy.insert(costs_pragmatic_receiver_h_array_densities, cost_num, levels_pragmatic_receiver_h_array_densities, axis = 0)
	choice_parameters_costs_pragmatic_receiver_h_array_densities = numpy.insert(choice_parameters_costs_pragmatic_receiver_h_array_densities, choice_parameter_num, costs_pragmatic_receiver_h_array_densities, axis = 0)

choice_parameters_costs_pragmatic_receiver_theta_1_array_densities = numpy.empty([0, len(cost_list), max_level/2 + 1, 3, num_states])
for choice_parameter_num in numpy.arange(len(choice_parameter_list)):
	costs_pragmatic_receiver_theta_1_array_densities = numpy.empty([0, max_level/2 + 1, 3, num_states])
	for cost_num in numpy.arange(len(cost_list)):
		levels_pragmatic_receiver_theta_1_array_densities = numpy.empty([0, 3, num_states])
		for pragmatic_receiver_level in sorted(pragmatic_receiver_memo):
			levels_pragmatic_receiver_theta_1_array_densities = numpy.insert(levels_pragmatic_receiver_theta_1_array_densities, pragmatic_receiver_level/2, numpy.sum(pragmatic_receiver_memo[pragmatic_receiver_level], axis = 2)/((upper_bound - lower_bound)/num_states), axis =0)
		costs_pragmatic_receiver_theta_1_array_densities = numpy.insert(costs_pragmatic_receiver_theta_1_array_densities, cost_num, levels_pragmatic_receiver_theta_1_array_densities, axis = 0)
	choice_parameters_costs_pragmatic_receiver_theta_1_array_densities = numpy.insert(choice_parameters_costs_pragmatic_receiver_theta_1_array_densities, choice_parameter_num, costs_pragmatic_receiver_theta_1_array_densities, axis = 0)

choice_parameters_costs_pragmatic_sender_h_array_densities = numpy.empty([0, len(cost_list), max_level/2, 3, num_states])
for choice_parameter_num in numpy.arange(len(choice_parameter_list)):
	costs_pragmatic_sender_h_array_densities = numpy.empty([0, max_level/2, 3, num_states])
	for cost_num in numpy.arange(len(cost_list)):
		levels_pragmatic_sender_h_array_densities = numpy.empty([0, 3, num_states])
		for pragmatic_sender_level in sorted(pragmatic_sender_memo):
			levels_pragmatic_sender_h_array_densities = numpy.insert(levels_pragmatic_sender_h_array_densities, (pragmatic_sender_level/2), numpy.sum(pragmatic_sender(pragmatic_sender_level) * theta_1_distribution_array_memo[pragmatic_sender_level - 1], axis = 1), axis =0)
		costs_pragmatic_sender_h_array_densities = numpy.insert(costs_pragmatic_sender_h_array_densities, cost_num, levels_pragmatic_sender_h_array_densities, axis = 0)
	choice_parameters_costs_pragmatic_sender_h_array_densities = numpy.insert(choice_parameters_costs_pragmatic_sender_h_array_densities, choice_parameter_num, costs_pragmatic_sender_h_array_densities, axis = 0)

# choice_parameters_costs_pragmatic_sender_theta_1_array_densities = numpy.empty([0, len(cost_list), max_level/2 + 1, 3, num_states])
# for choice_parameter_num in numpy.arange(len(choice_parameter_list)):
# 	costs_pragmatic_sender_theta_1_array_densities = numpy.empty([0, max_level/2 + 1, 3, num_states])
# 	for cost_num in numpy.arange(len(cost_list)):
# 		levels_pragmatic_sender_theta_1_array_densities = numpy.empty([0, 3, num_states])
# 		for pragmatic_sender_level in sorted(pragmatic_sender_memo):
# 			levels_pragmatic_sender_theta_1_array_densities = numpy.insert(levels_pragmatic_sender_theta_1_array_densities, (pragmatic_sender_level/2), numpy.sum(pragmatic_sender(pragmatic_sender_level) * theta_1_distribution_array_memo[pragmatic_sender_level - 1], axis = 2)/((upper_bound - lower_bound)/num_states), axis =0)
# 		costs_pragmatic_sender_theta_1_array_densities = numpy.insert(costs_pragmatic_sender_theta_1_array_densities, cost_num, levels_pragmatic_sender_theta_1_array_densities, axis = 0)
# 	choice_parameters_costs_pragmatic_sender_theta_1_array_densities = numpy.insert(choice_parameters_costs_pragmatic_sender_theta_1_array_densities, choice_parameter_num, costs_pragmatic_sender_theta_1_array_densities, axis = 0)

for pragmatic_receiver_level in sorted(pragmatic_receiver_memo):

	fig, ax = pyplot.subplots(1, 1, figsize = (12,5))
	color_list = ['k', 'b', 'r']
	linestyle_list = ['-', '--']
	lw_list = [1., 3.]
	pyplot.subplot(1, 1, 1)
	pyplot.grid(True)

	for choice_parameter_num in numpy.arange(len(choice_parameter_list)):
		for cost_num in numpy.arange(len(cost_list)):
			for signal in numpy.arange(3):
				pyplot.plot(array_0, choice_parameters_costs_pragmatic_receiver_h_array_densities[choice_parameter_num, cost_num, pragmatic_receiver_level/2, signal], color = color_list[signal], linestyle = linestyle_list[cost_num], lw = lw_list[choice_parameter_num])

	real_legend = pyplot.legend([r'$\rho_{%s}(h|u_{0})$' % pragmatic_receiver_level, r'$\rho_{%s}(h|u_{1})$' % pragmatic_receiver_level, r'$\rho_{%s}(h|\neg u_{1})$' % pragmatic_receiver_level], loc = 'lower left', bbox_to_anchor = (.5, .5))
	for legobj in real_legend.legendHandles:
		legobj.set_linewidth(1.5)

	ax = pyplot.gca().add_artist(real_legend)
	dummy_line_0, = pyplot.plot([], label = r'$C(u) = 1/3 * length(u)$', color = 'gray', lw = 2)
	dummy_line_1, = pyplot.plot([], label = r'$C(u) = 4/3 * length(u)$', color = 'gray', lw = 2, linestyle = '--')
	dummy_line_2, = pyplot.plot([], label = r'$\lambda = 4.$', color = 'gray', lw = 1)
	dummy_line_3, = pyplot.plot([], label = r'$\lambda = 8.$', color = 'gray', lw = 3)
	pyplot.legend(handles = [dummy_line_0, dummy_line_1, dummy_line_2, dummy_line_3], loc = 'lower right', bbox_to_anchor = (.5, .5), fontsize = 12)

	fig.text(0, 0, r'$Lassiter\ and\ Goodman\ Two\ Signals\ with\ Not\ Iterated\ Arbitrary\ Variable\ Cost\ Variable\ Choice\ Parameter$' + '\n', fontsize = 10)

	fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \mu = %s, \sigma = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter_list, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost_list, decimals = 2)), str(numpy.around(cost_of_not_list, decimals = 2)), mu, sigma, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 	fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \alpha = %s, \beta = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost_list, decimals = 2)), str(numpy.around(cost_of_not_list, decimals = 2)), alpha_parameter, beta_parameter, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 	fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, Uniform distribution, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost_list, decimals = 2)), str(numpy.around(cost_of_not_list, decimals = 2)), num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)

	pyplot.show()
	pyplot.close()

for choice_parameter_num in numpy.arange(len(choice_parameter_list)):
	for cost_num in numpy.arange(len(cost_list)):
		
		fig = pyplot.figure(figsize = (12,5))
		color_list = ['k', 'b', 'r']
		linestyle_list = ['-', '--']
		lw_list = [1., 3.]
		ax = fig.add_subplot(111)
		pyplot.grid(True)
		ax.set_ylim(0., 1.2)
	
		for pragmatic_receiver_level in sorted(pragmatic_receiver_memo):
			for signal in numpy.arange(3):
				pyplot.plot(array_0, choice_parameters_costs_pragmatic_receiver_h_array_densities[choice_parameter_num, cost_num, pragmatic_receiver_level/2, signal], color = color_list[signal], linestyle = linestyle_list[cost_num], lw = pragmatic_receiver_level/2 + 1)
# 				pyplot.plot(array_0, choice_parameters_costs_pragmatic_receiver_theta_1_array_densities[choice_parameter_num, cost_num, pragmatic_receiver_level/2, signal], color = color_list[signal], linestyle = linestyle_list[1], lw = pragmatic_receiver_level/2 + 1)

		ax.text(.5, .90, r'$\lambda = %s, C(u_{0}) \approx %s,$' % (choice_parameter_list[choice_parameter_num], str(numpy.around(cost_of_null_signal, decimals = 2))) + '\n' + r'$C(u_{n}) \approx %s, C(\neg u_{n}) \approx %s$' % (str(numpy.around(cost_list[cost_num], decimals = 2)), str(numpy.around(cost_list[cost_num] + cost_of_not_list[cost_num], decimals = 2))), bbox={'facecolor':'white', 'alpha':1., 'pad':10}, horizontalalignment = 'center', verticalalignment = 'center', transform = ax.transAxes)

		real_legend = pyplot.legend([r'$\rho_{n}(h|u_{0})$', r'$\rho_{n}(h|u_{1})$', r'$\rho_{n}(h|\neg u_{1})$'], loc = 'upper left', bbox_to_anchor = (.025, .975))
		for legobj in real_legend.legendHandles:
			legobj.set_linewidth(1.5)
		ax = pyplot.gca().add_artist(real_legend)

		dummy_line_list = []
		for pragmatic_receiver_level in sorted(pragmatic_receiver_memo):
			dummy_line, = pyplot.plot([], label = r'$\rho_{%s}$' % pragmatic_receiver_level, color = 'gray', lw = pragmatic_receiver_level/2 + 1)
			dummy_line_list.append(dummy_line)
		pyplot.legend(handles = dummy_line_list, loc = 'upper right', bbox_to_anchor = (.975, .975), fontsize = 14)
		
		fig.text(0, 0, r'$Lassiter\ and\ Goodman\ Two\ Signals\ with\ Not\ Iterated\ Arbitrary\ Variable\ Cost\ Variable\ Choice\ Parameter$' + '\n', fontsize = 10)

		fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \mu = %s, \sigma = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter_list[choice_parameter_num], str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost_list[cost_num], decimals = 2)), str(numpy.around(cost_of_not_list[cost_num], decimals = 2)), mu, sigma, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \alpha = %s, \beta = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost_list, decimals = 2)), str(numpy.around(cost_of_not_list, decimals = 2)), alpha_parameter, beta_parameter, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, Uniform distribution, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost_list, decimals = 2)), str(numpy.around(cost_of_not_list, decimals = 2)), num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)

		pyplot.show()
		pyplot.close()

for choice_parameter_num in numpy.arange(len(choice_parameter_list)):
	for cost_num in numpy.arange(len(cost_list)):
		
		fig = pyplot.figure(figsize = (12,5))
		color_list = ['k', 'b', 'r']
		linestyle_list = ['-', '--']
		lw_list = [1., 3.]
		ax = fig.add_subplot(111)
		pyplot.grid(True)
		ax.set_ylim(0., 1.0)
	
		for pragmatic_sender_level in sorted(pragmatic_sender_memo):
			for signal in numpy.arange(3):
				pyplot.plot(array_0, choice_parameters_costs_pragmatic_sender_h_array_densities[choice_parameter_num, cost_num, pragmatic_sender_level/2, signal], color = color_list[signal], linestyle = linestyle_list[cost_num], lw = pragmatic_sender_level/4. + 1)

		ax.text(.5, .70, r'$\lambda = %s, C(u_{0}) \approx %s,$' % (choice_parameter_list[choice_parameter_num], str(numpy.around(cost_of_null_signal, decimals = 2))) + '\n' + r'$C(u_{n}) \approx %s, C(\neg u_{n}) \approx %s$' % (str(numpy.around(cost_list[cost_num], decimals = 2)), str(numpy.around(cost_list[cost_num] + cost_of_not_list[cost_num], decimals = 2))), bbox={'facecolor':'white', 'alpha':1., 'pad':10}, horizontalalignment = 'center', verticalalignment = 'center', transform = ax.transAxes)

		real_legend = pyplot.legend([r'$\sigma_{n}(u_{0}|h)$', r'$\sigma_{n}(u_{1}|h)$', r'$\sigma_{n}(\neg u_{1}|h)$'], loc = 'upper left', bbox_to_anchor = (.025, .975))
		for legobj in real_legend.legendHandles:
			legobj.set_linewidth(1.5)
		ax = pyplot.gca().add_artist(real_legend)

		dummy_line_list = []
		for pragmatic_sender_level in sorted(pragmatic_sender_memo):
			dummy_line, = pyplot.plot([], label = r'$\sigma_{%s}$' % pragmatic_sender_level, color = 'gray', lw = pragmatic_sender_level/4. + 1)
			dummy_line_list.append(dummy_line)
		pyplot.legend(handles = dummy_line_list, loc = 'upper right', bbox_to_anchor = (.975, .975), fontsize = 14)
	
		fig.text(0, 0, r'$Lassiter\ and\ Goodman\ Two\ Signals\ with\ Not\ Iterated\ Arbitrary\ Variable\ Cost\ Variable\ Choice\ Parameter$' + '\n', fontsize = 10)
	
		fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \mu = %s, \sigma = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), mu, sigma, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, \alpha = %s, \beta = %s, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), alpha_parameter, beta_parameter, num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)
# 		fig.text(0, 0, r'$\lambda = %s, C(u_{0}) \approx %s, C(u_{n}) \approx %s, C(\neg) \approx %s, Uniform distribution, num\ states = %s, theta\ distribution\ type = %s,$' % (choice_parameter, str(numpy.around(cost_of_null_signal, decimals = 2)), str(numpy.around(cost, decimals = 2)), str(numpy.around(cost_of_not, decimals = 2)), num_states, theta_distribution_type) + r'$theta\ posterior\ source = %s, h\ posterior\ source = %s, pragmatic\ sender\ type = %s$' % (theta_posterior_source, h_posterior_source, pragmatic_sender_type), fontsize = 10)

		pyplot.show()
		pyplot.close()