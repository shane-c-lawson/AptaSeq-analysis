import scipy
from scipy.stats import kurtosis, skew
import gzip
import operator
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import islice
import sys
from scipy.spatial.distance import hamming

pos_control_1 = 'ATCCTCGTCCCGTCACGGCAGAACCACGTCAGGCCTTCAA'
pos_control_2 = 'CCCTAGTTACTACTACTCTTTTTAGCAAACGCCCTCGCTT' #full-length = TGCCCTAGTTACTACTACTCTTTTTAGCAAACGCCCTCGCTT
neg_control_1 = 'AGTCCATTTTATTCCTGAATATTTGTTAACCTCATGGACN' #requires equal lengths
neg_control_2 = 'TGGACAGGGTTAGGCGTAGGAGTTGAGTTTTTGAGACANN' #requires equal lengths
neg_control_3 = 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'

def Rev_comp(seq) :
	comp = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'N' :'N'}
	return ''.join([comp[x] for x in seq])[::-1]

with open('out.txt', 'r') as out :
	for line1 in out :
		if 'reads from' in line1 :
			input = line1.split(' ')[1]
		elif 'Of those' in line1 :
			passing = line1.split(' ')[2]
		elif 'Printed' in line1 :
			deduped = line1.split(' ')[1]

base_dir = "/home/FCAM/slawson/miseq/210610_M01212_0292_000000000-JPPMK/Data/Intensities/BaseCalls/AptaSeq/"
sample = sys.argv[1]; sample2 = sample.split('_')[0]
read_fields = ['read', 'seq', 'strand', 'qual']
openfile = open(base_dir + 'deduped/' + sample + '_R2_processed.fastq', 'r')
i = 0; j = 0; k = 0; l = 0; m = 0; dups = 0; read_info = {}; aptamers = {}; groups = {}; control_hits = []; control_hits_bad = []; c = 0; conts = 0; control_found = False
for line in openfile :
	line = line.strip()
	if i == 4 :
		j += 1; i = 0
		if read_info['seq'][16:22] == 'TAAAGC' :
			aptamer = Rev_comp(read_info['seq'][22:62])
			if 'N' in aptamer :
				if i == 0 :
					read_info['read'] = line
				if i == 1 :
					read_info['seq'] = line
				if i == 2 :
					read_info['strand'] = line
				if i == 3 :
					read_info['qual'] = line
				i = i + 1
				continue
			k += 1
			if aptamer in aptamers :
				dups += 1
				aptamers[aptamer] += 1
			else :
				aptamers[aptamer] = 1
				if aptamer == pos_control_1 :	#finds literal match instead of % match
					l += 1
					control_hits.append(aptamer)
					control_found = True
				if aptamer == pos_control_2 :
					l += 1
					control_hits.append(aptamer)
					control_found = True
				if aptamer == neg_control_1 :
					m += 1
					control_hits_bad.append(aptamer)
					control_found = True
				if aptamer == neg_control_2 :
					m += 1
					control_hits_bad.append(aptamer)
					control_found = True
				if aptamer == neg_control_3 :
					m += 1
					control_hits_bad.append(aptamer)
					control_found = True
				#dist = hamming(list(pos_control_1),list(aptamer))
				#if dist <= 0.1 :	#using 4 mismatches, or >= 90% identity, as the cutoff
				#	l += 1
				#	control_hits.append(aptamer)
				#	control_found = True
				#dist = hamming(list(pos_control_2),list(aptamer))
				#if dist <= 0.1 :
				#	l += 1
				#	control_hits.append(aptamer)
				#	control_found = True
				#dist = hamming(list(neg_control_1),list(aptamer))
				#if dist <= 0.1 :
				#	m += 1
				#	control_hits_bad.append(aptamer)
				#	control_found = True
				#dist = hamming(list(neg_control_2),list(aptamer))
				#if dist <= 0.1 :
				#	m += 1
				#	control_hits_bad.append(aptamer)
				#	control_found = True
				#dist = hamming(list(neg_control_3),list(aptamer))
				#if dist <= 0.1 :
				#	m += 1
				#	control_hits_bad.append(aptamer)
				#	control_found = True
	if i == 0 :
		read_info['read'] = line
	if i == 1 :
		read_info['seq'] = line
	if i == 2 :
		read_info['strand'] = line
	if i == 3 :
		read_info['qual'] = line
	i = i + 1
openfile.close()

aptamers_sorted = dict(sorted(aptamers.items(), key = operator.itemgetter(1), reverse = True))

if control_found == True :
	colors = []
	for aptamer in aptamers_sorted.keys() :
		if aptamer in control_hits :
			conts += aptamers_sorted[aptamer]
			colors.append('g')
		elif aptamer in control_hits_bad :
			colors.append('r')
		else:
			colors.append('b')
	plt.bar(range(200), list(aptamers_sorted.values())[0:200], color = colors)
	#plt.bar(range(len(aptamers_sorted)), list(aptamers_sorted.values()), color = colors)
else :
	plt.bar(range(200), list(aptamers_sorted.values())[0:200])
	#plt.bar(range(len(aptamers_sorted)), list(aptamers_sorted.values()))
#plt.xticks(range(len(aptamers_sorted)), list(aptamers_sorted.keys()), rotation=90)
plt.xticks([])
plt.title('Top 200 Aptamer Distribution - ' + sample2)
plt.xlabel('Aptamers')
plt.ylabel('Reads')
plt.savefig(base_dir + 'figures/' + sample2 + '_aptamer_distribution.png')
plt.clf()
data = list(aptamers_sorted.values())
apt_array = np.array(data).astype(np.float)
skew = skew(apt_array)

for aptamer in aptamers_sorted :
	if aptamers_sorted[aptamer] > 1 : # and c < 1 :
		groups[aptamer] = aptamers_sorted[aptamer]
		#c += 1
	#else :		#include for complete hamming check
	#	matched = False
	#	for group in groups.keys() :
	#		dist = hamming(list(group),list(aptamer))
	#		if dist <= 0.1 :
	#			groups[group] += 1
	#			matched = True
	#			break
	#	if matched == False and c < 11:
	#		groups[aptamer] = 1
	#		c += 1

if len(groups) > 0 :
	groups_sorted = dict(sorted(groups.items(), key = operator.itemgetter(1), reverse = True))
	groups_sorted_perc = {}
	for x, y in groups_sorted.items() :
		groups_sorted_perc[x] = round(float(y/k)*100, 10)
	groups_perc_sorted = dict(sorted(groups_sorted_perc.items(), key = operator.itemgetter(1), reverse = True))

	if control_found == True :
		colors = []
		for aptamer in groups_perc_sorted.keys() :
			if aptamer in control_hits :
				colors.append('g')
			elif aptamer in control_hits_bad :
				colors.append('r')
			else:
				colors.append('b')
		plt.bar(range(len(groups_perc_sorted)), list(groups_perc_sorted.values()), color = colors)
	else :
		plt.bar(range(len(groups_perc_sorted)), list(groups_perc_sorted.values()))
	plt.xticks(range(len(groups_perc_sorted)), list(groups_perc_sorted.keys()), rotation = 90, fontsize=3)
	plt.title('Aptamers with Duplicate Reads - ' + sample2)
	plt.xlabel('Aptamers')
	plt.ylabel('Percent of Reads')
	low = min(list(groups_perc_sorted.values()))
	high = max(list(groups_perc_sorted.values()))
	plt.ylim([0, (high+0.5*(high))])
	plt.tight_layout()
	plt.savefig(base_dir + 'figures/' + sample2 + '_best_aptamers.png')
	
print(sample2)
print('Total reads: ' + input)
print('Reads with passing UMI: ' + passing)
print('Total deduplicated reads: ' + deduped)	#this should be the same as j
#print('Total deduplicated reads: ' + str(j))
print('Usable reads: ' + str(k))
print('Percent usable reads: ' + str(round(float(k/int(input))*100, 2)) + '%')
print('Total unique aptamers: ' + str(len(aptamers)))
print('Aptamers with multi-hits: ' + str(dups))
#print('Total groups: ' + str(len(groups)))		#only with complete hamming check
print('Total positive control reads: ' + str(conts))
print('Percent positive control reads: ' + str(round(float(l/k)*100, 5)) + '%')
print('Total negative control reads: ' + str(m))
print('Skewness of aptamers: ' + str(round(skew, 2)) + '\n')

new = open('groups.txt', 'a+')
outline = sample2 + ' ' + str(groups_perc_sorted) + '\n'
new.write(outline)
new.close

new = open('skews.txt', 'a+')
outline = sample2 + ' ' + str(round(skew, 2)) + '\n'
new.write(outline)
new.close

new = open('unique.txt', 'a+')
outline = sample2 + ' ' + str(round(float(len(aptamers)/k), 5)) + '\n'
new.write(outline)
new.close

new = open('controls.txt', 'a+')
outline = sample2 + ' ' + str(round(float(conts/k)*100, 5)) + '\n'
new.write(outline)
new.close
