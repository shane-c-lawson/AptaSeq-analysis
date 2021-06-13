from scipy.spatial.distance import hamming

dicts = {}
with open('groups.txt', 'r') as out :
	for line in out :
		dicts[line.split(' ')[0]] = eval(' '.join(line.strip().split(' ')[1:]))

groups = {}; hit_ever = False
for dict in dicts.values() :
	for aptamer in dict.keys() :
		hit = False
		for group in groups.keys() :
			dist = hamming(list(group),list(aptamer))
			if dist <= 0.1 :
				hit = True; hit_ever = True
				groups[group] += 1
				break
		if hit == False :
			groups[aptamer] = 1

best_skew = 0; worst_skew = 9999999999; bs = False; ws = False
with open('skews.txt', 'r') as skews :
	for line in skews :
		if float(line.strip().split(' ')[1]) > best_skew :
			best_skew = float(line.strip().split(' ')[1])
			best_skew_sample = line.split(' ')[0]
			bs = True
		if float(line.strip().split(' ')[1]) < worst_skew :
			worst_skew = float(line.strip().split(' ')[1])
			worst_skew_sample = line.split(' ')[0]
			ws = True

if bs == True :
	print('Best distribution: ' + best_skew_sample + ' (' + str(best_skew) + ' skewness)')
if ws == True :
	print('Worst distribution: ' + worst_skew_sample + ' (' + str(worst_skew) + ' skewness)' + '\n')

best_cont = 0; worst_cont = 9999999999; bc = False; wc = False
with open('controls.txt', 'r') as conts :
	for line in conts :
		if float(line.strip().split(' ')[1]) > best_cont :
			best_cont = float(line.strip().split(' ')[1])
			best_cont_sample = line.split(' ')[0]
			bc = True
		if float(line.strip().split(' ')[1]) < worst_cont :
			worst_cont = float(line.strip().split(' ')[1])
			worst_cont_sample = line.split(' ')[0]
			wc = True
	
if bc == True :
	print('Best control performance: ' + best_cont_sample + ' (' + str(best_cont) + ' controls per read)')
if wc == True :
	print('Worst control performance: ' + worst_cont_sample + ' (' + str(worst_cont) + ' controls per read)' + '\n')

best_div = 9999999999 ; worst_div = 0; bd = False; wd = False
with open('unique.txt', 'r') as conts :
	for line in conts :
		if float(line.strip().split(' ')[1]) < best_div :
			best_div = float(line.strip().split(' ')[1])
			best_div_sample = line.split(' ')[0]
			bd = True
		if float(line.strip().split(' ')[1]) > worst_div :
			worst_div = float(line.strip().split(' ')[1])
			worst_div_sample = line.split(' ')[0]
			wd = True
	
if bd == True :
	print('Best (lowest) diversity: ' + best_div_sample + ' (' + str(best_div) + ' unique aptamers per read)')
if wd == True :
	print('Worst (highest) diversity: ' + worst_div_sample + ' (' + str(worst_div) + ' unique aptamers per read)' + '\n')

if hit_ever == True :
	print('Aptamers with duplicate hits in multiple samples')
	for group in groups :
		if groups[group] > 1 :
			print(group + ': ' + str(groups[group]) + ' samples')
