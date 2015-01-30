import sys, csv, random

seq_loci = ["CST3", "SERA2", "TEP", "Th2R", "Th3R", "TRAP", "Unnamed", "LD"]
count_loci = ["BEP"]

def collect_matched_sample(clinical_or_cross):
	"""Return a list of samples matching to clinical or cross-sectional"""
	samples = []
	reader = csv.DictReader(open("qdata/clinical/sample_data.tsv"), delimiter='\t')			# list of dicts
	for line in reader:
		sample = line['sample']
		type = line['cross(X)_or_clinical(C)']
		if type == clinical_or_cross:
			samples.append(sample)
	return samples
	
def collect_study_site():
	"""Return a dict that matches samples to study site"""
	study_site = {}
	reader = csv.DictReader(open("adata/RTSSclinicalData.csv"))								# list of dicts
	for line in reader:
		sample = line['id']
		site = line['site']
		study_site[sample] = site
	return study_site	
	
def collect_age_cohort():
	"""Return a dict that matches samples to age cohort"""
	age_cohort = {}
	reader = csv.DictReader(open("adata/RTSSclinicalData.csv"))								# list of dicts
	for line in reader:
		sample = line['id']
		cohort = line['ageCateg'].replace("[", "").replace("]", "")
		age_cohort[sample] = cohort
	return age_cohort		
	
def collect_vaccine_status():
	"""Return a dict that matches samples to vaccine status"""
	vaccine_status = {}
	reader = csv.DictReader(open("adata/RTSSclinicalData.csv"))								# list of dicts
	for line in reader:
		sample = line['id']
		status = line['vaccine']
		vaccine_status[sample] = status
	return vaccine_status	

def collect_data(subjects, seq_data, matched_samples):
	"""Return a dict that matches subject and locus to data"""
	all_subjects = []
	for locus in seq_loci:
		locus_reader = csv.DictReader(open("qdata/sequences/" + locus + ".tsv"), delimiter='\t')		# list of dicts
		for line in locus_reader:
			subject = line['subject'] + "\t" + line['sample']
			if line['sample'] in matched_samples:
				all_subjects.append(subject)
				if (subject, locus) in seq_data:
					seq_data[(subject, locus)].append(line)
				else:
					seq_data[(subject, locus)] = []
					seq_data[(subject, locus)].append(line)
	locus = "BEP"
	locus_reader = csv.DictReader(open("qdata/sequences/" + locus + ".tsv"), delimiter='\t')		# list of dicts
	for line in locus_reader:
		subject = line['subject'] + "\t" + line['sample']
		read_count = int(line['reads'])
		if read_count >= 20:		
			if line['sample'] in matched_samples:		
				all_subjects.append(subject)
				if (subject, locus) in seq_data:
					seq_data[(subject, locus)].append(line)
				else:
					seq_data[(subject, locus)] = []
					seq_data[(subject, locus)].append(line)			
	unique_subjects = set(all_subjects)		
	subjects.extend(list(unique_subjects))
	
def multiplicity(subjects, seq_data):
	"""Mark (subject, locus) with number of infections"""
	marks = {}
	for subject in subjects:
		multiplicities = []
		for locus in seq_loci:
			differences = []
			if (subject, locus) in seq_data:
				for line in seq_data[(subject, locus)]:
					differences.append(int(line['pep_3D7_hamming']))
			if len(differences) > 0:
				mark = len(differences)
				multiplicities.append(mark)
		if len(multiplicities) > 0:
			marks[(subject, "NA")] = max(multiplicities)
		else:
			marks[(subject, "NA")] = -1
	return marks		
	
def subsample_data(subjects, seq_data):
	"""Select a single haplotype if multiple are present in an individual"""
	for subject in subjects:
		for locus in seq_loci:	
			if (subject, locus) in seq_data:
				data = seq_data[(subject, locus)]
				seq_data[(subject, locus)] = [random.choice(data)]

def identical_3D7(subjects, seq_data):
	"""Mark (subject, locus) as 1 if infection matches 3D7 exactly"""
	mark_data = {}
	for subject in subjects:
		for locus in seq_loci:
			marks = []
			if (subject, locus) in seq_data:
				for line in seq_data[(subject, locus)]:
					hamming = int(line['pep_3D7_hamming'])
					mark = 0
					if hamming == 0:
						mark = 1
					marks.append(mark)
			mark_data[(subject, locus)] = marks
	return mark_data
		
def hamming_3D7(subjects, seq_data):
	"""Mark (subject, locus) Hamming distance to 3D7"""
	mark_data = {}
	for subject in subjects:
		for locus in seq_loci:
			marks = []
			if (subject, locus) in seq_data:
				for line in seq_data[(subject, locus)]:
					hamming = int(line['pep_3D7_hamming'])
					marks.append(hamming)
			mark_data[(subject, locus)] = marks
	return mark_data
	
def repeat_category(subjects, seq_data):
	"""Mark (subject, locus) as 0 if repeat count is 39 or less and 1 if repeat count is 40 or greater"""
	locus = "BEP"
	marks = {}	
	for subject in subjects:
		repeats_list = []
		if (subject, locus) in seq_data:
			for line in seq_data[(subject, locus)]:
				repeats_list.append(int(line['repeats']))
		mark = -1
		if len(repeats_list) > 0:
			if repeats_list[0] <= 39:
				mark = 0
			if repeats_list[0] >= 40:
				mark = 1
		marks[(subject, locus)] = mark
	return marks
	
def repeat_count(subjects, seq_data):
	"""Mark (subject, locus) with repeat count"""
	locus = "BEP"
	marks = {}	
	for subject in subjects:
		repeats_list = []	
		if (subject, locus) in seq_data:
			for line in seq_data[(subject, locus)]:
				repeats_list.append(int(line['repeats']))	
		mark = -1		
		if len(repeats_list) > 0:
			mark = repeats_list[0]
		marks[(subject, locus)] = mark
	return marks	
		
def identical_3D7_site(subjects, seq_data, locus, index):
	"""Mark (subject, locus) as 1 if infection matches 3D7 exactly"""
	mark_data = {}
	for subject in subjects:
		marks = []
		if (subject, locus) in seq_data:
			for line in seq_data[(subject, locus)]:
				hamming = int(line["pep_3D7_hamming_" + str(index)])
				mark = 0
				if hamming == 0:
					mark = 1
				marks.append(mark)
		mark_data[(subject, locus)] = marks
	return mark_data
	
def hamming_3D7_sig_sites(subjects, seq_data):
	"""Mark (subject, locus) as with the number of differences to 3D7 at 7 signature sites"""
	# sig sites are 299, 301, 317, 354, 356, 359, 361 in canonical CSP numbering
	# sig sites are 6, 8, 24, 61, 63, 66, 68 in amplicon numbering
	sig_sites = (6, 8, 24, 61, 63, 66, 68)
	locus = "TEP"
	mark_data = {}
	for subject in subjects:
		marks = []
		if (subject, locus) in seq_data:
			for line in seq_data[(subject, locus)]:
				mark = 0
				for site in sig_sites:
					hamming = int(line["pep_3D7_hamming_" + str(site)])
					mark += hamming
				marks.append(mark)
		mark_data[(subject, locus)] = marks
	return mark_data			
		
def print_marks(subjects, seq_data, mark_names, mark_data, study_site, age_cohort, vaccine_status):
	line = ["subject", "sample", "sample_count", "locus", "mark_name", "mark_value", "study_site", "age_cohort", "vaccine_status"]
	print "\t".join(line)
	for subject in subjects:
	
		subject_id = subject.split("\t")[0]
		
		# hamming_3D7 and match_3D7	
		for locus in seq_loci:
			for mark_name in mark_names:
				data = mark_data[mark_name]
				marks = data[(subject, locus)]
				for index, mark in enumerate(marks):
					line = [subject, str(index+1), locus, mark_name, str(mark), str(study_site[subject_id]), str(age_cohort[subject_id]), str(vaccine_status[subject_id])]
					print "\t".join(line)
					
		# site-specific match
		locus = "TEP"
		for index in range(294,388):
			mark_name = "match_3D7_" + str(index)
			data = mark_data[mark_name]
			marks = data[(subject, locus)]
			for index, mark in enumerate(marks):
				line = [subject, str(index+1), locus, mark_name, str(mark), str(study_site[subject_id]), str(age_cohort[subject_id]), str(vaccine_status[subject_id])]
				print "\t".join(line)	
				
		# SERA site-specific match
		locus = "SERA2"
		for index in range(36,120):
			mark_name = "match_3D7_" + str(index)
			data = mark_data[mark_name]
			marks = data[(subject, locus)]
			for index, mark in enumerate(marks):
				line = [subject, str(index+1), locus, mark_name, str(mark), str(study_site[subject_id]), str(age_cohort[subject_id]), str(vaccine_status[subject_id])]
				print "\t".join(line)						

		# hamming_3D7_sig_sites
		locus = "TEP"
		mark_name = "hamming_3D7_sig_sites"
		data = mark_data[mark_name]
		marks = data[(subject, locus)]
		for index, mark in enumerate(marks):
			line = [subject, str(index+1), locus, mark_name, str(mark), str(study_site[subject_id]), str(age_cohort[subject_id]), str(vaccine_status[subject_id])]
			print "\t".join(line)	
								
		# multiplicity
		mark_name = 'multiplicity'
		data = mark_data[mark_name]
		mark = str(data[(subject, "NA")])			
		line = [subject, "1", "NA", mark_name, mark, str(study_site[subject_id]), str(age_cohort[subject_id]), str(vaccine_status[subject_id])]
		print "\t".join(line)		
				
		# repeat_category		
		mark_name = 'repeat_category'
		data = mark_data[mark_name]
		mark = str(data[(subject, "BEP")])
		line = [subject, "1", "BEP", mark_name, mark, str(study_site[subject_id]), str(age_cohort[subject_id]), str(vaccine_status[subject_id])]
		print "\t".join(line)				
		
		# repeat_count
		mark_name = 'repeat_count'
		data = mark_data[mark_name]
		mark = str(data[(subject, "BEP")])
		line = [subject, "1", "BEP", mark_name, mark, str(study_site[subject_id]), str(age_cohort[subject_id]), str(vaccine_status[subject_id])]
		print "\t".join(line)						
			
def main(argv):
	"""Construct mark data"""
	clinical_or_cross = "C"
	if len(argv) > 0:
		clinical_or_cross = argv[0]
	matched_samples = collect_matched_sample(clinical_or_cross)
	study_site = collect_study_site()	
	age_cohort = collect_age_cohort()	
	vaccine_status = collect_vaccine_status()	
	subjects = []
	seq_data = {}	
	collect_data(subjects, seq_data, matched_samples)									# sets subjects and data
	mark_names = ['match_3D7', 'hamming_3D7']
	mark_data = {}
	mark_data['multiplicity'] = multiplicity(subjects, seq_data)
	mark_data['match_3D7'] = identical_3D7(subjects, seq_data)
	mark_data['hamming_3D7'] = hamming_3D7(subjects, seq_data)
	mark_data['repeat_category'] = repeat_category(subjects, seq_data)
	mark_data['repeat_count'] = repeat_count(subjects, seq_data)
	for index in range(1,95):
		mark_data["match_3D7_" + str(index+293)] = identical_3D7_site(subjects, seq_data, "TEP", index)
	for index in range(1,85):
		mark_data["match_3D7_" + str(index+35)] = identical_3D7_site(subjects, seq_data, "SERA2", index)
	mark_data['hamming_3D7_sig_sites'] = hamming_3D7_sig_sites(subjects, seq_data)				
	print_marks(subjects, seq_data, mark_names, mark_data, study_site, age_cohort, vaccine_status)

if __name__ == "__main__":
    main(sys.argv[1:])