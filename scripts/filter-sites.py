# Requires files:
#	adata/marks_data_c_sites.tsv
#	adata/marks_data_x_sites.tsv

# TEP sites are numbered 294 to 387
# SERA2 sites are numbered 36 to 119

import sys, csv

def get_data(endpoint):
	"""Return list of dicts for entire dataset"""
	data = []
	reader = csv.DictReader(open("adata/marks_data_" + endpoint.lower() + "_sites.tsv"), delimiter='\t')
	for line in reader:
		data.append(line)
	return data

def subset_data(data, locus, sites, cohort):
	"""Get match data subsetted by locus, study sites, age cohort and endpoint"""
	subset = []
	for line in data:
		if line["locus"].lower() == locus.lower() and line["age_cohort"].lower() == cohort.lower() and line["study_site"] in sites:
			subset.append(line)
	return subset

def count_match_mismatch(data, index):
	"""Count subjects with match in at least 1 parasite"""
	subjects = set()
	new_data = []
	for line in data:
		if line["mark_name"] == "match_3D7_" + str(index):
			subjects.add(line["subject"])
			new_data.append(line)
	match_map = {}
	mismatch_map = {}	
	for subject in subjects:
		match_map[subject] = 0
		mismatch_map[subject] = 0		
		for line in new_data:
			if line["subject"] == subject and line["mark_value"] == "1":
				match_map[subject] = 1
			if line["subject"] == subject and line["mark_value"] == "0":
				mismatch_map[subject] = 1				
	match_count = 0
	mismatch_count = 0	
	for subject, value in match_map.items():
		match_count += value
	for subject, value in mismatch_map.items():
		mismatch_count += value		
	return (match_count, mismatch_count)	
		
def get_filtered_sites(data, start, end, threshold):
	"""Filter sites"""
	sites = set()
	for index in range(start, end+1):
		(match_count, mismatch_count) = count_match_mismatch(data, index)	
		if match_count < mismatch_count:
			if match_count >= threshold:		
				print str(index) + "\t" + str(match_count) + "\t" + str(mismatch_count) + "\t"					
				sites.add(index)
		if mismatch_count < match_count:
			if mismatch_count >= threshold:		
				print str(index) + "\t" + str(match_count) + "\t" + str(mismatch_count) + "\t"					
				sites.add(index)
	return sites
	
def print_filtered_sites(locus, site_group, cohort, aa_sites):
	"""Print filtered sites"""
	for aa in aa_sites:
		print "\t".join([locus, site_group, cohort, str(aa)]) 

def main(argv):
	"""Filter amino acid sites"""
	endpoint = "C"
	if len(argv) > 0:
		endpoint = argv[0]		
	
	site_map = {}
	site_map["5 Major Pooled"] = ["Agogo", "Kintampo", "Kombewa", "Nanoro", "Siaya"]
	site_map["11 Pooled"] = ["Agogo", "Kintampo", "Kombewa", "Nanoro", "Siaya", "Bagamoyo", "Kilifi", "Korogwe", "Lambarene", "Lilongwe", "Manhica"]
	start_map = {}
	start_map["TEP"] = 294
	start_map["SERA2"] = 36
	end_map = {}
	end_map["TEP"] = 387
	end_map["SERA2"] = 119
	
	threshold_map = {}
	threshold_map["TEP_5 Major Pooled_6-12 Weeks_C"] = 4
	threshold_map["TEP_11 Pooled_6-12 Weeks_C"] = 4
	threshold_map["TEP_5 Major Pooled_5-17 Months_C"] = 4
	threshold_map["TEP_11 Pooled_5-17 Months_C"] = 4	
	threshold_map["TEP_5 Major Pooled_6-12 Weeks_X"] = 3
	threshold_map["TEP_11 Pooled_6-12 Weeks_X"] = 3
	threshold_map["TEP_5 Major Pooled_5-17 Months_X"] = 4
	threshold_map["TEP_11 Pooled_5-17 Months_X"] = 4	
	
	threshold_map["SERA2_5 Major Pooled_6-12 Weeks_C"] = 4
	threshold_map["SERA2_11 Pooled_6-12 Weeks_C"] = 4
	threshold_map["SERA2_5 Major Pooled_5-17 Months_C"] = 4
	threshold_map["SERA2_11 Pooled_5-17 Months_C"] = 4	
	threshold_map["SERA2_5 Major Pooled_6-12 Weeks_X"] = 4
	threshold_map["SERA2_11 Pooled_6-12 Weeks_X"] = 3
	threshold_map["SERA2_5 Major Pooled_5-17 Months_X"] = 4
	threshold_map["SERA2_11 Pooled_5-17 Months_X"] = 4			
		
	data = get_data(endpoint)
	
	print "locus\tsite\tcohort\taa"

	for locus in ["SERA2"]:
		for site_group in ["5 Major Pooled", "11 Pooled"]:
			for cohort in ["6-12 Weeks", "5-17 Months"]:		
				subset = subset_data(data, locus, site_map[site_group], cohort)
				threshold = threshold_map["_".join([locus, site_group, cohort, endpoint])]
				aa_sites = get_filtered_sites(subset, start_map[locus], end_map[locus], threshold)
				print_filtered_sites(locus, site_group, cohort, aa_sites)

if __name__ == "__main__":
    main(sys.argv[1:])
