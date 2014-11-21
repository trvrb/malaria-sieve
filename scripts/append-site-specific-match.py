import sys, csv, random

seq_loci = ["TEP", "SERA2", "TRAP"]

def get_site_specific_match(locus, comparison_sequence):
	"""Get data with site specific match appended"""
	reader = csv.DictReader(open("qdata/sequences/" + locus + ".tsv"), delimiter='\t')		# list of dicts
	data = []
	for line in reader:
		data.append(line)
	for line in data:
		seq = line["pep_sequence"]
		for index, (aa1, aa2) in enumerate(zip(comparison_sequence, seq)):
			name = "pep_3D7_hamming_" + str(index+1)
			hamming = 0 if aa1 == aa2 else 1
	#		print name + "\t" + aa1 + aa2 + "\t" + str(hamming)
			line[name] = hamming		
	return data

def print_match(data, site_count):
	"""Print data with site specific match appended"""
	header = ["subject", "sample", "sample_count", "site", "nuc_haplo", "nuc_sequence", "nuc_3D7_hamming", "pep_haplo", "pep_sequence", "pep_3D7_hamming", "reads"]
	extra_header = ["pep_3D7_hamming_" + str(index) for index in range(1, site_count)]
	print "\t".join(header + extra_header)
	for line in data:
		out1 = [line[h] for h in header]
		out2 = [str(line["pep_3D7_hamming_" + str(index)]) for index in range(1, site_count)]
		print "\t".join(out1 + out2)
			
def main(argv):
	"""Construct amended sequence data"""
	
	locus = "TEP"
	if len(argv) > 0:
		locus = argv[0]		
	
	match_3d7_locus = {}
	match_3d7_locus["TEP"] = "DENANANSAVKNNNNEEPSDKHIKEYLNKIQNSLSTEWSPCSVTCGNGIQVRIKPGSANKPKDELDYANDIEKKICKMEKCSSVFNVVNSSIGLI"
	match_3d7_locus["SERA2"] = "LSSDGSRVTTQARIEKPKQQPTLPTLAQETQPQQQQQQKEVGSGIGAEQKVESARPGAEVSQSDVERAGRSSGTGGSVGTKISPG"	
	match_3d7_locus["TRAP"] = "QNNLPNDKSDRYIPYSPLSPKVLDNERKQSDPQSQDNNGNRHVPNSEDRETRPHGRNNENRSYNRKHNNTPKHPEREEHEKPDNNKKKAGS-DNKYKIAGGIAGGL"
	seq_3d7 = match_3d7_locus[locus]
		
	data = get_site_specific_match(locus, seq_3d7)
	print_match(data, len(seq_3d7))
	
if __name__ == "__main__":
    main(sys.argv[1:])