import sys

nuc_prio = {'C' : 0, 'A' : 1, 'T' : 2, 'G' : 3}
nuc_comp = {'C' : 'G', 'A' : 'T', 'T' : 'A', 'G' : 'C'}


def orient_file(old_path, new_path):
	#l = 0
	with open(old_path, 'r') as rfile, open(new_path, 'w') as wfile:
		for line in rfile:
			#l+= 1
			#if l%10000==0:
			#	print(l)
			oriented_line = orient(clean_line(line))
			wfile.write(oriented_line+"\n")
			
	
def clean_line(line):
	if line:
		return line.split(" ")[0].strip()
	else:
		return ""


def orient(kmer):
	if len(kmer) == 0:
		return ""
	complement = ""
	for nuc in kmer:
		if nuc not in nuc_prio:
			return ""
		complement = nuc_comp[nuc] + complement
	for i in range(len(kmer)):
		if nuc_prio[kmer[i]] == nuc_prio[complement[i]]:
			continue
		if nuc_prio[kmer[i]] > nuc_prio[complement[i]]:
			return kmer
		if nuc_prio[kmer[i]] < nuc_prio[complement[i]]:
			return complement
	return kmer


def main(old_path, new_path):
	print("Orienting k-mer file\n")
	orient_file(old_path, new_path)


if __name__ == "__main__":
   main(sys.argv[1], sys.argv[2])


