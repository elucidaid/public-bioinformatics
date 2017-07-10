with open("/Users/ebrahimih/Desktop/Schizosaccharomyces_pombe.ASM294v2.30.gff3", 'r') as file, open('/Users/ebrahimih/Desktop/new.gff', 'w') as output_file:
	for line in file:
		columns = line.split("\t")
		if columns[0]=="II":
			if int(columns[3]) >= 2137121:
				columns[3] = str(int(columns[3]) + diff)
				columns[4] = str(int(columns[4]) + diff)
		output_file.write('\t'.join(columns))