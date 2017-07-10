import pysam, numpy, argparse, os, math, sys
import pybedtools
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--bam", help="bam input file to normalize", required=True)
parser.add_argument("--bw", help="bw output file name", required=True)
args = parser.parse_args()
# generate BAM index if absent
if not os.path.isfile(args.bam+'.bai'):
	cmd = "/anaconda/bin/samtools index %s"%args.bam
	sys.stderr.write("[%s]  Indexing BAM file: %s\n"%(datetime.ctime(datetime.now()), cmd))
	os.system(cmd)
sam_file = pysam.AlignmentFile(args.bam)
totalNumReads = sam_file.mapped
print totalNumReads
control_array=list()
bkgnd_array=list()
binsize = 1
stepsize = 1
bin_index = 0
binned_control_reads = numpy.zeros(1744)
sizes = zip(sam_file.references, sam_file.lengths)
regions = [(name, 0, length) for name, length in sizes]
for col in sam_file.pileup("III", 1304000,1380000):
	n = float(col.n)
	bkgnd_array.append(n)
for col in sam_file.pileup("III", 1316291, 1318035):
	n = float(col.n)
	control_array.append(n)
#startpos=1315499
#for endpos in xrange( 1316291+stepsize, 1318035, stepsize ):
#	hits = 0
#	for col in sam_file.pileup("III", startpos,endpos):
#		hits += col.n
#	#print ("%s-%s: %s" % (str(startpos),str(endpos),str(hits)))
#	#binned_control_reads[bin_index] = (hits*1e6/total)
#	binned_control_reads[bin_index] = hits
#	startpos += stepsize
#	bin_index += 1
#total = 0
#for col in sam_file.pileup("III", 1259763,1299991):
#	total += col.n
#print total
bkgnd_median = numpy.median(numpy.array(bkgnd_array).astype(numpy.float))
ade6_mapped = numpy.array(control_array).astype(numpy.float)
ade6_median = numpy.median(ade6_mapped)
ade6_normalized =  ade6_median / 1744
bckgrnd_normalized = bkgnd_median / 76000
#_scale = 1 / ade6_normalized
_scale =  10 / (ade6_median - bkgnd_total)
print _scale
a = pybedtools.BedTool(args.bam).genome_coverage (d=True, scale=(_scale), output= args.bw)
#.saveas(args.bw)
#	bw.write("track %s\n" % " ".join(["type=wiggle_0",
#	        "name=%s" % os.path.splitext(os.path.split(args.bam)[-1])[0],
#	        "visibility=full",
#	        ]))
#	for chrom, start, end in regions:
#		if end is None and chrom in sam_file.references:
#			 end = sam_file.lengths[sam_file.references.index(chrom)]
#		assert end is not None, "Could not find %s in header" % chrom
#		wig.write("variableStep chrom=%s\n" % chrom)
#		for col in sam_file.pileup(chrom, start, end):
		#	n = ((float(col.n) / total * 1e6)-bkgnd) * factor
#			n = math.log (float(col.n) / control_median, 2)
		#	n = (float(col.n) * factor)
		#	if n >= 20:
		#		n = 20
#			wig.write("%s %.1f\n" % (col.pos+1, n))
#			is_valid = True
