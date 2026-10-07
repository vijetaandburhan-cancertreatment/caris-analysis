"""C-backed samtools regional filter for an unsorted BAM stream; no indexing required."""
import sys,pysam
if __name__=='__main__':
 bed,source,destination=sys.argv[1:]
 pysam.samtools.view('--no-PG','-b','-L',bed,'-o',destination,source,catch_stdout=False)
