
from somacaller import model
import io 
import tempfile
import os
import shutil

BED_FILE = "examples/hotspot.bed" 
BAM_FILE = "examples/tpos.bam"


def test_read_bedfile():

	assert os.path.exists(BED_FILE)

	df = model.read_bedfile(BED_FILE)
	assert list(df.columns) == ["chrom","start", "end", "name"]
	with open(BED_FILE) as file:
		line_count = len([i for i in file])
		assert len(df) == line_count


def test_slop_bedfile():

	assert os.path.exists(BED_FILE)

	for slop in range(1, 10):
		newfile = model.slop_bedfile(BED_FILE, slop)

		with open(BED_FILE) as file:
			line_count = len([i for i in file])

		with open(newfile) as file:
			new_line_count = len([i for i in file])
		assert line_count * ( slop * 2 + 1) == new_line_count

def test_binary_exist():
	assert shutil.which("samtools")
	assert shutil.which("bam-readcount")
	

def test_bam_read_count():
	assert "HG19" in os.environ, "please set $HG19 global environement variable with hg19.fa path"
	assert os.path.exists(BAM_FILE)
	HG19_FILE = os.environ["HG19"]

	bed_df = model.read_bedfile(BED_FILE)
	# Testing with hotspot 
	df = model.bam_read_count(BAM_FILE, HG19_FILE, BED_FILE)

	assert len(df) == len(bed_df)
	assert int(df.query("chrom == 'chr7' & pos == '55241707'")["depth"]) == 1776

	assert list(df.columns) == ["chrom","pos","name","depth","A","C","G","T","file","sample","first_vaf","second_vaf", "first", "second"]

	print(df)

	# chr7   55241707    EGFR  1776.0   328.0     0.0  1448.0     0.0  tpos.bam  SM:T POS ATC           1448.0             328.0


def test_async_bam_read_count():
	assert "HG19" in os.environ, "please set $HG19 global environement variable with hg19.fa path"
	assert os.path.exists(BAM_FILE)
	HG19_FILE = os.environ["HG19"]	
	bed_df = model.read_bedfile(BED_FILE)

	first_bam = BAM_FILE
	second_bam = tempfile.mkstemp(suffix = ".bam")[1]
	files = [first_bam, second_bam]
	shutil.copyfile(BAM_FILE, second_bam)
	shutil.copyfile(BAM_FILE+".bai", second_bam + ".bai")
	df = model.async_bam_read_counts(files, HG19_FILE, BED_FILE, threads=4)

	assert len(df) == len(files) * len(bed_df)




