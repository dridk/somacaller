
from somacaller.model import *
import io 
import tempfile
import os
import shutil
import pytest

BED_FILE = "examples/hotspot.bed" 
BAM_FILE = "examples/tpos.bam"

HG19_FILE = os.environ["HG19"]  




def create_bam_files():
    first_bam = BAM_FILE
    second_bam = tempfile.mkstemp(suffix = ".bam")[1]
    shutil.copyfile(first_bam, second_bam)
    shutil.copyfile(first_bam+".bai", second_bam + ".bai")
    return [first_bam]    


@pytest.fixture(scope="session")
def model():
    print("test")
    HG19_FILE = os.environ["HG19"]  
    model = SomaModel(BED_FILE, HG19_FILE)
    model.fit(create_bam_files())
    return model
    


def test_read_bedfile():

    assert os.path.exists(BED_FILE)

    df = read_bedfile(BED_FILE)
    assert list(df.columns) == ["chrom","start", "end", "name"]
    with open(BED_FILE) as file:
        line_count = len([i for i in file])
        assert len(df) == line_count


def test_slop_bedfile():

    assert os.path.exists(BED_FILE)

    for slop in range(1, 10):
        newfile = slop_bedfile(BED_FILE, slop)

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

    bed_df = read_bedfile(BED_FILE)
    # Testing with hotspot 
    df = bam_read_count(BAM_FILE, HG19_FILE, BED_FILE)

    assert len(df) == len(bed_df)
    assert int(df.query("chrom == 'chr7' & pos == '55241707'")["depth"]) == 342

    assert list(df.columns) == ["chrom","pos","name","depth","A","C","G","T","file","sample","first_vaf","second_vaf", "first", "second"]

    print(df)


def test_entropy():

    assert entropy(pd.Series([0.25, 0.25, 0.25, 0.25])) == 2
    assert entropy(pd.Series([1, 0, 0, 0])) == 0




def test_serialization(model):
    _, save_file = tempfile.mkstemp(suffix=".h5")
    model.to_hdf(save_file)

    # create new model 




def test_info(model):

    items = model.info()

    assert items["hotspot_file"] == os.path.abspath(BED_FILE)
    assert items["reference_file"] == os.path.abspath(HG19_FILE)
    assert items["bam_files"] == create_bam_files()



# def test_model(model):

#     HG19_FILE = os.environ["HG19"]  



#     my_model = SomaModel(BED_FILE, HG19_FILE)
#     my_model.fit([BAM_FILE])

#     _, save_file = tempfile.mkstemp(suffix=".h5")
#     my_model.to_hdf(save_file)

#     new_model = SomaModel.from_hdf(save_file)

#     assert new_model.raw_data == my_model.raw_data





