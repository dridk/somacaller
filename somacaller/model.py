import pandas as pd 
import os
import subprocess 
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor
from io import StringIO
import functools
import time 
from glob import glob
import csv
import tempfile
import os 
import seaborn as sns
import logging
import numpy as np

from scipy import stats

import altair as alt
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.linear_model import BayesianRidge
from sklearn.neighbors import LocalOutlierFactor, NearestNeighbors, KNeighborsRegressor
from sklearn.pipeline import make_pipeline

import matplotlib.pyplot as plt

HG19 =  "/mnt/R60/bioinformatique/hg19/hg19.fasta"
HOTSPOT = "../../hotspot.bed"
BAMFILE = "../../IonXpress_032_R_2019_05_16_08_36_56_user_OUE-855-2019-05-15_P2_Cancero_ampliseq_318_mnk_Auto_user_OUE-855-2019-05-15_P2_Cancero_ampliseq_318_mnk_254.bam"
BAMFILE2 = "../../IonXpress_096_R_2019_12_12_14_10_53_user_OUE-968-2019-12-11_P2-S50_Cancero_ampliseq_318_Auto_user_OUE-968-2019-12-11_P2-S50_Cancero_ampliseq_318_377.bam"

def read_bedfile(bed_file: str) -> pd.DataFrame:
    """Read bed file and the corresponding dataframe 
    
    Bedfile must be a tab separated with 4 columns:
        chr, start, end, gene

    Note: 
        Start position will be converted to 1-based 
    
    Args:
        bed_file (str): Path to bedfile
    
    Returns:
        pd.DataFrame: A dataframe with 4 columns (chrom,start,end,name)
    """
    df = pd.read_csv(bed_file, sep="\t", header=None)

    df.columns = ["chrom","start", "end", "name"]
    df["start"] = df["start"] + 1
    df["id"] = df["chrom"] + ":" + df["start"].astype(str) 
    df = df.set_index("id").drop_duplicates()
    return df

def slop_bedfile(bed_file:str, slop = 5) -> str:
    """Slop a bedfile from right and left 
    
    Args:
        bed_file (str): bedfile 
        slop (int, optional): slop , by default 30
    
    Returns:
        str: return a new bedfile path
    """
    large_bedfile = tempfile.mkstemp(suffix = ".bed")[1]
    with open(bed_file) as file_in, open(large_bedfile,"w") as file_out :
        reader =  csv.reader(file_in, delimiter="\t")
        writer = csv.writer(file_out, delimiter="\t")

        for line in reader:
            start = int(line[1])
            end = int(line[2])

            for i in range(-slop, slop + 1):
                new_line = [line[0], start + i , end + i, line[3]]
                writer.writerow(new_line)

    return large_bedfile


def bam_read_count(bamfile:str, reference_file: str, bed_file:str) -> pd.DataFrame:
    """Call bam-readcount and return a data frame 
    
    Args:
        bamfile (str): an indexes bam file  *.bam ( with *.bai)
        reference_file (str): a reference genom *.fa
        bed_file (str): a bed file *.bed
    
    Returns:
        pd.DataFrame : Return a table with chr, pos, ref, depth, A,C,G,T
    """

    # Run bam-readcount
 
    logging.debug("process bamfile: " + bamfile)
    cmd = f"samtools view -H {bamfile}|grep -oE \"SM:.+\""
   
    sample = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
    logging.debug("process sample: " + sample)
    cmd = f"bam-readcount -f  {reference_file} {bamfile} -l {bed_file}|tr \":\" \"\t\"|cut -f1,2,3,4,20,34,48,62 "
    data = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
 
    # Capture stdout 
    try:
        data = StringIO(data)
        df = pd.read_csv(data,sep="\t", header=None)
        df.columns = ["chrom","pos","ref","depth","A","C","G","T"]

    except:
        print("something is wrong with", bamfile)
        df = pd.DataFrame({"chrom":[], "pos":[],"ref":[], "depth":[], "A":[],"C":[],"G":[],"T":[]}) 
    
    # Set ids 
    df["id"] = df["chrom"] +":"+ df["pos"].astype(str)
    df = df.set_index("id")

    #  Fill uncalled variant from bedfile 
    bed_df = read_bedfile(bed_file).copy()
    final_df = bed_df.join(df[["depth","A","C","G","T"]]).fillna(0)
    final_df = final_df.drop("end", axis=1)
    
    # Compute extra columns 
    final_df["file"] = os.path.basename(bamfile)
    final_df["sample"] = os.path.basename(sample)

    final_df["first_vaf"] =  final_df.filter(regex="(A|C|G|T)").apply(lambda x: x.drop_duplicates().max(), axis=1)
    final_df["second_vaf"] = final_df.filter(regex="(A|C|G|T)").apply(lambda x: x.nlargest(2)[1], axis=1)

    final_df["first"]  = final_df.apply(lambda s: s[["A","C","G","T"]].astype(int).nlargest(2).index[0], axis=1)
    final_df["second"] = final_df.apply(lambda s: s[["A","C","G","T"]].astype(int).nlargest(2).index[1], axis=1)

    final_df = final_df.rename({"start":"pos"}, axis=1)

    return final_df.reset_index(drop=True).drop_duplicates()



def async_bam_read_counts(bamfiles: list, reference_file: str, bed_file, threads = 10):
    """Execute bam_read_count asynchronly on multiple bam files 
    
    Args:
        bamfiles (list): Description
        reference_file (str): Description
        hotspot_file (TYPE): Description
        threads (int, optional): Description
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        fct = functools.partial(bam_read_count, reference_file=reference_file, bed_file=bed_file)
        # bamfiles = glob("../../*.bam")[:10]

        results = list(executor.map(fct, bamfiles))

    all_df = pd.concat(results,ignore_index = True).reset_index() 
    all_df["pos"]  = all_df["pos"].astype(int)
    all_df["depth"]  = all_df["depth"].astype(int)
    all_df["A"]  = all_df["A"].astype(int)
    all_df["C"]  = all_df["C"].astype(int)
    all_df["G"]  = all_df["G"].astype(int)
    all_df["T"]  = all_df["T"].astype(int)

    return all_df







class SomaModel(object):

    def __init__(self, hotspot_file:str, reference_file:str):
        self.hotspot_file = os.path.abspath(hotspot_file)
        self.reference_file = os.path.abspath(reference_file)

        #  Compute hotspot df 
        self.hotspot = pd.read_csv( self.hotspot_file, sep="\t", header=None)
        self.hotspot.columns = ["chr","start","end", "gene"]
        self.hotspot["start"] = self.hotspot["start"] + 1 
        self.hotspot["id"] = self.hotspot["chr"] + ":" + self.hotspot["start"].astype(str)
        self.hotspot= self.hotspot.set_index("id")

        self.linear_models = {}
        self.outlier_models = {}
        self.scalers = {}

    def fit(self, bamlist, threads = 15):


        # Compute new slop bed 
        slop_bedfile = slop_bedfile(self.hotspot_file, 5)

        # Execute thread map/reduce
        self.raw_data = async_bam_read_counts(bamlist, reference_file = self.reference_file, bed_file = slop_bedfile, threads = threads)


    def _create_models(self):

        target_data = self.target_data()

        for index in self.hotspot.index:
            #chrom, pos = index.split(":")
            logging.info(f"create linear model: {index}")
            # subset df selection 
            df = target_data[target_data["id"] == index][["depth", "second_vaf"]]
            df.columns = ["x","y"]
            
            # Linear regression
            X = df["x"].values.reshape(-1,1)
            Y = df["y"].values
            lr =  BayesianRidge()
            lr.fit(X,Y)
            self.linear_models[index] = lr

            #Outlier model 
            clf = make_pipeline(MinMaxScaler(), NearestNeighbors(10))
            clf.fit(df)
            self.outlier_models[index] = clf

    def plot(self, filename,position):
        
        data = self.target_data().query(f"id == '{position}'")

        x = np.arange(data["depth"].min(), data["depth"].max(), 1000).reshape(-1,1)
        y = self.linear_models[position].predict(x)

        line = pd.DataFrame({"x":x.reshape(-1), "y":y.reshape(-1)})

        c1 = alt.Chart(data).mark_point().encode(x="depth", y="second_vaf", tooltip = ["file","sample"])
        c2 = alt.Chart(line).mark_line().encode(x="x",y="y")

        return c1 + c2

        # sns.scatterplot(x="depth", y="second_vaf", data=data)

        # plt.savefig(filename)



            # data = hotspot.query("id == 'chr12:25398280'")

            # X = data["depth"].values.reshape(-1,1)
            # Y = data["second_vaf"].values




    def to_hdf(self, filename):
        self.raw_data.to_hdf(filename, key='raw_data', mode='w')
        #self.hotspot_data.to_hdf(filename, key="hotspot_data", mode="a")

        parameters = pd.Series({
            "reference_file": self.reference_file,
            "hotspot_file": self.hotspot_file
            })
        
        parameters.to_hdf(filename, key='parameters', mode='a')


    @classmethod
    def from_hdf(cls,filename):

        raw_data = pd.read_hdf(filename, 'raw_data')
        #hotspot_data = pd.read_hdf(filename, 'hotspot_data')
        parameters = pd.read_hdf(filename, "parameters")

        reference_file = parameters["reference_file"]
        hotspot_file = parameters["hotspot_file"]

        model =  SomaModel(hotspot_file, reference_file)   
        model.raw_data = raw_data
        return model



    def test(self, bamfile):
        df = bam_read_count(bamfile, self.reference_file, self.hotspot_file)
        target_df = self._create_target(df)
        target_df = target_df[target_df["depth"] > 10].copy()

        target_df["af"] = target_df["second_vaf"] / target_df["depth"] * 100

        # Compute regression score 
        target_df["z_linear"] = target_df.apply(lambda s : self._predict_regression_score(s["depth"], s["second_vaf"], s["id"]), axis=1)
        target_df["outlier"] = target_df.apply(lambda s : self._predict_outlier_score(s["depth"], s["second_vaf"], s["id"]), axis=1)


        # Compute outliers 


        return target_df


    def _predict_regression_score(self, x, y, position):
        #df = sample_df.query(f"id == '{position}'")
        all_df = self.target_data().query(f"id == '{position}'")
        
        std = all_df["second_vaf"].std()

        x = np.array([x]).reshape(-1,1)
        y_pred = self.linear_models[position].predict(x).reshape(-1)[0]
        
        return abs(y_pred - y) / std 


    def _predict_outlier_score(self, x, y, position):
        #df = sample_df.query(f"id == '{position}'")
        
        df = pd.DataFrame({"x":[x], "y":[y]})
        df = self.outlier_models[position]["minmaxscaler"].transform(df)
    
        distance, indices = self.outlier_models[position]["nearestneighbors"].kneighbors(df)

        mean_distance = np.mean(distance[:,1:], axis=1)[0]

        return mean_distance
        #return self.outlier_models[position].fit(df, novelty=True)[0]

        
        


        #return (1 - stats.norm(0,1).cdf(abs(y - y_pred)/std)) * 100

  
        # index = f"{chrom}:{pos}"


        # data = self.hotspot_data.loc[index][["depth", "second_vaf"]]
        # data = pd.DataFrame(MinMaxScaler().fit_transform(data))
        # data.columns = ["x","y"]

        # sns.scatterplot("x","y", data=data)
        # plt.savefig("test.png")

        # df.loc[f"{chrom}:{pos}"]





    def data(self):
        return self.raw_data
      
    def target_data(self):
     # Create index 
        return self._create_target(self.raw_data)

    def _create_target(self, dataframe):
        dataframe["id"] = dataframe["chr"] + ":" + dataframe["pos"].astype(str)
        target_data = dataframe[dataframe["id"].isin(self.hotspot.index)].copy()
        return target_data.reset_index(drop=True)
       




#bam_read_count(BAMFILE, HG19, HOTSPOT)

if __name__ == '__main__':
    



    from glob import glob 

    model = SomaModel(HOTSPOT, HG19)

    model.fit(glob("../../*.bam"))

    model.to_hdf("data.h5")


    # model = SomaModel.from_hdf("data.h5")

    # testfile = "../../IonXpress_082_R_2019_07_11_09_01_37_user_OUE-885-2019-07-10_S28_P2_Cancero_ampliseq_318_Auto_user_OUE-885-2019-07-10_S28_P2_Cancero_ampliseq_318_285.bam"
    # model.test_position(testfile,"chr7", 55241707)

#  

#  somacaller fit -i *.bam -target hotspot.bed > model.h5 
#  somacaller test -i file.bam -m model.h5 > output.csv 
#  somacaller plot -r output.csv -m model.h5 -r "chr3" > output.png 
