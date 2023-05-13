import sys,os
import pandas as pd
import json
import re
#import requests
import json

sys.path.append("..")
from utils.basic import Basic

"""
this script is created on 20211002
"""

#data collection
def super_enhancer_data_collection(meta_data_dir = "/home/mwshi/project/CEdb/cell_line/metajson/", result_dir = "/home/zhluo/Project/TF_enrichment/test/work_dir_all_data", super_enhancer_dir=None):
    sample_list = []
    for one_json in os.listdir(meta_data_dir):
        sample_list.append(one_json.split(".")[0])
        json_file = os.path.join(meta_data_dir, one_json)
        #case_data = json.load(open(json_file, 'r'))
        
    for one_dir in os.listdir(result_dir):
        if one_dir == "tmp":
            continue
        if one_dir.split("_")[1] not in sample_list:
            continue

        dir_path = os.path.join(result_dir, one_dir)
        file_path = os.path.join(dir_path, "super_enhancer", one_dir + "_SuperEnhancers.table.txt")
        if not os.path.exists(file_path):
            continue
        else:
            cmd = "cp %s %s;" % (file_path, super_enhancer_dir)
            Basic.run(cmd, wkdir=super_enhancer_dir)


def super_enhancer_peak_collection(super_enhancer_dir=None, output_bed=None):

    df_total = pd.DataFrame(columns=["REGION_ID", "CHROM", "START", "STOP", "enhancerRank", "isSuper", "sample_id"])
    for one_file in os.listdir(super_enhancer_dir):
        file_path = os.path.join(super_enhancer_dir, one_file)
        df_one = pd.read_csv(file_path, sep="\t", comment='#')
        df_one["sample_id"] = one_file.split("_")[1]
        if "enhancerRank" not in df_one.columns:
            continue
        df_one = df_one[["CHROM", "START", "STOP", "REGION_ID", "enhancerRank", "isSuper", "sample_id"]]
        df_total = df_total.append(df_one)
        print(df_one[0:10])
    df_total = df_total[["CHROM", "START", "STOP", "REGION_ID", "enhancerRank", "isSuper", "sample_id"]]
    df_total.to_csv(output_bed, index=False, sep="\t", header=None)

def cancer_normalize(phenotype=None):
    #cancer_dict_file = "/home/zhluo/Project/TF_enrichment/test/analysis/phenotype.txt"
    cancer_dict_file = phenotype
    handle = open(cancer_dict_file, "r")
    cancer_dict = {}
    for one_line in handle:
        array = one_line.strip("\n\r").split("\t")
        orignal_cancer = array[0]
        cancer_normalized = array[5]
        cancer_id = array[4]
        cancer_dict[orignal_cancer] = [cancer_normalized, cancer_id]    
    return(cancer_dict)

def creat_sample_info(super_enhancer_file=None, original_json_dir = None, output_file = "/home/zhluo/Project/TF_enrichment/test/analysis/all_sample_info.csv", phenotype="/home/zhluo/Project/TF_enrichment/cell_line/peak_clean/phenotype.txt"):
    #super_enhancer_file = "/home/zhluo/Project/TF_enrichment/test/analysis/all_super_enhancer.csv"
    #original_json_dir = "/home/mwshi/project/CEdb/metajson/"  
    
    result_dir = "/home/zhluo/Project/TF_enrichment/test/work_dir_all_data"
    GSM_list = []
    sample_info = []
    cancer_dict = cancer_normalize(phenotype=phenotype)
    
    for one_dir in os.listdir(result_dir):
        if one_dir == "tmp":
            continue
        dir_path = os.path.join(result_dir, one_dir)
        file_path = os.path.join(dir_path, "super_enhancer", one_dir + "_SuperEnhancers.table.txt")
        if not os.path.exists(file_path):
            continue
        GSM_list.append(one_dir.split("_")[1])
        
        #sample info
        GSM_json_file = os.path.join(original_json_dir, one_dir.split("_")[1] + ".json")
        if not os.path.exists(GSM_json_file):
            print("json file is missing!!!")
            continue
            
            #this part is cell line
            GSM_json_file = os.path.join("/home/mwshi/project/CEdb/cell_line/metajson/", one_dir.split("_")[1] + ".json")
            if not os.path.exists(GSM_json_file):
                print("json file is really missing!!!")
                continue          
        case_data = json.load(open(GSM_json_file, 'r'))
        
        
        
        tissue_type = case_data["layout"]
        data_type = case_data["layout"].lower()
        GSM_id = case_data["GSM"]
        GEO_id = case_data["GEO"]
        Biosample_type = "Primary tissue" if case_data["TUSSUE"] !=  "cell line" else "Cell line"
        Cancer_type = case_data["CANCER_TYPE"]
        has_input = "Y" if case_data["INPUT"] is not None else "NA"
        has_normal = "Y" if case_data["normal"] is not None else "NA"
        input_id = case_data["INPUT"] if case_data["INPUT"] != "" else "No input"
        normal_id = case_data["normal"] if case_data["normal"] != "" else "No adjacent"
        pass_coverage = "Y"
        pass_fastqc = "Y"
        
        handle = open(file_path, "r")
        super_enhancer_number = len(handle.readlines()) - 6
        handle.close()
        #print(cancer_dict)
        #print(case_data)
        if Cancer_type not in cancer_dict:
            print(GSM_json_file)
            continue
        #print(Cancer_type)
        sample_info.append({"GEO_id": GEO_id, "GSM_id": GSM_id, "Biosample_type": Biosample_type, "Cancer_type": Cancer_type, "cancer_abbr": cancer_dict[Cancer_type][0], "cancer_nor": cancer_dict[Cancer_type][1],  "has_input": has_input, "has_normal": has_normal, "pass_coverage": pass_coverage, "pass_fastqc": pass_fastqc, "super_enhancer_number": super_enhancer_number})
    df_sample = pd.DataFrame(sample_info)
    df_sample.to_csv(output_file, index=False)
    
def define_enhancer_gene(reference="hg38", super_enhancer_file=None, annotated_peaks_file=None):
    enahcner_bed_path = super_enhancer_file
    annotate_peaks = annotated_peaks_file #"/home/zhluo/Project/TF_enrichment/test/analysis/super_enhancer_annotate_peak.csv"
    cmd = "annotatePeaks.pl %s %s > %s" %(enahcner_bed_path, reference, annotate_peaks)
    Basic.run(cmd, wkdir= "/home/zhluo/Project/TF_enrichment/test/pbs/annotate_peak/")  
    
def filter_enahcner(annotated_peaks_file=None, filterd_file_path=None):
    #this function used to filter out some enahcners
    peak_anno_file = annotated_peaks_file
    #filterd_file = open("/home/zhluo/Project/TF_enrichment/test/analysis/super_enhancer_annotate_peak_remove_chrM.csv", "w")
    filterd_file = open(filterd_file_path, "w")
    handle = open(peak_anno_file, "r")
    header = handle.readline()
    filterd_file.write(header)
    for line in handle:
        chr_str = line.split("\t")[1]
        if not chr_str.startswith("c"):
            continue
        if chr_str == "chrM":
            continue
        filterd_file.write(line)
    filterd_file.close()
    handle.close()
    
def intersect_super_enhancer_cancer_snp(super_enhancer_file = None, cancer_snp_file=None, output_file=None):
    df = pd.read_csv(super_enhancer_file, sep="\t", skiprows=[0], names=["PeakID", "Chr", "Start", "End", "Strand", "Peak Score", "Focus Ratio", "Size Annotation", "Detailed Annotation", "Distance to TSS", "Nearest PromoterID", "Entrez ID", "Nearest Unigene", "Nearest Refseq",  "Nearest Ensembl", "Gene Name", "Gene Alias", "Gene Description", "Gene Type"])
    print(df[0:5])
    
    df = df[["Chr", "Start", "End", "PeakID"]]
    super_enhancer_tmp = super_enhancer_file + ".tmp"
    df.to_csv(super_enhancer_tmp, index=False, header=False, sep="\t")
    cmd = "bedtools intersect -a %s -b %s -wa -wb > %s" %(super_enhancer_tmp, cancer_snp_file, output_file)
    Basic.run(cmd)

def add_cancer_snp(SE_cancer_snp_file=None, super_enhancer_file=None, sample_info=None, output_file=None):
    #snp
    enhancer_snpid = {}
    cancer_snp_file = SE_cancer_snp_file
    snp_handle = open(cancer_snp_file, "r")
    for line in snp_handle:
        array = line.strip("\r\n").split("\t")
        enhancer_id = array[3]
        snp_id = array[8]
        if enhancer_id not in enhancer_snpid:
            enhancer_snpid[enhancer_id] = [snp_id]
        else:
            enhancer_snpid[enhancer_id].append(snp_id)        
    snp_handle.close()    
    
    
    peak_anno_file = super_enhancer_file
    peak_handle = open(peak_anno_file, "r")
    header = peak_handle.readline()
    name_dict = {}
    
    """
    #name
    json_file = open("/home/zhluo/Project/TF_enrichment/test/analysis/sample_info.json","r")
    mesh_id = {}
    for line in json_file:
        case_data = json.loads(line)
        mesh_id[case_data["GSM_id"]] = case_data["sample_id"]
    json_file.close()
    """
    df_sample = pd.read_csv(sample_info)
    cancer_type_id = dict(zip(df_sample["GSM_id"], df_sample["cancer_abbr"]))
    
    
    new_peak_handle = open(output_file, "w")
    header = header.strip("\n\r")
    header = header + "\tCancer_eQTL_numer\tCancer_eQTL_rsid\n"
    new_peak_handle.write(header)
    for line in peak_handle:
        array = line.strip("\r\n").split("\t")
        enhancer_id = array[0]
        #name_dict[enhancer_id]  = enhancer_id.split(".")[0].split("_")[2] + enhancer_id.split(".")[1].split("_")[2]
        new_name  = enhancer_id.split(".")[0].split("_")[2] +"_" + enhancer_id.split(".")[1].split("_")[2]
        #array[0] = new_name
        GSM_id = enhancer_id.split(".")[0].split("_")[2]
        if GSM_id not in cancer_type_id:
            print(GSM_id)
            continue
        array[0] = GSM_id + "_" + cancer_type_id[enhancer_id.split(".")[0].split("_")[2]]+"_" + enhancer_id.split(".")[1].split("_")[2]
        if enhancer_id in enhancer_snpid:
            snps = enhancer_snpid[enhancer_id]
            snps = list(set(snps))
            snp_str = ",".join(snps)
            num = len(snps)
        else:
            num = 0
            snp_str = ""                    
        array.append(str(num))
        array.append(snp_str)
        new_peak_handle.write("\t".join(array) + "\n") 
    peak_handle.close()
    new_peak_handle.close()
    
#GTEx
def get_lookup_file():
    #get lookup file
    snp_id_df = "/home/zhluo/Project/TF_enrichment/test/analysis/GETx_eQTL/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt"
    df_snp = pd.read_csv(snp_id_df, sep="\t")
    df_snp = df_snp[["variant_id", "rs_id_dbSNP151_GRCh38p7"]]
    #df_snp["rs_id_dbSNP147_GRCh37p13"].str.startswith("rs")
    df_snp['has_id'] = list(map(lambda x: x.startswith('rs'), df_snp["rs_id_dbSNP151_GRCh38p7"])) 
    df_snp = df_snp[df_snp["has_id"] == True]
    df_snp = df_snp[["variant_id", "rs_id_dbSNP151_GRCh38p7"]]
    df_snp.to_csv("/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/GTEx_snp_name_and_id.txt", sep="\t", index=False)
    print(df_snp[0:10]) 

def gtex_eQTL():
    df_snp = pd.read_csv("/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/GTEx_snp_name_and_id.txt", sep="\t") 
    #snp_dict = dict(zip(df_snp["variant_id"], df_snp["rs_id_dbSNP147_GRCh37p13"]))
    #exit(1)
    human_gene_table = "/home/zhluo/Project/TF_enrichment/test/analysis/GETx_eQTL/human_ensembl_name.txt"
    df_ensembl = pd.read_csv(human_gene_table, sep="\t")
    df_ensembl.columns = ["gene_id", "gene_symbol"]
    #print(df_ensembl[0:10])
    #exit(1)
    
    
    output_dir = "/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/bed_GTEX"
    Basic.mkdir(output_dir)
    gtex_dir = "/home/zhluo/Project/TF_enrichment/test/analysis/GETx_eQTL/GTEx_Analysis_v8_eQTL"
    for one_file in os.listdir(gtex_dir):
        if  ".egenes.txt" in one_file:
            continue
        file_path = os.path.join(gtex_dir, one_file)
        df = pd.read_csv(file_path, sep="\t")
        df = df[["variant_id", "gene_id", "pval_nominal"]]
        #df[['First','Last']] = df.Name.str.split(" ",expand=True,)
        
        df_merge = pd.merge(df, df_snp, how="inner", left_on="variant_id", right_on = "variant_id")
        df_merge[["chr", "pos", "maf", "mif", "ver"]] = df_merge.variant_id.str.split("_",expand=True)
        df_merge["chr"] = df_merge["chr"]
        df_merge["start"] = df_merge["pos"]
        df_merge["end"] = df_merge["pos"]
        df_merge["tissue"] = one_file.split(".")[0]
        df_merge['gene_id'] = df_merge['gene_id'].map(lambda x: x.split('.')[0])
        df_merge = pd.merge(df_merge, df_ensembl, how="inner", left_on="gene_id", right_on = "gene_id")
        df_merge = df_merge[["chr", "start", "end", "gene_symbol", "rs_id_dbSNP151_GRCh38p7", "pval_nominal", "tissue"]]   
        df_merge.to_csv(os.path.join(output_dir, one_file), sep="\t", index=False, header=False)

def total_gtex_eQTL(gtex_file_dir="/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/bed_GTEX", filterd_file_path=None, output=None):
    total_gtex_file = "/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/gtex_total_snp.bed"
    cmd = "cat %s/* >%s;" %(gtex_file_dir, total_gtex_file)
    Basic.run(cmd)
    
    super_enhancer_tmp = filterd_file_path + ".tmp"
    cmd = "bedtools intersect -a %s -b %s -wa -wb >%s"%(super_enhancer_tmp, total_gtex_file, output)
    Basic.run(cmd)

def add_GTEx_snp(SE_cancer_snp_file=None, super_enhancer_file=None, sample_info=None, output_file=None):
    #snp
    enhancer_snpid = {}
    cancer_snp_file = SE_cancer_snp_file
    snp_handle = open(cancer_snp_file, "r")
    for line in snp_handle:
        array = line.strip("\r\n").split("\t")
        enhancer_id = array[3]
        snp_id = array[8]
        if enhancer_id not in enhancer_snpid:
            enhancer_snpid[enhancer_id] = [snp_id]
        else:
            enhancer_snpid[enhancer_id].append(snp_id)        
    snp_handle.close()    
    
    
    peak_anno_file = super_enhancer_file
    peak_handle = open(peak_anno_file, "r")
    header = peak_handle.readline()
    name_dict = {}
    
    """
    #name
    json_file = open("/home/zhluo/Project/TF_enrichment/test/analysis/sample_info.json","r")
    mesh_id = {}
    for line in json_file:
        case_data = json.loads(line)
        mesh_id[case_data["GSM_id"]] = case_data["sample_id"]
    json_file.close()
    """
    df_sample = pd.read_csv(sample_info)
    cancer_type_id = dict(zip(df_sample["GSM_id"], df_sample["cancer_abbr"]))
    
    
    new_peak_handle = open(output_file, "w")
    header = header.strip("\n\r")
    header = header + "\tGTEx_eQTL_numer\tGTEx_eQTL_rsid\n"
    new_peak_handle.write(header)
    for line in peak_handle:
        array = line.strip("\r\n").split("\t")
        enhancer_id = array[0]
        #name_dict[enhancer_id]  = enhancer_id.split(".")[0].split("_")[2] + enhancer_id.split(".")[1].split("_")[2]
        new_name  = enhancer_id.split(".")[0].split("_")[2] +"_" + enhancer_id.split(".")[1].split("_")[2]
        #array[0] = new_name
        GSM_id = enhancer_id.split(".")[0].split("_")[2]
        if GSM_id not in cancer_type_id:
            print(GSM_id)
            continue
        array[0] = GSM_id + "_" + cancer_type_id[enhancer_id.split(".")[0].split("_")[2]]+"_" + enhancer_id.split(".")[1].split("_")[2]
        if enhancer_id in enhancer_snpid:
            snps = enhancer_snpid[enhancer_id]
            snps = list(set(snps))
            snp_str = ",".join(snps)
            num = len(snps)
        else:
            num = 0
            snp_str = ""                    
        array.append(str(num))
        array.append(snp_str)
        new_peak_handle.write("\t".join(array) + "\n") 
    peak_handle.close()
    new_peak_handle.close()
    
def get_GWAS_snp():
    """
    the full file is downloaded from: https://www.ebi.ac.uk/gwas/api/search/downloads/full
    """
    gwas_file  = "/home/zhluo/Project/TF_enrichment/test/analysis/gwas_snp/full"
    df = pd.read_csv(gwas_file, sep="\t")
    df_snp = df[['SNPS', 'MAPPED_GENE', 'DISEASE/TRAIT', 'CONTEXT', 'P-VALUE', 'PUBMEDID']]
    df_snp['has_id'] = list(map(lambda x: x.startswith('rs'), df_snp["SNPS"])) 
    df_snp = df_snp[df_snp["has_id"] == True]
    #df_snp['has_more'] = list(map(lambda x: len(x) < 17 and "-" not in x and "_" not in x and "d" not in x and ";" not in x and bool(re.match(r"rs[0-9]+$", x)), df_snp["SNPS"])) 
    df_snp['has_more'] = list(map(lambda x: bool(re.match(r"^rs[0-9]+$", x)), df_snp["SNPS"])) 
    df_snp = df_snp[df_snp["has_more"] == True]
    df_snp.to_csv("/home/zhluo/Project/TF_enrichment/test/analysis/gwas_snp/gwas_snp.txt", sep="\t", index=False)
    
def select_cancer_gwas_snp():
    gwas_snp_file = "/home/zhluo/Project/TF_enrichment/test/analysis/gwas_snp/gwas_snp.txt"
    gwas_cancer_phenotype = "/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/gwas_cancer_phenotype.txt"
    
    phenotype_df = pd.read_csv(gwas_cancer_phenotype, names=["rank1", "rank2", "disease", "mesh_id", "tree_code"], sep="\t")
    print(phenotype_df[0:10])
    
    gwas_df = pd.read_csv(gwas_snp_file, sep="\t")
    print(gwas_df[0:10])
    
    gwas_selected = gwas_df[gwas_df['DISEASE/TRAIT'].isin(phenotype_df["disease"].to_list())]
    print(gwas_selected[0:30])
    gwas_selected.to_csv("/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/cancer_gwas_snp.txt", sep="\t", index=False)

def gwas_intersect(filterd_file_path=None, output=None):
    super_enhancer_tmp = filterd_file_path + ".tmp"
    gwas_bed = "/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/cancer_gwas_snp_finall.bed"
    cmd = "bedtools intersect -a %s -b %s -wa -wb >%s" %(super_enhancer_tmp, gwas_bed, output)
    Basic.run(cmd)
    
def count_GWAS_SNP_number(SE_cancer_snp_file=None, super_enhancer_file=None, sample_info=None, output_file=None):
    #snp
    enhancer_snpid = {}
    cancer_snp_file = SE_cancer_snp_file
    snp_handle = open(cancer_snp_file, "r")
    for line in snp_handle:
        array = line.strip("\r\n").split("\t")
        enhancer_id = array[3]
        snp_id = array[8]
        if enhancer_id not in enhancer_snpid:
            enhancer_snpid[enhancer_id] = [snp_id]
        else:
            enhancer_snpid[enhancer_id].append(snp_id)        
    snp_handle.close()    
    
    
    peak_anno_file = super_enhancer_file
    peak_handle = open(peak_anno_file, "r")
    header = peak_handle.readline()
    name_dict = {}
    
    """
    #name
    json_file = open("/home/zhluo/Project/TF_enrichment/test/analysis/sample_info.json","r")
    mesh_id = {}
    for line in json_file:
        case_data = json.loads(line)
        mesh_id[case_data["GSM_id"]] = case_data["sample_id"]
    json_file.close()
    """
    df_sample = pd.read_csv(sample_info)
    cancer_type_id = dict(zip(df_sample["GSM_id"], df_sample["cancer_abbr"]))
    
    
    new_peak_handle = open(output_file, "w")
    header = header.strip("\n\r")
    header = header + "\tGWAS_snp_numer\tGWAS_snp_rsid\n"
    new_peak_handle.write(header)
    for line in peak_handle:
        array = line.strip("\r\n").split("\t")
        enhancer_id = array[0]
        #name_dict[enhancer_id]  = enhancer_id.split(".")[0].split("_")[2] + enhancer_id.split(".")[1].split("_")[2]
        new_name  = enhancer_id.split(".")[0].split("_")[2] +"_" + enhancer_id.split(".")[1].split("_")[2]
        #array[0] = new_name
        GSM_id = enhancer_id.split(".")[0].split("_")[2]
        if GSM_id not in cancer_type_id:
            print(GSM_id)
            continue
        array[0] = GSM_id + "_" +  cancer_type_id[enhancer_id.split(".")[0].split("_")[2]]+"_" + enhancer_id.split(".")[1].split("_")[2]
        if enhancer_id in enhancer_snpid:
            snps = enhancer_snpid[enhancer_id]
            snps = list(set(snps))
            snp_str = ",".join(snps)
            num = len(snps)
        else:
            num = 0
            snp_str = ""                    
        array.append(str(num))
        array.append(snp_str)
        new_peak_handle.write("\t".join(array) + "\n") 
    peak_handle.close()
    new_peak_handle.close()
    
def merge_3_table(cancer_eqtl=None, gtex_eqtl=None, gwas_snp=None, output=None):
    cancer_eqtl = cancer_eqtl
    gtex_eqtl = gtex_eqtl
    gwas_snp = gwas_snp
    
    
    df_cancer_eqtl = pd.read_csv(cancer_eqtl, sep="\t")
    gtex_eqtl = pd.read_csv(gtex_eqtl, sep="\t")
    gwas_snp = pd.read_csv(gwas_snp, sep="\t")
    
    df_cancer_eqtl["GTEx_eQTL_numer"] = gtex_eqtl["GTEx_eQTL_numer"]
    df_cancer_eqtl["GTEx_eQTL_rsid"] = gtex_eqtl["GTEx_eQTL_rsid"]
    df_cancer_eqtl["GWAS_snp_numer"] = gwas_snp["GWAS_snp_numer"]
    df_cancer_eqtl["GWAS_snp_rsid"] = gwas_snp["GWAS_snp_rsid"]
    
    df_cancer_eqtl.to_csv(output, sep="\t", index=False)
    
def redundent_snp():
    gtex_snp = "/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/gtex_total_snp.bed"
    handle = open(gtex_snp, "r")
    
    
    snp_dict = {}
    for line in handle:
        array = line.strip("\n\r").split("\t")
        loc = ":".join(array[0:3] + [array[4]])
        qtl = "\t".join([array[3]] + array[5:])
        if loc in snp_dict:
            snp_dict[loc].append(qtl)
        else:
            snp_dict[loc] = [qtl]
    handle.close()
        
    handle = open("/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/primary_tissue_used_gtex_snp_unique_redundent.txt", "w")
    for one_key in snp_dict.keys():
        qtls = snp_dict[one_key]
        handle.write("%s:%s\n" %(one_key, ",".join(qtls)))
    handle.close()

def run():
    result_prefix = "/home/zhluo/Project/TF_enrichment/test/result_20210912"
    cancer_total_eqtl = "/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/cancer_eQTL_total_merged.bed"
    
    #primary
    json_dir_path = "/home/mwshi/project/CEdb/metajson/"
    orignal_output_dir = "/home/zhluo/Project/TF_enrichment/test/work_dir_all_data"
    super_enhancer_dir = os.path.join(result_prefix, "super_enhancer_peaks")
    Basic.mkdir(super_enhancer_dir)
    analysis_dir = os.path.join(result_prefix, "analysis")
    Basic.mkdir(analysis_dir)
    super_enhancer_bed = os.path.join(analysis_dir, "all_super_enhancer.bed")
    all_sample_info = os.path.join(analysis_dir ,"all_sample_info.csv")
    phenotype = "/home/zhluo/Project/TF_enrichment/cell_line/peak_clean/phenotype.txt"
    annotated_peaks_file = os.path.join(analysis_dir ,"super_enhancer_annotate_peak.csv")
    filterd_file_path = os.path.join(analysis_dir ,"super_enhancer_annotate_peak_remove_chrM.csv")
    #snp_map_file = os.path.join(analysis_dir ,"hg19_snp_to_grch38.txt")
    super_enahcner_cancer_snp_file = os.path.join(analysis_dir, "super_enhancer_cancer_eQTL.txt")
    super_enhancer_cancer_eqtl_summary = os.path.join(analysis_dir, "super_enhancer_annotate_peak_remove_chrM_add_cancerEQTL.csv")
    super_enhancer_gtex_eqtl = os.path.join(analysis_dir, "super_enhancer_GTEx_eQTL.txt")
    super_enhancer_gtex_eqtl_summary = os.path.join(analysis_dir, "super_enhancer_annotate_peak_remove_chrM_add_GTExeQTL.csv")
    super_enhancer_gwas_snp = os.path.join(analysis_dir, "super_enhancer_GWAS_SNP.txt")
    super_enhancer_gwas_snp_summary = os.path.join(analysis_dir, "super_enhancer_annotate_peak_remove_chrM_add_GWAS_SNP.csv")
    merged_table = os.path.join(analysis_dir, "super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv")
    
    ##1. collect bed
    #super_enhancer_data_collection(meta_data_dir = json_dir_path, result_dir = orignal_output_dir, super_enhancer_dir=super_enhancer_dir)
    ##2. collect peak
    #super_enhancer_peak_collection(super_enhancer_dir=super_enhancer_dir, output_bed=super_enhancer_bed)
    ##3. case info
    #creat_sample_info(super_enhancer_file=super_enhancer_bed, original_json_dir = json_dir_path, output_file = all_sample_info, phenotype=phenotype)
    ##4. annotate peaks
    #define_enhancer_gene(reference="hg38", super_enhancer_file=super_enhancer_bed, annotated_peaks_file=annotated_peaks_file)
    ##5. filter enhancers
    #filter_enahcner(annotated_peaks_file=annotated_peaks_file, filterd_file_path=filterd_file_path)
    ##convert snp hg19 to hg38, step2_snp_convert.py
    ##6. add cancer eQTL
    #intersect_super_enhancer_cancer_snp(super_enhancer_file = filterd_file_path, cancer_snp_file=cancer_total_eqtl, output_file=super_enahcner_cancer_snp_file)
    ##7. summary cancer eQTL
    #add_cancer_snp(SE_cancer_snp_file=super_enahcner_cancer_snp_file, super_enhancer_file=filterd_file_path, sample_info=all_sample_info, output_file=super_enhancer_cancer_eqtl_summary)
    ##8. GTEx
    #get_lookup_file()
    #gtex_eQTL()
    #total_gtex_eQTL(gtex_file_dir="/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/bed_GTEX", filterd_file_path=filterd_file_path, output=super_enhancer_gtex_eqtl)
    #add_GTEx_snp(SE_cancer_snp_file=super_enhancer_gtex_eqtl, super_enhancer_file=filterd_file_path, sample_info=all_sample_info, output_file=super_enhancer_gtex_eqtl_summary)
    ##9. add GWAS snp
    #select_cancer_gwas_snp()
    #gwas_intersect(filterd_file_path=filterd_file_path, output=super_enhancer_gwas_snp)
    #count_GWAS_SNP_number(SE_cancer_snp_file=super_enhancer_gwas_snp, super_enhancer_file=filterd_file_path, sample_info=all_sample_info, output_file=super_enhancer_gwas_snp_summary)
    ##10. merge
    #merge_3_table(cancer_eqtl=super_enhancer_cancer_eqtl_summary, gtex_eqtl=super_enhancer_gtex_eqtl_summary, gwas_snp=super_enhancer_gwas_snp_summary, output=merged_table)
    ##database table
    redundent_snp()
    
    ###################################################cell line###################################################################################
    #cell line
    json_dir_path = "/home/mwshi/project/CEdb/cell_line/metajson/"
    orignal_output_dir = "/home/zhluo/Project/TF_enrichment/test/work_dir_all_data"
    super_enhancer_dir = os.path.join(result_prefix, "cell_line_super_enhancer_peaks")
    Basic.mkdir(super_enhancer_dir)
    analysis_dir = os.path.join(result_prefix, "analysis")
    Basic.mkdir(analysis_dir)
    super_enhancer_bed = os.path.join(analysis_dir, "cell_line_all_super_enhancer.bed")
    all_sample_info = os.path.join(analysis_dir ,"cell_line_all_sample_info.csv")
    phenotype = "/home/zhluo/Project/TF_enrichment/cell_line/peak_clean/cell_line_phenotype.txt"
    annotated_peaks_file = os.path.join(analysis_dir ,"cell_line_super_enhancer_annotate_peak.csv")
    filterd_file_path = os.path.join(analysis_dir ,"cell_line_super_enhancer_annotate_peak_remove_chrM.csv")
    super_enahcner_cancer_snp_file = os.path.join(analysis_dir, "cell_line_super_enhancer_cancer_eQTL.txt")
    super_enhancer_cancer_eqtl_summary = os.path.join(analysis_dir, "cell_line_super_enhancer_annotate_peak_remove_chrM_add_cancerEQTL.csv")
    super_enhancer_gtex_eqtl = os.path.join(analysis_dir, "cell_line_super_enhancer_GTEx_eQTL.txt")
    super_enhancer_gtex_eqtl_summary = os.path.join(analysis_dir, "cell_line_super_enhancer_annotate_peak_remove_chrM_add_GTExeQTL.csv")
    super_enhancer_gwas_snp = os.path.join(analysis_dir, "cell_line_super_enhancer_GWAS_SNP.txt")
    super_enhancer_gwas_snp_summary = os.path.join(analysis_dir, "cell_line_super_enhancer_annotate_peak_remove_chrM_add_GWAS_SNP.csv")
    merged_table = os.path.join(analysis_dir, "cell_line_super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv")
    ##1. collect bed
    #super_enhancer_data_collection(meta_data_dir = json_dir_path, result_dir = orignal_output_dir, super_enhancer_dir=super_enhancer_dir)
    ##2. collect peak
    #super_enhancer_peak_collection(super_enhancer_dir=super_enhancer_dir, output_bed=super_enhancer_bed)
    ##3. case info
    #creat_sample_info(super_enhancer_file=super_enhancer_bed, original_json_dir = json_dir_path, output_file = all_sample_info, phenotype=phenotype)
    ##4. annotate peaks
    #define_enhancer_gene(reference="hg38", super_enhancer_file=super_enhancer_bed, annotated_peaks_file=annotated_peaks_file)
    ##5. filter enhancers
    #filter_enahcner(annotated_peaks_file=annotated_peaks_file, filterd_file_path=filterd_file_path)
    ##6. add cancer eQTL
    #intersect_super_enhancer_cancer_snp(super_enhancer_file = filterd_file_path, cancer_snp_file=cancer_total_eqtl, output_file=super_enahcner_cancer_snp_file)
    ##7. summary cancer eQTL
    #add_cancer_snp(SE_cancer_snp_file=super_enahcner_cancer_snp_file, super_enhancer_file=filterd_file_path, sample_info=all_sample_info, output_file=super_enhancer_cancer_eqtl_summary)
    ##GTEx
    #total_gtex_eQTL(gtex_file_dir="/home/zhluo/Project/TF_enrichment/test/result_20210912/analysis/bed_GTEX", filterd_file_path=filterd_file_path, output=super_enhancer_gtex_eqtl)
    #add_GTEx_snp(SE_cancer_snp_file=super_enhancer_gtex_eqtl, super_enhancer_file=filterd_file_path, sample_info=all_sample_info, output_file=super_enhancer_gtex_eqtl_summary)
    ##9. add GWAS snp
    #gwas_intersect(filterd_file_path=filterd_file_path, output=super_enhancer_gwas_snp)
    #count_GWAS_SNP_number(SE_cancer_snp_file=super_enhancer_gwas_snp, super_enhancer_file=filterd_file_path, sample_info=all_sample_info, output_file=super_enhancer_gwas_snp_summary)
    #merge_3_table(cancer_eqtl=super_enhancer_cancer_eqtl_summary, gtex_eqtl=super_enhancer_gtex_eqtl_summary, gwas_snp=super_enhancer_gwas_snp_summary, output=merged_table)
if __name__ == "__main__":
    run()
    """
    gene="VEGFA"
    grep $gene super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv   |wc -l
    grep $gene super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv  |grep "COAD" |wc -l
    grep $gene cell_line_super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv  |wc -l
    grep $gene cell_line_super_enhancer_annotate_peak_remove_chrM_add_gwas_gtex_cancer_snp.csv  |grep "COAD" |wc -l
    
    """