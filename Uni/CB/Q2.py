# Python code to extract the sample name from these files ignoring any files which do not match the format given below
# 1. Written lane number
# 2. Barcode
# 3. Sample name
# 4. Numeric lane number (starting with L)
# 5. Read number (R1/2/3/4)
# 6. File extension

import re


def extract_sample_name(file_name):
    pattern = r'^lane\d+_([A-Z\d]+)_([A-Za-z\d_]+)_L\d+_([R]\d+)\.fastq\.gz$'
    match = re.match(pattern, file_name)

    if match:
        return match.group(2)


file_names = [
    "lane1_NewCode_L001_R1.fastq.gz",
    "lane1_NoIndex_L001_R1.fastq.gz",
    "lane1_NoIndex_L001_R2.fastq.gz",
    "pipeline_processing_output.log",
    "lane7027_ACTGAT_JH25_L001_R1.fastq.gz",
    "lane7027_ACTTGA_E30_1_2_Hap4_24h_L001_R1.fastq.gz",
    "lane7027_AGTTCC_JH14_L001_R1.fastq.gz",
    "lane7027_CGGAAT_JH37_L001_R1.fastq.gz",
    "lane7027_GCCAAT_E30_1_2l_Hap4_log_L001_R1.fastq.gz",
    "lane7127_GGCTAC_E30_1_4_Hap4_48h_L001_R1.fastq.gz",
    "lane8127_GCCAAT_S30_1_2l_Hap4_log_L001_R1.fastq.gz"
]

for file_name in file_names:
    sample_name = extract_sample_name(file_name)
    if sample_name is not None:
        print(f"File: {file_name}, Sample Name: {sample_name}")
