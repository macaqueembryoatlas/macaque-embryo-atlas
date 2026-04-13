import os

input_dir = '/data/work/blood_archr_orignal/mac2hg19/'

cell_types = []

for filename in os.listdir(input_dir):
    if filename.endswith('.bed'):
        celltype = filename.replace('.bed', '')
        cell_types.append(celltype)

# 定义文件路径模板
base_path = "/data/work/02.LDScore/"
background_file = f"{base_path}all_hg19."
output_file = "/data/work/blood_cts.ldcts"

# 打开输出文件
with open(output_file, 'w') as f:
    for cell_type in cell_types:
        # 构造LD得分文件路径
        ldscore_file = f'{base_path}{cell_type}.'
        # 生成cts.ldcts中的每一行，添加背景文件路径
        f.write(f'{cell_type}\t{ldscore_file},{background_file}\n')

print(f"blood_cts.ldcts has been generated with the specified cell types and background file.")

# 导入必要的模块
import os

# 输入和输出文件路径
input_file = "/data/work/ldsc_disease_blood.txt"  # 原始文件，包含疾病名称和统计文件路径,(这里需要改)
output_script = "/data/work/03.DiseaseEnrichment_blood.sh"  # 生成的 Shell 脚本文件(这里需要改)

# 配置路径模板

ldsc_script = "/data/users/changan/changan_eaa4dc72964b415089b45ea85dfe72d1/online/storage/LDSC_GWAS/ldsc-master/ldsc.py"
ref_ld_chr = "/data/users/changan/changan_eaa4dc72964b415089b45ea85dfe72d1/output/NB2025101917283071341064/storage/storage/LDSC_GWAS/LDSC_reference/1000G_EUR_Phase3_baseline/baseline."
output_dir = "/data/work/03.DiseaseEnrichment/"
ref_ld_chr_cts = "/data/work/blood_cts.ldcts"#改这里
weights_ld_chr = "/data/users/changan/changan_eaa4dc72964b415089b45ea85dfe72d1/output/NB2025101917283071341064/storage/storage/LDSC_GWAS/LDSC_reference/weights_hm3_no_hla/weights."

# 读取输入文件并生成命令
with open(input_file, "r") as infile, open(output_script, "w") as outfile:
    for line in infile:
        # 跳过空行
        if not line.strip():
            continue
        
        # 分割行中的两个字段
        disease, sumstats_path = line.strip().split()
        
        # 构造命令
        command = (
            f"python2 {ldsc_script} "
            f"--h2-cts {sumstats_path} "
            f"--ref-ld-chr {ref_ld_chr} "
            f"--out {output_dir}{disease} "
            f"--ref-ld-chr-cts {ref_ld_chr_cts} "
            f"--w-ld-chr {weights_ld_chr};"
        )
        
        # 写入到输出脚本
        outfile.write(command + "\n")

print(f"Shell script '{output_script}' has been successfully generated.")
