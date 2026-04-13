import os

input_dir = '/data/work/blood_archr_orignal/mac2hg19/'

cell_types = []

for filename in os.listdir(input_dir):
    if filename.endswith('.bed'):
        celltype = filename.replace('.bed', '')
        cell_types.append(celltype)

import os

input_dir = '/data/work/blood_archr_orignal/mac2hg19/'

cell_types = []

for filename in os.listdir(input_dir):
    if filename.endswith('.bed'):
        celltype = filename.replace('.bed', '')
        cell_types.append(celltype)

        
sh_file = "/data/work/01.AnnotFile_blood.sh"

# 打开文件准备写入
with open(sh_file, "w") as f:
    # 遍历细胞类型
    for cell_type in cell_types:

        # 遍历染色体1到22
        for chrom in range(1, 23):
            # 生成命令
            bed_file = f"/data/users/changan/changan_8e27def345ec42198941baae9c80c122/online/blood_archr_orignal/mac2hg19/{cell_type}.bed"#一定要改online后面，改为需要读入的dapeak文件夹
            bim_file = f"/data/users/changan/changan_eaa4dc72964b415089b45ea85dfe72d1/output/NB2025101917283071341064/storage/storage/LDSC_GWAS/LDSC_reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}.bim"
            annot_file = f"/data/work/01.AnnotFile/{cell_type}.{chrom}.annot.gz"  # 修改为根据染色体号生成文件名
            
            # 将命令写入shell脚本
            command = f"python2 /data/users/changan/changan_891358097e434bc6a3f4719ea5062ea4/online/make_annot.py --bed-file {bed_file} --bimfile {bim_file} --annot-file {annot_file}\n"
            f.write(command)

print(f"Shell脚本 {sh_file} 已生成。")
