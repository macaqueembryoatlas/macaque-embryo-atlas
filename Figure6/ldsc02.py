import os

input_dir = '/data/work/blood_archr_orignal/mac2hg19/'

cell_types = []

for filename in os.listdir(input_dir):
    if filename.endswith('.bed'):
        celltype = filename.replace('.bed', '')
        cell_types.append(celltype)


# 定义染色体范围
chromosomes = range(1, 23)

# 输出文件名
output_file = "/data/work/02.LDScore_blood.sh"

# 打开文件并写入命令
with open(output_file, "w") as f:
    for cell_type in cell_types:
        for chrom in chromosomes:
            cmd = (
                f"python '/data/users/changan/changan_eaa4dc72964b415089b45ea85dfe72d1/output/NB2025101917283071341064/storage/storage/LDSC_GWAS/ldsc-master/ldsc.py' "
                f"--l2 "
                f"--bfile '/data/users/changan/changan_eaa4dc72964b415089b45ea85dfe72d1/output/NB2025101917283071341064/storage/storage/LDSC_GWAS/LDSC_reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}' "
                f"--ld-wind-cm 1 "
                f"--annot /data/work/01.AnnotFile/{cell_type}.{chrom}.annot.gz "
                f"--thin-annot "
                f"--out /data/work/02.LDScore/{cell_type}.{chrom} "
                f"--print-snps '/data/users/changan/changan_eaa4dc72964b415089b45ea85dfe72d1/online/storage/LDSC_GWAS/LDSC_reference/hapmap3_snps/hm.{chrom}.snp';\n"
            )
            f.write(cmd)

print(f"Shell script {output_file} has been generated.")
