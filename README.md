

# create_plant_cisTarget_databases

Documenting my process of building a cistarge database for pyscenic analysis of other species



## Install create_cisTarget_databases
For more detiled instructions see [installation](https://github.com/aertslab/create_cisTarget_databases?tab=readme-ov-file#installation).

## Prepare input files

Extract promoter sequence.

### Download TF information
For plants, we can download TF and motif information from [PlantTFDB](https://planttfdb.gao-lab.org/), [JASPAR](https://jaspar.elixir.no/) and [CIS-BP Database](https://cisbp.ccbr.utoronto.ca/).

For example, I download *Solanum tuberosum* TF from CIS-BP Database. Then unzip file.

```shell
unzip Solanum_tuberosum_2025_07_04_7:10_am.zip
logos_all_motifs
pwms_all_motifs
README.txt
SignalIntensities.txt
TF_Information_all_motifs_plus.txt
TF_Information_all_motifs.txt
TF_Information.txt
```

### create `.cb` files

```shell
cut -f4,6 TF_Information.txt | grep PGSC | awk '{print $1"_"$2}' | grep 3.00 > motifs_id.txt
python generate_cb_files.py
# Generated: motif_dir/M02635_3.00_PGSC0003DMG401010558.cb
# Generated: motif_dir/M01407_3.00_PGSC0003DMG401031196.cb
# Generated: motif_dir/M07782_3.00_PGSC0003DMG402006935.cb
# ...
# Generated: motif_dir/M02630_3.00_PGSC0003DMG402007388.cb
# Generated: motif_dir/M02649_3.00_PGSC0003DMG402028822.cb
# 
# Summary:
# Successfully processed: 935 files
# Missing PWM files: 0
# All .cb files generated in motif_dir/
```



### Extract promoter sequenc

Input file:

- gene.gtf

- genome.fa

Output file:

- Stu_dm4_gene_3kpromoter.fasta
- Stu_dm4_CDS_3kpromoter.fasta

```R
library(rtracklayer)
library(tidyverse)
library(Biostrings)

gtf <- import("/lustre/home/jianghao/01_database/01_ref/02_dm4/dm4.strand_fixed.gtf",format = "gtf") %>% 
  as.data.frame()

# get gene promoters regions
genes <- gtf %>% 
  dplyr::filter(type == "gene" & gene_name != "NA") %>% 
  dplyr::select(seqnames, start, end, strand, gene_name) %>% 
  dplyr::mutate(start2 = ifelse(strand == "+",start - 3000,end + 1),
                end2 = ifelse(strand == "+",start - 1,end + 3000))

head(genes)
#     seqnames  start    end strand            gene_name start2   end2
# 1 ST4.03ch00  63411  66816      + PGSC0003DMG400039401  60411  63410
# 2 ST4.03ch00  70051  73227      + PGSC0003DMG400013996  67051  70050
# 3 ST4.03ch00 164907 174068      - PGSC0003DMG400044441 174069 177068
# 4 ST4.03ch00 174833 179264      + PGSC0003DMG400013995 171833 174832
# 5 ST4.03ch00 192677 194558      - PGSC0003DMG400038530 194559 197558
# 6 ST4.03ch00 275470 276001      - PGSC0003DMG400046907 276002 279001

genome <- readDNAStringSet(filepath = "/lustre/home/jianghao/01_database/01_ref/02_dm4/dm4.fa")

library(Biostrings)

########## Extract 3k bp sequence upstream of TSS #############
seq_list <- lapply(1:nrow(genes), function(x) {
  tmp <- genes[x, ]

  # 尝试提取序列（使用 tryCatch 防止越界等错误）
  out <- tryCatch({
    seq <- genome[[as.character(tmp$seqnames)]][tmp$start2:tmp$end2]
    
    # 如果是负链，取反向互补
    if (tmp$strand == "-") {
      seq <- reverseComplement(seq)
    }
    return(seq)
  }, error = function(e) {
    # 提取失败返回空 DNAString
    return(DNAString())
  })

  return(out)
})

# 转换为 DNAStringSet 对象
seq_list <- DNAStringSet(seq_list)

# 添加名称
names(seq_list) <- genes$gene_name

writeXStringSet(fasta_filtered,filepath = "Stu_dm4_gene_3kpromoter.fasta",format="fasta")

########## Extract 3k bp sequence upstream of CDS #############
cds_seqs <- lapply(1:nrow(cds_promoter), function(i) {
  tmp <- cds_promoter[i, ]
  
  chr <- as.character(tmp$seqnames)
  
  # 确保染色体存在
  if (!chr %in% names(genome)) return(DNAString())
  
  # 获取上下游坐标（确保 start 不小于 1，end 不超过 chr 长度）
  chr_seq <- genome[[chr]]
  start <- max(tmp$start, 1)
  end <- min(tmp$end, length(chr_seq))

  # 尝试截取序列
  seq <- tryCatch({
    subseq(chr_seq, start = start, end = end)
  }, error = function(e) {
    DNAString()
  })

  # 负链要 reverseComplement
  if (tmp$strand == "-") {
    seq <- reverseComplement(seq)
  }

  return(seq)
})

# 封装成 DNAStringSet
cds_seqs <- DNAStringSet(cds_seqs)

# 设置名称为 gene_id 或 gene_name
names(cds_seqs) <- cds_promoter$gene_id  # or $gene_name

writeXStringSet(fasta_filtered,filepath = "Stu_dm4_CDS_3kpromoter.fasta",format="fasta")
```

### create_cisTarget_databases

```shell
python create_cisTarget_databases/create_cistarget_motif_databases.py \
	-f Stu_dm4_3kpromoter.fasta \
	-M cis_bp/motif_dir/ -m cis_bp/motifs_id.filtered.txt \
	-t 15 -o Stu_dm4_gene
```

Output files:

- Stu_dm4_gene.motifs_vs_regions.scores.feather
- Stu_dm4_gene.regions_vs_motifs.rankings.feather
- Stu_dm4_gene.regions_vs_motifs.scores.feather



### Create `motif2TF.tbl`

````R
library(dplyr)
motif2TF <- read.table("/lustre/home/jianghao/02_workspace/08_spatial_omics/18_GRN/pyscenic_databases/cis_bp/TF_Information_all_motifs_plus.txt",
    header = TRUE, sep = "\t")
head(motif2TF)
````

![image-20250707230053229](D:\06_GitHub\create_plant_cisTarget_databases\image-20250707230053229.png)

````R

# Assuming motif2TF is a three-column data frame with the information of motifs and their corresponding TFs like this:
# motif TF source
# MP00120   AT1G01250   PlantTFDB
# MP00100   AT1G01260   PlantTFDB

# Then use the code:
motif2TF <- motif2TF %>% dplyr::transmute(`#motif_id`=paste0(Motif_ID,'_', TF_Name), motif_name=paste0(Motif_ID,'_', TF_Name), motif_description=TF_Name,
                                          source_name='CIS_BP', source_version=1.1, gene_name=TF_Name,
                                          motif_similarity_qvalue=0.000000, similar_motif_id="None", 
                                          similar_motif_description="None", orthologous_identity=1.000000,
                                          orthologous_gene_name=Family_Name, orthologous_species="None", 
                                          description=DBDs)
write.table(motif2TF, file = "motif2TF.tbl", sep = "\t", row.names = F, quote = F)
````



## Run pyscenic

```shell
#!/bin/bash

source ~/.bashrc
mamba activate pyscenic

INPUT_LOOM="/lustre/home/jianghao/02_workspace/08_spatial_omics/21_dm4_count/bin50_data/s123_run.loom"

TFS_PATH="/lustre/home/jianghao/02_workspace/08_spatial_omics/18_GRN/pyscenic_databases/cis_bp/tf_list.txt"
FEATHER_PATHS="/lustre/home/jianghao/02_workspace/08_spatial_omics/18_GRN/pyscenic_databases/Stu_dm4_gene.regions_vs_motifs.rankings.feather"
TABLE_PATH="/lustre/home/jianghao/02_workspace/08_spatial_omics/18_GRN/pyscenic_databases/cis_bp/motif2TF.tbl"

# Run PySCenic GRN analysis with 20 threads
pyscenic grn \
    --num_workers 20 \
    --output grn.tsv \
    --method grnboost2 \
    "${INPUT_LOOM}" "${TFS_PATH}"

# Run PySCenic Cistarget analysis with 20 threads
pyscenic ctx \
    grn.tsv "${FEATHER_PATHS}" \
    --annotations_fname "${TABLE_PATH}" \
    --expression_mtx_fname "${INPUT_LOOM}" \
    --mode "dask_multiprocessing" \
    --output ctx.csv \
    --num_workers 20 \
    --mask_dropouts

# Run PySCenic AUCell analysis with 20 threads
pyscenic aucell \
    "${INPUT_LOOM}" \
    ctx.csv \
    --output aucell.loom \
    --num_workers 20
```

## Reference

[pyscenic 构建自己的 cisTarget 数据库](https://mp.weixin.qq.com/s/7-vKrLiFS4Tlkt-rHxEGeQ)

[pySCENIC分析全流程相关代码](https://mp.weixin.qq.com/s/xbmyAiLyWBJk0E1x6PSLWw)

Cao, S., He, Z., Chen, R., Luo, Y., Fu, L.-Y., Zhou, X., He, C., Yan, W., Zhang, C.-Y., & Chen, D. (2023). scPlant: A versatile framework for single-cell transcriptomic data analysis in plants. *Plant Communications*, *4*(5), 100631. https://doi.org/10.1016/j.xplc.2023.100631
