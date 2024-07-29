# _A systems biology approach to elucidate the association between osteoporosis and sarcopenia_



## Overall Workflow

![workflow figure](./Figure/Figure1_overall_workflow.pdf "Overall workflow")





## Tutorial

### Module 1. Differential gene expression analysis

Please run the `` script in R or Rstudio.

### Module 2. Robust Rank Aggregation (RRA)

Please run the `` script in R or Rstudio.

### Module 3. Network analysis

Please use [STRING](https://string-db.org) and [Cytoscape](https://cytoscape.org) for network analysis.

### Module 4. Functional enrichment analysis

### Module 5. Machine learning framework



- Option 1. Single sample 

```bash
sample_name='example_sample'
bam_cram_file=/path/to/example_sample.cram
windows_file=${project_dir}/windows.bed
mkdir -p ${output_rd_dir}

mosdepth -n --fasta ${reference_file} --by ${windows_file} --mapq 30 ${output_rd_dir}/${sample_name} ${bam_cram_file}
zcat ${output_rd_dir}/${sample_name}.regions.bed.gz | cut -f1-3,5 | bgzip -c > ${output_rd_dir}/${sample_name}.cov.bed.gz
tabix -f -p bed ${output_rd_dir}/${sample_name}.cov.bed.gz
```

This is a relatively time-consuming step, however, you can use the following commands to accelerate the speed if you have a cluster. Otherwise, run the above commands in a loop.

- Option 2. Multiple samples (by cluster)

```bash
bam_cram_file_path_list=${project_dir}/sample_rd_file_list.txt
num_tasks=`wc -l ${bam_cram_file_path_list} | cut -f1 -d" "`
windows_file=${project_dir}/windows.bed

# Suppose we want to calculate RD for all samples using cluster
## By SGE cluster
qsub -t 1-${num_tasks} ${script_dir}cluster_mosdepth.sh \
    ${bam_cram_file_path_list} ${reference_file} ${windows_file} ${output_rd_dir}      

## By slurm workload manager
sbatch -a 1-${num_tasks} ${script_dir}cluster_mosdepth.sh \
    ${bam_cram_file_path_list} ${reference_file} ${windows_file} ${output_rd_dir}
```

Collect all the coverage files for downstream usage. 

```bash
find ${output_rd_dir} -maxdepth 1 -name "*.cov.bed.gz" > ${project_dir}/sample_raw_rd.txt
```

### Step 3. GC normalization

- Option 1. Single sample

```bash
python ${script_dir}cnv_espresso.py normalization \
    --windows ${project_dir}/windows.bed \
    --input   ${output_rd_dir}/example_sample.cov.bed.gz \
    --output  ${project_dir}/norm/
```

- Option 2. Multiple samples (input a simple list) 

```bash
python ${script_dir}cnv_espresso.py normalization \
    --windows    ${project_dir}/windows.bed \
    --input_list ${project_dir}/sample_raw_rd.txt \
    --output     ${project_dir}/norm/
```

- Option 3. Multiple samples (by cluster) 

```bash
windows_file=${project_dir}/windows.bed
RD_list=${project_dir}/sample_raw_rd.txt
output_dir=${project_dir}/norm/
num_tasks=`wc -l ${RD_list} | cut -f1 -d" "`

## By SGE cluster
qsub -t 1-${num_tasks} ${script_dir}cluster_gc_norm.sh \
    ${script_dir} ${windows_file} ${RD_list} ${output_dir}

## By Slurm workload manager
sbatch -a 1-${num_tasks} ${script_dir}cluster_gc_norm.sh \
    ${script_dir} ${windows_file} ${RD_list} ${output_dir}
```

### Step 4. Select reference samples

```bash
ls ${project_dir}/norm/*.gz > ${project_dir}/sample_norm_rd.txt
```

```bash
python ${script_dir}cnv_espresso.py reference \
    --project_dir ${project_dir} \
    --norm_list   ${project_dir}/sample_norm_rd.txt \
    --num_ref     100 \
    --corr_threshold -1 
```

### Step 5. Encode images

Here, we will take the CNV caller's output (e.g. `xhmm.xcnv`) as our input and we will use the following commands to encode those CNV predictions into images.

```bash
RD_norm_dir=${project_dir}/norm/
ref_samples_dir=${project_dir}/ref_samples/
output_dir=${project_dir}
cnv_list=${project_dir}/xhmm.xcnv # or other CNV caller's output
```

- Option 1. Generate images via single thread

```bash
python ${script_dir}cnv_espresso.py images \
    --rd_norm_dir ${RD_norm_dir} \
    --ref_dir     ${ref_samples_dir} \
    --cnv_list    ${cnv_list} \
    --output      ${output_dir} \
    --overwrite_img False
```

- Option 2. Generate images by a cluster
  Note: please modify the path of script in the `cluster_images.sh` file at first.

```bash
num_tasks=`wc -l ${cnv_list} | cut -f1 -d" "`
overwrite_img=False

# By SGE cluster
qsub -t 1-${num_tasks} ${script_dir}cluster_images.sh \
    ${script_dir} ${RD_norm_dir} ${ref_samples_dir} ${cnv_list} ${output_dir} ${overwrite_img} 

# By Slurm workload manager
sbatch -a 1-${num_tasks} ${script_dir}cluster_images.sh \
    ${script_dir} ${RD_norm_dir} ${ref_samples_dir} ${cnv_list} ${output_dir} ${overwrite_img}

# Collect info of CNVs and images processed by cluster into a single file.
wc -l ${output_dir}cnv_info_w_img.csv

find ${output_dir}/single_img_info/ -maxdepth 1 -type f -name '*.csv' -print0 |
sort -zV |
xargs -0 cat >${output_dir}cnv_info_w_img_collected.csv
```

A few example images are located [here](https://github.com/ShenLab/CNV-Espresso/tree/main/example/images). 

### Step 6. Training (Optional)

In general, you can directly use our pretrained CNN model ([MobileNet_v1_fine-tuning_3classes_log_transformed.h5](https://github.com/ShenLab/CNV-Espresso/blob/main/model/MobileNet_v1_fine-tuning_3classes_log_transformed.h5)) for your *in silico* validation (Skip step 6). However, if you have a bunch of validated or high-confidence CNVs, you can also train the CNN model from scratch. If so, please follow the tutorial below:

1. Please use `images` function in `cnv_espresso.py` as **Step5** to generate images for your prepared true deletion, true duplication, false deletion and false duplication. Note that false deletion and false duplication will be treated as diploid together in the downstream steps.

2. We recommend using a GPU server to handle this step. You can train your custom model by the following commands:

```bash
true_del_img=${project_dir}/train/true_del.list
true_dup_img=${project_dir}/train/true_dup.list
false_del_img=${project_dir}/train/false_del.list
false_dup_img=${project_dir}/train/false_dup.list
output_dir=${project_dir}/model/

python ${script_dir}cnv_espresso.py train \
    --true_del  ${true_del_img} \
    --true_dup  ${true_dup_img} \
    --false_del ${false_del_img} \
    --false_dup ${false_dup_img} \
    --use_gpu   True \
    --output    ${output_dir}
```

Alternatively, we also prepared a jupyter notebook (**[train.ipynb](https://github.com/ShenLab/CNV-Espresso/blob/main/src/train.ipynb)**) for tracking and debugging the entire training process.

### Step 7. Confirming CNV predictions *in silico*

```bash
cnv_w_img_file=${project_dir}/cnv_info_w_img.csv (or cnv_info_w_img_collected.csv)
model_file='/path/to/model/MobileNet_v1_fine-tuning_3classes_log_transformed.h5'
output_file=${project_dir}/cnv_espresso_prediction.csv

python ${script_dir}cnv_espresso.py predict \
    --cnv_list ${cnv_w_img_file} \
    --model    ${model_file} \
    --output   ${output_file} \
    --use_gpu  False
```

## Utilities

1. Collect the absolute path of images and annotate this information to ‘cnv_info_w_img.csv’

```bash
python ${script_dir}cnv_espresso.py collect_images \
    --cnv_file ${cnv_w_img_file} \
    --output_dir ${output_dir}
```

## Citation

[TBD] Submitted to a peer-reviewed journal.

## Contact

Renjie Tan (renjie.tan at hrbmu.edu.cn) 
