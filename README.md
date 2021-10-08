# Pooled CRISPR Screen Data Processing

This guide will explain our workflow for processing *in vivo* pooled CRISPR screens...

## Generating Guide Fold Changes

### Installing Requirements 

##### with Anaconda:

1) Install anaconda, following documentation [here](https://docs.anaconda.com/anaconda/install/index.html)

2) Install requirements from yml file.   

    ```shell script
    conda env create -f envs/environment.yml
    ```
   
3) Activate conda environment 

    ```shell script
    conda activate screening_data_processing
    ```
    

##### with Docker:

1) Install docker following documentation [here](https://docs.docker.com/get-docker/)

2) Run jupyter notebook using [Jupyter Docker Stacks scipy-notebook](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html) docker container.
From the protocol directory, run:

    ```shell script
    docker run -v "${PWD}":/home/jovyan/ -p 8888:8888 jupyter/scipy-notebook 
    ```

3) Copy and paste the displayed URLs into any browser
 and open the jupyter notebook /notebooks/calc_guide_foldchanges.ipynb
 
### Running Workflow

notebooks/calc_guide_foldchanges.ipynb shows how to generate 
Quantile normalized fold changes for the *in vivo* screen performed in Dubrot et. al. 2021<sup>1</sup>

The general workflow is: 

1) Read counts are converted converted to RPMs with pseudocount of 1, and then log2 transformed

2) Gene-targeting guides are then z-normalized based on control sgRNA distribution. 

3) Guide fold changes are calculated as residuals fit to a natural cubic spline with 4 degrees of freedom. 

4) Fold changes from all four pools are quantile normalized with gene-wise mean
 imputation. Quantile normalized fold changes are used for downstream analysis and hit calling
 
 
## Downstream analysis of guide fold changes
 
 
### Installing Requirements 

##### with Anaconda:


##### with Docker:

### Hit calling with STARS:

The STARS algorithm<sup> 2</sup> can be used to find enriched an depleted genes from a screen.

to run: 

1) First generate STARS null distribution using /scripts/stars_null.py. 

    ```shell script
    python scripts/stars_null.py \
    --input-file /data/Renca_zresid_guides.txt \
    --chip-file /data/MS_tide_data/uniq_ID_chip.txt \
    --thr 20 \
    --max 4 \
    --num-ite 250 \
    --output-dir data/outputs/
    ```
2) Run stars hit calling using the null file generated in step 1 

    To find depleted genes 
     
    ```shell script
    python scripts/stars.py \
    --input-file ./data/Renca_zresid_guides.txt \
    --chip-file ./data/uniq_ID_chip.txt \
    --dir N \
    --null data/outputs/Null_STARSOutput8_thr20_max4.txt \
    --max 4 --thr 20 \
    --output-dir data/outputs/STARS_N
    ```
    
    to find enriched genes 
    
    ```shell script
    python scripts/stars.py \
    --input-file ./data/Renca_zresid_guides.txt \
    --chip-file ./data/uniq_ID_chip.txt \
    --dir P \
    --null data/outputs/Null_STARSOutput8_thr20_max4.txt \
    --max 4 --thr 20 \
    --output-dir data/outputs/STARS_P
    ```
     

### Running hypergeom analysis:

```shell script
python scripts/hypergeom.py \
--input-file data/Renca_zresid_guides.txt \
--output-dir data/outputs/hypergeom \
--chip-file data/uniq_ID_chip.txt \
--max 4
```


citations: 
 1) Dubrot J, Lane-Reticker SK, Kessler EA, Ayer A, Mishra G, Wolfe CH, et al. In vivo screens using a selective CRISPR antigen removal lentiviral vector system reveal immune dependencies in renal cell carcinoma. Immunity [Internet]. 2021 Jan 20; Available from: http://dx.doi.org/10.1016/j.immuni.2021.01.001
 2) Doench JG, Hartenian E, Graham DB, Tothova Z, Hegde M, Smith I, et al. Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation. Nat Biotechnol. 2014 Dec;32(12):1262–7.
 