# tetraCall - convert diploid variant calling results to tetraploid information

Do you have the problem that tetraploid capable variant callers like *Freebayes* do perform less good than other variant callers? (like *samtools / bcftools*, *GATK*, *DeepVariant*,..) 

Here comes the solution. We comvert the diplod calles to tetraploid information by extracting the allelic depth informaation *AD* from your variant calling files and "recalculate" the tetraploid information - just as *Freebayes* does it. 

## Requirements:

1. you need to call the variants with **allelic depth** format *AD* 
    - bcftools does it like:

        ```
        # define the input variables
        REF="/RefGenome.fa" 
        mydir="StoreFolder" # e.g. ~/Documents

        bcftools mpileup -Ov -q 25 -Q 30 -a FORMAT/AD -I -f $REF $mydir/*filtered_sorted.bam | bcftools call -vm -Ov > $mydir/variants.vcf
        ```

2. you need to install `julia` and `vcftools`
    - `vcftools` needs to be added to the *path*, so that the script can access it
        - download from here: http://vcftools.sourceforge.net/
        - if not attached to the *path*, please add the path to *line 16* of `diploid_TO_tetra.jl`
    - `julia` depends on the packages `CSV and DataFrames`, so please install these accordingly
        - install packages as shown below:
        ![](https://github.com/mischn-dev/popRR/blob/docs/install_juliaPackages.gif)


After this is done, execute one bash line and the job will be executed for you. Easy? Easy!

## Run the transformation

Here is not much to explain - run it as any other bash script. Simply type:

```
# t = number of threads; more threads == faster progress
../../julia -t 10 ./diploid_TO_tetra.jl ./path_to/AD_VCFfile.vcf
```

NOTE - vcf files are not supported in *.gz* format. If you have a gunzipped file, change code of line 16 to: 
```
run(`vcftools --gzvcf $vcffile --out diploids --max-alleles 2 --min-alleles 2 --extract-FORMAT-info AD`)
```

## Data output

Some intermediate files will be created and deleted after the job is done. Do not mind these. 

All files will be generated in the location where you stored the original *VCF* file.

3 files are generated:

1. The actual alleles with the tetraploid information `Filename_Alleles.txt`
    - besides the Chromosome and the Position information, the *GT* information is presented for each sample included in the variant calling file

        | Chr | Pos  | Sample1  | Sample2 | Sample3 | .. |
        | ------- | --- | --- | --- | --- | --- |
        | chr2 |	268269	|0/0/1/1 | 	0/0/0/1	|1/1/1/1 | .. |
        | chr2 |	352795	| 0/1/1/1	| 0/0/0/0	| 0/1/1/1 | .. |
        .. | ..| ..| ..| ..|

2. The information of the read depth for this variant locus for each sample includes `Filename_ReadDepth.txt`
    - same structure as above

        | Chr | Pos  | Sample1  | Sample2 | Sample3 | .. |
        | ------- | --- | --- | --- | --- | --- |
        | chr2 |	268269	|12 | 	15	|5 | .. |
        | chr2 |	352795	| 22	| 18	| 23 | .. |
        .. | ..| ..| ..| ..|

3. The Reference allele frequency per variant locus / sample `Filename_RefAF.txt`
    - this file provides the "raw" allele frequency of the reference allele. So you can create the number of reads observed for reference / alternative allele in combintion with file *Filename_ReadDepth.txt*

        | Chr | Pos  | Sample1  | Sample2 | Sample3 | .. |
        | ------- | --- | --- | --- | --- | --- |
        | chr2 |	268269	| 0.48 | 	0.2	| 1.0 | .. |
        | chr2 |	352795	| 0.71	| 0.01	| 0.81 | .. |
        .. | ..| ..| ..| ..|

