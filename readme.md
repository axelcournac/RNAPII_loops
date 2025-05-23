#  Cohesin and RNA Pol II loops 

Codes and scripts to reproduce the analyses and plots from the article _RNA Pol II-based regulations of chromosome folding_ by
Christophe Chapard, Nathalie Bastié, Axel Cournac, Laura Chaptal, Henri Mboumba, Sophie Queille, Agnès Thierry, Olivier Gadal, Armelle Lengronne, Romain Koszul, Frédéric Beckouët. 

### Dependencies

Scripts and codes can be run on OS X and other Unix-based systems. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python (>=3.6)
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Pandas
* scikit-misc 

#### External programs

* `SRA tools` / [SRA tools](https://github.com/ncbi/sra-tools)
* `Bowtie2` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* `hicstuff` / [hicstuff](https://github.com/koszullab/hicstuff)
* `cooler` / [cooler](https://github.com/open2c/cooler)
* `tinyMapper` / [tinyMapper](https://github.com/js2264/tinyMapper)
* `deepTools` / [deepTools](https://deeptools.readthedocs.io/en/develop/)
* `chromosight` / [chromosight](https://github.com/koszullab/chromosight)

### Raw data extraction and alignment
#### Data extraction
Data can be dowloaded on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.

A SRA executable called fastq-dump from SRA can be used to extract and split both mates of a library (to use it, you can go with your terminal to the directory containg the executables files by using the bash command cd).Then the program can be used like this:  /fastq-dump library_identification --split-3 -O /path_to_a_directory
 
```bash
./fastq-dump SRR639031 --split-3 -O /home/axel/data/
```

#### Alignment of the Hi-C libraries
To align the reads and generate the contact files in cool format, we used hicstuff pipeline: 
```bash
f1="LCH13_nxq_R1.fq.gz"
f2="LCH13_nxq_R2.fq.gz"

i="LCH13"
title2=$i

hicstuff pipeline -t 8 --mapping="iterative" --read-len=50 -f -D -a bowtie2 -e DpnII,HinfI --matfmt bg2 --no-cleanup -F -o out_$i -g /media/axel/RSG5/diverse_yeast_data_copy/cohesion_paper/genomes/W303_2micron $f1 $f2
```
and convert into cool file:
```bash
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 sacCer3.chr_sizes.txt:200 valid_idx_pcrfree.pairs valid_idx_pcrfree.pairs.cool
```

### Processing of genomic data like ChIP-seq or RNA-seq
We used tinyMapper: 
```bash
./tinyMapper.sh -m ChIP -s chIP-seq/LCH17_nxq -i chIP-seq/LCH16_nxq -g genomes/W303_2micron -c genomes/glabrata -o results_G2_arrested_cells_25C_Rpb1_ChIP

./tinyMapper.sh -m RNA -s SRR8503057.1 -g SC288_with_micron_SC88 -o results_RNAseq
```

### Visualisation of pileups and quantification 

####  Plot of the agglomerated plot with home made code
python3 agglomerated_pileup.py $contact_data $name_exp $pair_file $name_pair

####   Quantification with Chromosight 
chromosight quantify --pattern=loops  --perc-undetected=100 --perc-zero=100   $pair_file $contact_data $name_exp"_"$name_pair

####   Plot of the agglomerated plot from Chromosight (a little different from the home made plot) and plot of the distribution of pattern scores
python3 agglo_distrib.py $name_exp"_"$name_pair.json $name_exp"_"$name_pair.tsv


### Computation and visualisation of loop spectrum 

To compute the spectrum of cohesin loops at differents sizes, we use the python code **analysis_spectrum_cohesin.py** after quantifiying the loop scores with Chromosight with the following command: 
```bash
chromosight quantify --pattern=loops_small --perc-zero=100 --perc-undetected=20  /home/axel/Bureau/YEAST/pairs_peaks_cohesins3.txt.bg2.all_sizes.2  /media/axel/RSG5/disk/copy_diverse_yeast/data_LChaptel/out_LCH10/tmp/valid_idx_pcrfree.pairs.2000.S288C.cool  LCH10
```
To compute the spectrum of cohesin loops at differents sizes, we use the python code **analysis_spectrum_RNAPolII.py** after quantifiying the loop scores with Chromosight with the following command: 
```bash
chromosight quantify --pattern=loops --perc-undetected=100  --perc-zero=100  /home/axel/Bureau/YEAST/pairs_peaks_1.5_PolII_Log_1.5.txt.bg2 SRR13736654_55_tmp_valid_idx_pcrfree.pairs.2000.cool  SRR13736654_55
```


### Generation and visualisation of 4C like plot 
To generated 4C like plot at various genomic positions, we use the home made python code 4C_smoothed_pol2.py







