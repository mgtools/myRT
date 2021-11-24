# MyRT: identification and classification of reverse transcriptases in bacterial genomes and metagenomes

## To make sure myRT pipeline works properly please follow these steps: 

All the files needed for executing the pipelines are included in this github repository except cdd-pfamA.hmm files: 

### Download cdd-pfamA.hmm files

Use one of these commands to download the hmm files: 

`wget https://omics.sice.indiana.edu/myRT/Models/cdd-pfamA.tar.gz -P Models/`

`curl -L https://omics.sice.indiana.edu/myRT/Models/cdd-pfamA.tar.gz --output Models/cdd-pfamA.tar.gz`

`wget https://omics.informatics.indiana.edu/myRT/share/cdd-pfamA.tar.gz -P Models/`

To decompress: 

`tar xvzf Models/cdd-pfamA.tar.gz -C Models/`

Add read/execute permission if needed

`chmod -R 777 Models/cdd-pfamA/`  OR  `chmod a+rx Models/cdd-pfamA/ -R`

After Successful uncompress delete the Zip file: `rm -f Models/cdd-pfamA.tar.gz`

### Check the file permissions in myRT folder:

Use `pwd` to make sure you're in the main folder (myRT).

Make sure that you have permission to read and execute files in this folder.

You can use `chmod -R 777 *` or `chmod a+rx * -R` to add permission for all.

### Install FragGeneScan 

To install FGS Run these commands

`current=$(pwd)`

`cd bin/FragGeneScan*`

`make clean`

`make fgs`

`cd ${current}`


### Test the Executable Files (if needed)

You can try these commands one by one to make sure all the executable files work

`Rscript`

`bin/pplacer-Linux-v1.1.alpha19/pplacer`

`bin/FragGeneScan1.31/run_FragGeneScan.pl`

`bin/blast+/bin/blastp`

`bin/hmmer3.2/src/hmmscan`

If **hmmscan** doesn't work, download the latest version of hmmer into bin folder and follow their installation instructions.

If you install/use a new version of these into bin folder, make sure you update the folder name in `myRT-local-py3.py` too.  

### Test Jplace Parser Rscript

Try this command to see if Jplace Parser (Scripts/JplaceParser.R) works 
```
Rscript Scripts/JplaceParser.R QUERY/CP016020.1/CP016020.1-combo.jplace
```

If you already have treeio and castor package, you don't need to install them.

Please install a package only if you get an error that the package wasn't found. 

If you need to install **treeio** then use the commands below, and modify `Scripts/JplaceParser.R` to include the path to the library: 

```
mkdir ~/r-packages
R CMD INSTALL --library=~/r-packages bin/treeio_1.14.0.tar.gz
R CMD INSTALL --library=~/r-packages bin/castor_1.6.4.tar.gz
```

Next, modify the Scripts/JplaceParser.R use these packagee, for example:

```
library("treeio",lib.loc="~/r-packages")
library("castor",lib.loc="~/r-packages")
```

Then, try running JplaceParser.R again: 

`Rscript Scripts/JplaceParser.R QUERY/CP016020.1/CP016020.1-combo.jplace`

You can check **myRT-R.log** to see what the output of Jplaceparser looks like. 

You should also see this line if the jplace file is successfully parsed: 

CP016020.1_3624763_3626067_+ RVT-UG3 1.00 {1887} leaf

This shows 
RT gene, prediction, confidence, {edgenum}, leaf/non-leaf

The result of parsing the jplace file will be saved in QUERY/CP016020.1/CP016020.1-pred.txt. 

These files only exist if phylogenetic information were needed for a more accurate prediction, otherwise, these steps will be skipped.

You can modify the script to enforce the phylogenetic analysis, but the default is to skip if not needed. 

### Run myRT

Execute this command
```
./rerun.sh CP008696.1 QUERY
```

Make sure there's no "/" after QUERY.

Check the log file in myRT folder to see what to expect from these two commands.

### Check the output files
If the run is successful you should have these two files, along with other files in `QUERY/CP008696.1` folder:

```
ls QUERY/CP008696.1/CP008696.1-gn.gff

ls QUERY/CP008696.1/CP008696.1-RT.gff
```

You can compare it with `QUERY/CP016020.1` folder if you want to make sure you have all the output files.

`QUERY/BA000019.2` has less files in this folder compared to `CP016020.1`, because for this query the pplacer step was skipped (See Figure1 of our manuscript for more info).

See **QUERY/README** for a description of output files. 

### Download a reference genome from ncbi (if needed)

When all the dependencies are installed, try this command in the terminal to download this reference genome from ncbi: 

```
./Scripts/download_ncbi.sh CP008696.1
```

If successful you should have a new folder (CP008696.1) in the QUERY folder.

To download to a new folder (use any name you'd like): 

```
mkdir Newfolder
./Scripts/download_ncbi.sh CP008696.1  Newfolder
```

If successful you should have a new folder (CP008696.1) in the Newfolder folder.

If you'd like to download a list of reference genomes prepare a file like list_genomes and run the following command:

list_genomes should contain one accession number per line. 

For example, CP008696.1 in line 1 and CP016020.1 in line 2 with no empty lines,and no space after accession number.

```
cat list_genomes | while read line ; do ./Scripts/download_ncbi.sh ${line} ; done
```

##  Other materials

### The training data:
https://omics.informatics.indiana.edu/myRT/share/download.php

### List of RT classes and compiled collection of RTs: 

myRT-collection:
https://omics.sice.indiana.edu/myRT/collection.php
