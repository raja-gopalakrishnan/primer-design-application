# Primer design application
Primer designing tool developed in R shiny for deleting and tagging *S. cerevisiae* genes using the pFA6a cassettes developed by [Longtine et al](https://onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291097-0061%28199807%2914%3A10%3C953%3A%3AAID-YEA293%3E3.0.CO%3B2-U).

## Running the script

### 1. Install required packages
```
#Install the packages - Shiny, XML, stringr, stringi
install.packages("shiny")
install.packages("XML")
install.packages("stringr")
install.packages("stringi")

#Install Biostrings packages from Bioconductor
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
```
### 2. Change the ```path``` variable in the script
Open the app.R script. Change the ```path``` variable to the directory which has the app.R script.
```
path="path/to/primer-design-application/"
```
Remember to save the file after you make the changes!

### 3. Run the application
```
runApp("path/to/primer-design-application/")
```

## Using the application
#### 1. Enter the gene name
The gene name can be the common name or the systematic name.
#### 2. Select the experiment
Select if you want primers for deleting or tagging the gene
#### 3. Enter number of homologous nucleotides
This number specifies the number of overhanging homologous nucleotides for directing the deletion/tagging cassette to the right location in the genome.
#### 4. Select the tag (if tagging gene) and selection cassette (optional)
These values are used to return the plasmid number from our lab database that can be used for the deletion/selection.
#### 5. Click Submit!

## Interpreting the results

- The FB number indicates the plasmid number to be used from the Winston lab collection
- The 5' and 3' primers for tagging or deletion are displayed next
- The checking primer sequence contains a primer that anneals just outside the site where the deletion or tagging cassette is inserted. This primer is selected just base on GC content and melting temperature. The primer sequence has not been run through BLAST to check if it binds elsewhere in the yeast genome.
- The gene sequence displays the sequence of the gene chosen for deletion or tagging along with the 1 kb of 5' and 3' flanking sequences. Text in red indicates the checking primer sequence. Text in blue indicates the 40 bp 5' and 3' homology sequence. Text in yellow indicates the portion of the gene that will get replaced by the insertion cassette.

## Troubleshooting
If you get the following error, it is due to an invalid gene name that has been entered.
```
Error: arguments imply differing number of rows: 1, 0
```
Check that you have entered the gene name correctly. Try entering the systematic name for the gene instead of the common name. You can find the systematic name for your favorite gene on [SGD](https://www.yeastgenome.org/). The application will identify any of the 6819 genes present in the sgd_ids.tsv file.
