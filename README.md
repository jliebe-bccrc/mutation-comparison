# mutation-comparison
The MutCompare script (created by Jenna Liebe, Oct. 2022, BCCRC - Aparicio Lab) takes 2 or more cBioPortal mutation data files (MAF or TXT ) as input, and provides options for analyzing the data within them:
1. Create a Venn diagram showing the overlap of gene mutations in the samples, based on the name of the mutated gene (Hugo_Symbol) and the altered allele (Tumor_Seq_Allele2),
2. Create a bar chart showing the types of mutations present in each sample (SNPs, insertions, and deletions, grouped by the classification of IGR, RNA, or intron),
3. Create both the Venn diagram and the bar chart.

A sample_stats text file will also be created, which contains the sample IDs of each inputted file, along with their total number of mutations and their number of unknown genes. It is outputted using the month/date and time from the system, to avoid appending to existing files, and has been formatted to avoid using illegal characters in the filename. For example, running the script at 3:02:32 PM on October 24th will produce the sample_stats file sample_stats-10-24-15-02-32.txt.

Beyond 4 samples total, a Venn diagram becomes uninterpretable. Therefore, if the user inputs more than 4 samples, an upset plot will automatically be generated instead.

---

### Running the Script
The script is run via the command line, from wherever you generally run R (Rscript.exe, R.exe, etc.):

<Rscript.exe or R.exe> "path/to/mutCompare.R" "path/to/sample1.min.maf" "path/to/sample2.min.maf"

The input files will be validated to ensure they are of the proper file format; the script will throw an error if any of them are incorrect.
The user will then be prompted to enter the full filepath to their working directory - this directory will be where the graphs are written to.
In the script, there is a check to make sure that the specified directory exists. If not, a new directory is created.
If you get an error relating to creating the output file(s), make sure you have the correct permissions for writing to your specified working directory.

If your input files were validated earlier but are missing any of the five required columns (Tumor_Sample_Barcode, Hugo_Symbol, Tumor_Seq_Allele2, Variant_Classification, or Variant_Type) the script will fail. Please ensure that these columns are present in each of your input files.

The Venn diagram/upset plot and/or bar chart will be created and saved under the names “venn_plot.pdf”/"upset_plot.pdf" or “bar_plot.pdf”, respectively, in the working directory you specified earlier. Please make sure you rename these files if you want to rerun the script (or move them to a different directory), as any file with those names will be overwritten the next time the script is run.
