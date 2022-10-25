if (!require(tools)) install.packages("tools")
if (!require(data.table)) install.packages("data.table")
if (!require(ggvenn)) install.packages("ggvenn")
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(UpSetR)) install.packages("UpSetR")

library(tools)
library(data.table)
library(ggvenn)
library(dplyr)
library(ggplot2)
library(UpSetR)

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("\nAt least 2 samples must be compared. Use the argument 'help' for more information.", call.=FALSE)
} else if (args[1] == "help") {
  stop("\nThis script takes a minimum of 2 MAF files as input to compare mutations.
       Inputs should be given with full paths. Example:
       'script.R <full/filepath/to/sample_1.maf> <full/filepath/to/sample_2.maf> <full/filepath/to/sample_3.maf>'\n", call.=FALSE)
} else if (length(args) == 1 & args[1] != "help") {
  stop("\nAt least 2 samples must be compared. Use the argument 'help' for more information.", call.=FALSE)
} else if (length(args) >= 2) {
  # Check to make sure files are valid (either MAF or TXT); quit if they aren't
  for (i in seq_along(args)) {
    if (file_ext(args[i]) != "maf" && file_ext(args[i]) != "txt") {
      stop("\nInvalid file entered; all sample files must be either .maf or .txt format.", call.=FALSE)
    }
  }
  
  # Set the working directory
  cat("Please enter your working directory's full filepath: ")
  read_filepath <- readLines("stdin", n = 1)
  if (file.exists(read_filepath)) {
    setwd(read_filepath)
  } else {
    dir.create(read_filepath)
    setwd(read_filepath)
  }

  # All files and working directory validated, so we can continue
  number_of_files <- length(args)
  
  sample_names <- list()

  for (i in seq_along(args)) {
    sample_names[i] <- basename(gsub("\\..*", "", args[i]))
  }
}


sample_data <- list()

print("Sample stats file will be outputted with the naming format 'sample_stats-<month>-<day>-<hour>-<minute>'")
sys_time <- format(Sys.time(), "%m-%d-%H-%M")
stats_filename <- paste0("sample_stats-", sys_time)
stats_filename <- paste0(stats_filename, ".txt")

### Stats on the samples ###

for (i in seq_along(args)) {
  sample_data[[i]] <- fread(file = args[i], header = TRUE, sep = "\t", data.table = FALSE, select = c("Tumor_Sample_Barcode", "Hugo_Symbol", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type"))
  current_sample <- sample_names[i]
  num_mutations <- nrow(sample_data[[i]])
  num_unknown <- sum(sample_data[[i]]$Hugo_Symbol == "Unknown")

  sink(file = stats_filename, append = TRUE)
  print(sprintf("Sample: %s", current_sample))
  print(sprintf("Total mutations: %d", num_mutations))
  print(sprintf("Total unknown genes: %d", num_unknown))
  sink()

  #print(sprintf("There are %d total mutations, %d of which are unknown genes, in sample %s", y, x, sample_names[i]))
}

check_cond <- 0
while (check_cond == 0) {
  cat("Enter 1 to compare gene mutation overlap between samples, enter 2 to compare variant types between samples, or enter 3 to do both: ")
  read_var <- readLines("stdin", n = 1)
  if (read_var != 1 && read_var != 2 && read_var != 3) {
    print("That is not a valid option. Please try again.")
    next
  } else {
    check_cond <- 1
  }
}

combined_sample_data <- bind_rows(sample_data)




### Venn diagram of gene mutation overlap between samples ###
if (read_var == 1 || read_var == 3) {
  combined_sample_data_venn <- data.frame(combined_sample_data)
  combined_sample_data_venn <- combined_sample_data_venn[(combined_sample_data_venn$Hugo_Symbol != "Unknown"),]
  combined_sample_data_venn$GeneAllele <- paste(combined_sample_data_venn$Hugo_Symbol, combined_sample_data_venn$Tumor_Seq_Allele2, sep = "_")

  venn_list <- list()
  for (i in seq_along(sample_names)) {
    current_sample <- sample_names[i]
    current_set <- assign(paste0("set_", i), subset(combined_sample_data_venn, Tumor_Sample_Barcode == sample_names[i]))

    venn_list[i] <- list(current_set$GeneAllele)
    names(venn_list)[i] <- current_sample
  }
  
  if (number_of_files <=4) {
    output_file <- "venn_plot.pdf"
    ggvenn(data = venn_list,
           show_percentage = FALSE,
           stroke_linetype = 2,
           stroke_size = 0.25,
           set_name_size = 4.8,
           text_size = 6
    )
    ggsave(output_file)
    print(sprintf("***** Finished creating Venn diagram file. Location: %s/%s. *****", getwd(), output_file))
  } else {
    pdf(file = "upset_plot.pdf")
    upsetPlot <- upset(fromList(venn_list), order.by = "freq")
    print(upsetPlot, newpage = FALSE)
    dev.off()
    print(sprintf("***** Finished creating upset plot file. Location: %s/upset_plot.pdf. *****", getwd()))
  }
}


### Bar graph of mutation type comparison ###
if (read_var == 2 || read_var == 3) {
  combined_sample_data_bar <- subset(combined_sample_data, select = c(Tumor_Sample_Barcode, Variant_Classification, Variant_Type))
  
  combined_sample_data_bar <- combined_sample_data_bar[(combined_sample_data_bar$Variant_Classification == "IGR" | 
                                                          combined_sample_data_bar$Variant_Classification == "Intron" |
                                                          combined_sample_data_bar$Variant_Classification == "RNA"),]
  combined_sample_data_bar <- combined_sample_data_bar[(combined_sample_data_bar$Variant_Type == "SNP" | 
                                                          combined_sample_data_bar$Variant_Type == "INS" |
                                                          combined_sample_data_bar$Variant_Type == "DEL"),]

  output_file <- "bar_plot.pdf"
  ggplot(combined_sample_data_bar, aes(x = Variant_Type, fill = Variant_Classification)) +
    geom_bar(position = position_dodge(width = 0.9)) +
    geom_text(aes(label = ..count..), stat = "count", position = position_dodge(width = 0.9), colour = "black", size = 3) +
    theme_minimal(base_size = 12) +
    ylab("Number of Variants") +
    xlab("Variant Type") +
    scale_y_continuous(trans = 'log10') +
    facet_grid(.~Tumor_Sample_Barcode, scales = "free_x", space = "free_x",switch = "x")
  ggsave(output_file, width = 14, height = 14)
  print(sprintf("***** Finished creating bar diagram file. Location: %s/%s *****", getwd(), output_file))
}
