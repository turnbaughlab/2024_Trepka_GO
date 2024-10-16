This repository contains all processed data and code for generating figures in the manuscript "Expansion of a bacterial operon during cancer treatment ameliorates drug toxicity."
 - The SequencingDataREADME.rtf file contains instructions for adding metadata to files downloaded from the NCBI Sequence Read Archive (PRJNA1169175)
 - The GO_Figures_241016.Rmd file contains modular code for generating all figures. When running the .Rmd file, a new folder will be generated for each main or supplemental figure, with relevant PDFs saved into this folder
 - The AmpliconFunctions.R file contains several helper functions for 16S rRNA amplicon sequencing analysis
 - All other folders represent processed and raw data
 - To recapitulate original file structure
     - The .Rmd file should be in a layer with a superfolder entitled "CompiledData"
     - All zipped data folders should be unzipped and moved into the "CompiledData" folder
     - The AmpliconFunctions.R file should be placed inside the "CompiledData" folder
     - The UniRef90...csv.zip file should be unzipped and moved into the "CompiledData/GO_MGS" folder
  
  
