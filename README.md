This repository contains all processed data and code for generating figures in the manuscript "Expansion of a bacterial operon during cancer treatment ameliorates drug toxicity."
 - The SequencingDataREADME.rtf file contains instructions for adding metadata to files downloaded from the NCBI Sequence Read Archive (PRJNA1169175)
 - The .Rmd file contains modular code for generating all figures. When running the .Rmd file, a new folder will be generated for each main or supplemental figure, with relevant PDFs saved into this folder
 - The .R file contains several helper functions for 16S rRNA amplicon sequencing analysis
 - All other folders represent processed and raw data
 - To recapitulate original file structure
     - The .Rmd file should be in a layer with a folder entitled "CompiledData"
     - All other data folders and the .R amplicon seq file should go in the "CompiledData" folder
     - The UniRef90... .csv.zip file should be moved into the GO_MGS folder
  
  
