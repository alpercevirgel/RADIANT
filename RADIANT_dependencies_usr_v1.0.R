# Dependencies for RADIANT_usr_v1.0

# I. Pre-processing ----

  ## Chunk 1: setting wd, installing & loading required packages and checking fcs files present. ----
  
  reqPackages = c("readxl", "flowCore", "flowViz", "ggplot2", "tidyverse", "flowStats", "DelayedArray", "xlsx", "Rcpp", "limma", 
                  "RColorBrewer", "CATALYST", "SummarizedExperiment", "BiocParallel", "matrixStats", "dplyr", "Biobase", "GenomicRanges",
                  "GenomeInfoDb", "IRanges", "S4Vectors", "ExperimentHub", "AnnotationHub", "BiocFileCache", "dbplyr", 
                  "BiocGenerics", "diffcyt", "cowplot", "knitr", "BiocStyle", "flowAI", "pdftools", "metR", "gridExtra", "reshape2",
                  "PeacoQC", "ggpubr", "pheatmap", "tidyr", "FlowSOM", "openCyto", "BiocManager", "Spectre", "CytoNorm", "premessa", "rstatix", "digest", "umap")
  
  cranPackages = c("readxl", "ggplot2", "tidyverse", "xlsx", "Rcpp", "RColorBrewer", "matrixStats", "dplyr", "dbplyr", "cowplot", "knitr", 
                   "pdftools", "metR", "gridExtra", "reshape2", "ggpubr", "pheatmap", "tidyr", "BiocManager", "rstatix")
  
  biocPackages = c("flowCore", "flowViz", "flowStats", "DelayedArray", "limma", "CATALYST", "SummarizedExperiment", "BiocParallel", 
                   "Biobase", "GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors", "ExperimentHub", "AnnotationHub", "BiocFileCache",
                   "BiocGenerics", "diffcyt", "BiocStyle", "flowAI", "PeacoQC", "FlowSOM", "openCyto")
  
  # updating gitPackages requires adding lines in the autoInstallpkg function as well
  gitPackages = c("Spectre", "CytoNorm", "premessa")

  autoInstallpkg = function(reqPackages. = reqPackages, cranPackages. = cranPackages, biocPackages. = biocPackages, gitPackages. = gitPackages){
    
    message("Indexing missing packages...")
    new.packages = reqPackages.[!(reqPackages. %in% installed.packages()[,"Package"])]
    
    if(length(new.packages)>0){
      
      new.biocPackages = c(biocPackages.[which(biocPackages. %in% new.packages)])
      new.cranPackages = c(cranPackages.[which(cranPackages. %in% new.packages)])
      new.gitPackages = c(gitPackages.[which(gitPackages. %in% new.packages)])
      
      if (length(new.cranPackages)>0){cat("Missing CRAN packages: ", new.cranPackages, "\n")}
      if (length(new.biocPackages)>0){cat("Missing Bioconductor packages: ", new.biocPackages, "\n")}
      if (length(new.gitPackages)>0){cat("Missing gitHub packages: ", new.gitPackages, "\n")}    
      
      userInput = readline(prompt = "Do you want to install missing packages automatically? (y/n) ")
      
      if(userInput == "y"){
        
        #install CRAN packages
        if (length(new.cranPackages)>0){
          cat("Installing CRAN packages...\n")
          install.packages(new.cranPackages, dependencies = TRUE)
          cat("\n> CRAN packages successfully installed.\n\n")
        }
        
        #install Bioconductor packages
        if(length(new.biocPackages)>0){
          print("Installing Bioconductor packages...\n")
          library(BiocManager)
          for(bpkg in new.biocPackages){
            tryBioc = tryCatch(BiocManager::install(bpkg), error = function(e)e)
            if(!inherits(tryBioc, "error")){
              BiocManager::install(bpkg, ask = FALSE)
              
            }else{
              cat("\n")
              stop("Failed to install ", bpkg, ' from Bioconductor. Please install manually by typing BiocManager::install("',bpkg,'", ask = FALSE) in the console')
              cat("\n")
            }
          }
          
          cat("> Bioconductor packages successfully installed.\n")
        }
        
        #install gitHub packages
        if(length(new.gitPackages)>0){
          cat("Installing github packages...")
          
          if(!require('devtools')) {install.packages('devtools')}
          library('devtools')
          install_github("immunedynamics/spectre")
          install_github('saeyslab/CytoNorm')
          install_github("ParkerICI/premessa")
          
          cat("> Github packages successfully installed.\n")
        }
        
      }else{
        stop("Please install missing packages.")
      } 
      
    }else {
      cat("> All required packages are already installed.\n\n")
    }
  }
  
  autoLoadpkg = function(.reqPackages = reqPackages){
    message("Loading required packages...")
    lapply(.reqPackages, library, character.only = TRUE)
    cat("> Packages successfully loaded\n\n")
  }
  
  mapDirectories <- function() {
    
    # Map locations of the different files and documents
    assign("fcs", file.path(wd, cellType, "01 FCS files", "00 raw files"), envir = .GlobalEnv)
    assign("trns", file.path(wd, cellType, "01 FCS files", "01 transformed files"), envir = .GlobalEnv)
    assign("rfs", file.path(wd, cellType, "02 Reference files", "01 FCS files"), envir = .GlobalEnv)
    assign("cms", file.path(wd, cellType, "02 Reference files", "00 Compensation matrices"), envir = .GlobalEnv)
    
    if (RMD %in% c("Pre_processing", "FlowSOM_clustering", "Downstream_analysis")) {
      rsl <- file.path(wd, cellType, "03 Results", sprintf("%02d %s", match(RMD, c("Pre_processing", "FlowSOM_clustering", "Downstream_analysis")), RMD))
      assign("rsl", rsl, envir = .GlobalEnv)
      
    }
    
    assign("pnl",  file.path(wd, cellType, sprintf("%s Catalyst %s panel.xlsx", expID, cellType)), envir = .GlobalEnv)
    assign("mtd", file.path(wd, cellType, sprintf("%s Catalyst %s metadata.xlsx", expID, cellType)), envir = .GlobalEnv)
    
    # Create folder structure
    dirs <- c(fcs, rfs, cms, rsl)
    for (dir in dirs) {
      if (!file.exists(dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      }
    }
  }
  
  checkPreviousRuns <- function() {
    fcs_dir <- file.path(wd, cellType, "01 FCS files")
    
    # Check and print files in "00 raw files" directory
    raw_files_dir <- file.path(fcs_dir, "00 raw files")
    raw_files <- dir(raw_files_dir)
    
    if (length(raw_files) == 0) {
      cat("This is the first time you've run this chunk. Upload the FCS files to the 'cellType/00 FCS files/00 raw files' folder and set the 'to_import' variable in the next chunk to a.\n")
    } else {
      
      message("The following fcs files were detected:")
      
      check_and_print_files(raw_files_dir, "raw files")
      check_and_print_files(file.path(fcs_dir, "01 transformed files"), "transformed files")
      check_and_print_files(file.path(fcs_dir, "02 cleaned and transformed files/Preprocessed/QC/PeacoQC_results/fcs_files"), "cleaned and transformed files")
      check_and_print_files(file.path(fcs_dir, "03 normalized files"), "normalized files")
      
      if(RMD != "Downstream_analysis"){
      message("\nIn the next chunk, the only options for the 'to_import' variable are the options mentioned above!\n")
      }
      if(RMD == "Downstream_analysis"){
        message("\nIn the next chunk, variables of the clustering document are loaded and printed on screen!\n")
      }
    }
  }
  
  check_and_print_files <- function(path, label) {
    if (dir.exists(path)) {
      file_list <- list.files(path)
      file_count <- length(file_list)
      if (file_count > 0) {
        cat("> ", label, " (", file_count, " files)\n", sep = "")
      }
    }
  }
  
  ## Chunk 2: loading data into the working environment. ----
  
  import_fcs_files = function(){
    
    # imports raw data  
    if(to_import == "a") {
      
      path_fcs = fcs
      fs = read.flowSet(path=path_fcs,transformation=FALSE, truncate_max_range=FALSE)
      cleanFile = FALSE
      assign("cleanFile", cleanFile, envir = .GlobalEnv)
      
      cat("Raw files are loaded into a flowset (fs).\n")
      
    } 
    
    # imports transformed data
    if(to_import == "b") {
      
      path_fcs = paste0(wd, "/", cellType, "/01 FCS files/01 transformed files")
      assign("trns", path_fcs, envir = .GlobalEnv)
      fs = read.flowSet(path=path_fcs,transformation=FALSE, truncate_max_range=FALSE)
      cleanFile = FALSE     
      assign("cleanFile", cleanFile, envir = .GlobalEnv)
      
      cat("> transformed FCS files have successfully been loaded into a flowset (fs).\n")
      
    }
    
    # imports cleaned and transformed data
    if(to_import == "c") {
      
      path_fcs = file.path(wd, cellType, "01 FCS files", "02 cleaned and transformed files", "Preprocessed", "QC", "PeacoQC_results", "fcs_files")
      assign("pqcd", path_fcs, envir = .GlobalEnv)
      fs = read.flowSet(path=path_fcs,transformation=FALSE, truncate_max_range=FALSE)
      cleanFile = TRUE  
      assign("cleanFile", cleanFile, envir = .GlobalEnv)
      
      cat("> Cleaned FCS files have successfully been loaded into a flowset.\n")
      
    }
    
    if(to_import == "d") {
      
      path_fcs = file.path(wd, cellType, "01 FCS files", "03 normalized files")
      assign("cn", path_fcs, envir = .GlobalEnv)
      fs = read.flowSet(path=path_fcs,transformation=FALSE, truncate_max_range=FALSE)
      
      cat("> Normalized FCS files have successfully been loaded into a flowset. \n")
    }
    
    
    
    # gives output if no a, b or c was the input for the function
    if (to_import != "a" & to_import != "b" & to_import != "c" & to_import != "d") {
      cat("No valid imput was was given for the 'to_import' variable.")
    }
    
    
    assign("fs", fs, envir = .GlobalEnv)
    assign("file_name", list.files(path_fcs), envir = .GlobalEnv)
    assign("markers", c(colnames(fs)), envir = .GlobalEnv)
    
    # print available channels on screen
      pData(parameters(fs[[1]]))[,1:2]
  }
  
  ## Chunk 3: cleaning and transformation of the data. ----
  
  listFiles = function(){
    
    ## listing reference files and compensation matrices
      assign("ref_file_name", c(list.files(rfs)), envir = .GlobalEnv)
      assign("comp_file_name", c(list.files(cms)), envir = .GlobalEnv)
      assign("ref_fs", read.flowSet(path=rfs,transformation=FALSE, truncate_max_range=FALSE), envir = .GlobalEnv)
    
    assign("cmf", paste0(wd, "/", cellType, "/", expID, " Catalyst ", cellType, " Comp_Ref_mapping_file.xlsx"), envir = .GlobalEnv)
    assign("rmf", paste0(wd, "/", cellType, "/", expID, " Catalyst ", cellType, " Reference_specification_file.xlsx"), envir = .GlobalEnv)
  }
  
  variableChecks = function(){
    
    stopifnot(cleaningMethod %in% c("PeacoQC", "None"))
    stopifnot(transformationMethod %in% c("Logicle", "None"))
    
    if(to_import %in% c("b", "c", "d") & transformationMethod == "Logicle"){stop("Data was already transformed, please alter input or put transformationMethod to None")}
    if(to_import %in% c("c", "d") & cleaningMethod == "PeacoQC"){stop("Data was already cleaned, please alter input or put cleaningMethod to None")}
    
    if (length(comp_file_name) > length(ref_file_name)){
      stop("There are more compensation matrices than reference files. You need at least 1 reference file for every compensation matrix")
    }
    
    # checking previous files
    if(cleaningMethod == "PeacoQC" & file.exists(file.path(wd, cellType, "01 FCS files","02 cleaned and transformed files"))){
       
       userInput = readline(prompt = "You already started and/or finnished cleaning before, do you want to override these files? (y/n) ")
        
        if(userInput == "n"){
          stop("Either input the cleaned files, or put cleaningMethod to None.")
         
      }
    }
    
    if(writePanel & file.exists(pnl)){
      userInput = readline(prompt = "Do you want to overwrite the existing panel file? (y/n) ")
      
      if(userInput == "n"){
        stop("Please manually remove the panel file, or put writePanel to None.")
      }
    }
    
    if(write_files){
      if (file.exists(cmf) | file.exists(rmf)){
        userInput = readline(prompt = "Compensation mapping file and/or reference mapping file already exist. Overwrite? (y/n) ")
        if (userInput == "n"){
          stop("Process aborted by user. Please set write_files to FALSE")
        } else if (userInput != "y"){
          stop ("Invalid user input")
        }
      }
    }
    
    if (multComp == FALSE & multRef == TRUE){
    if (length(comp_file_name) < 1 || length(ref_file_name) < 1) {
      stop("Please provide compensation matrices and reference files. Upload more files or set multComp and/or multRef to FALSE.")}
    }
    if (multComp == FALSE & multRef == TRUE){
      if (length(comp_file_name) > 1 | length(ref_file_name) < 1){
        stop("You supplied too many compensation matrices and/or you have not supplied enough reference files. Alternatively set multComp to TRUE and/or set multRef to FALSE")
      }
    if (multComp == FALSE & multRef == FALSE){
      if (length(comp_file_name) > 1 | length(ref_file_name) > 1){
        stop("You have supplied multiple compensation matrices and/or reference files. Please remove these, or set multComp and/or multRef to TRUE")}
      }
    }
  }
    
  fileNameCheck = function(){
    
    # makes sure the order of file names in file_name list matches order of file names in flowSet. If so the internal fcs filename is changed to the external file name in the file_name list
    filenameVars = c("TUBE NAME", "$FIL"
                     #, "$TBNM"
    )
    if(!is.na(customFileNameVar)){
      filenameVars[length(filenameVars)+1] = customFileNameVar
      assign("filenameVars", filenameVars, envir = .GlobalEnv)
    }
    
    for(i in seq_along(fs)){
      assign("fnCheck", TRUE, envir = .GlobalEnv)
      if(length(which(filenameVars %in% names(keyword(fs[[i]])))) == 0){
        userInput = readline(prompt = "Could not  match fcs file to external file name based on internal parameters. This could could lead to misidentification of fcs files. Do you want to skip this check? (y/n) ")
        if(userInput == "y"){
          message("File name check skipped. Adjusting internal fcs file name to external fcs file name")
          keyword(fs[[i]])$FILENAME = identifier(fs[[i]]) = file_name[i]
        }else if(userInput == "n"){
          stop("Please include 'TUBE NAME', '$FIL', or '$TBNM' in the file name when exporting from FlowJo. Alternatively enter one of the parameters used in the file name in the customFileNameVar variable above")
          cat("\n")
        }else{
          stop("Invalid argument")
        }
      }else{
        for(fnVar in which(filenameVars %in% names(keyword(fs[[i]])))){
          if(grepl(str_remove(keyword(fs[[i]])[[filenameVars[fnVar]]], ".fcs"), str_remove(file_name[i], ".fcs"))){
            keyword(fs[[i]])$FILENAME = identifier(fs[[i]]) = file_name[i]
            assign("fnCheck", TRUE, envir = .GlobalEnv)
            break
          }else{
            assign("fnCheck", FALSE, envir = .GlobalEnv)
          }
        }
        if(fnCheck == FALSE){
          userInput = readline(prompt = "Could not  match fcs file to external file name based on internal parameters. This could could lead to misidentification of fcs files. Do you want to skip this check? (y/n) ")
          if(userInput == "y"){
            message("File name check skipped. Adjusting internal fcs file name to external fcs file name")
            keyword(fs[[i]])$FILENAME = identifier(fs[[i]]) = file_name[i]
          }else if(userInput == "n"){
            stop("Please include 'TUBE NAME', '$FIL', or '$TBNM' in the file name when exporting from FlowJo. Alternatively enter one of the parameters used in the file name in the customFileNameVar variable above")
            cat("\n")
          }else{
            stop("Invalid argument")
          }
        }
      }
    }
    
    message("\nFile name order validated and internal file name is aligned.\n ")
  }
  
  writeMappingFiles = function(){
    if (write_files == TRUE){
      
      dateList = c()
      refDateList = c()
      
      #indexing dates of the sample files
      for (i in 1:length(fs)){
        dateList[i] = as.character(keyword(fs[[i]])$"$DATE")
      }
      
      #Indexing dates of the reference files
      for (i in 1:length(ref_fs)){
        t = as.character(keyword(ref_fs[[i]])$"$DATE")
        if (is.null(keyword(ref_fs[[i]])$"$DATE")){
          t = "NA"
        }
        refDateList[i] = t
      }
      
      if(multComp == TRUE & multRef == TRUE){
        write.xlsx(data.frame(index_number = 1:length(file_name), date = dateList, file_name = file_name, compensation_matrix = "example.csv", reference_file = "example.fcs"), cmf, row.names = FALSE)
        write.xlsx(data.frame(index_number = 1:length(ref_file_name), date = refDateList, file_name = ref_file_name, compensation_matrix = "example.csv"), rmf, row.names = FALSE)
        stop("Compensation and Reference Mapping File and Reference Specification File successfully written. You can now map your compensation matrices and reference files to your sample files and specify which compensation matrix to use for the reference file. When finished set write_files to FALSE and rerun the chunk")
      }else if (multComp == FALSE & multRef == TRUE){
        write.xlsx(data.frame(index_number = 1:length(file_name), date = dateList, file_name = file_name, compensation_matrix = comp_file_name, reference_file = "example.fcs"), cmf, row.names = FALSE)
        write.xlsx(data.frame(index_number = 1:length(ref_file_name), date = refDateList, file_name = ref_file_name, compensation_matrix = comp_file_name), rmf, row.names = FALSE)
        stop("Compensation and Reference Mapping File and Reference Specification File successfully written. You can now map your reference files to your sample files. When finished set write_files to FALSE and rerun the chunk")
      }else if (multComp == FALSE & multRef == FALSE){
        write.xlsx(data.frame(index_number = 1:length(file_name), date = dateList, file_name = file_name, compensation_matrix = comp_file_name, reference_file = ref_file_name), cmf, row.names = FALSE)
        write.xlsx(data.frame(index_number = 1:length(ref_file_name), date = as.character(refDateList), file_name = ref_file_name, compensation_matrix = comp_file_name), rmf, row.names = FALSE)
        message("Compensation and Reference Mapping File and Reference Specification File successfully listed Continuing with data cleaning.")
      }
    }
  }

  cleanAndTransform = function(){
    
    if (cleaningMethod == "None" & transformationMethod == "None"){
      message("Data has neither been cleaned, nor transformed as it has been performed previously.")
      cat(paste0("> Continuing with to_import = ", to_import, ".\n"))
      assign("internalPQCcheck", FALSE, envir = .GlobalEnv)
    } else{
      
      ##Map files and directories
      #set locations of comp_ref_mapping and reference_specification files
      assign("cmf", paste0(wd, "/", cellType, "/", expID, " Catalyst ", cellType, " Comp_Ref_mapping_file.xlsx"), envir = .GlobalEnv)
      assign("rmf" , paste0(wd, "/", cellType, "/", expID, " Catalyst ", cellType, " Reference_specification_file.xlsx"), envir = .GlobalEnv)
      
      #import the mapping files
      compMap = read.xlsx(cmf, 1, stringsAsFactors=FALSE)
      refCompMap = read.xlsx(rmf, 1, stringsAsFactors=FALSE)
      
      #create output directories
      setwd(paste0(wd, "/", cellType, "/01 FCS files"))
      dir.create("02 cleaned and transformed files", showWarnings = FALSE)
      assign("pqcd", paste0(wd, "/", cellType, "/01 FCS files/02 cleaned and transformed files"), envir = .GlobalEnv)
      assign("dir_prepr", paste0(pqcd, "/Preprocessed/"), envir = .GlobalEnv)
      assign("dir_QC", paste0(pqcd, "/Preprocessed/QC"), envir = .GlobalEnv)
      assign("dir_RDS", paste0(pqcd, "/RDS"), envir = .GlobalEnv)
      assign("dir_results", paste0(pqcd, "/ResultsPeacoQC"), envir = .GlobalEnv)
      assign("dir_raw", paste0(fcs, "/"), envir = .GlobalEnv)
      
      for (path in c(dir_prepr, dir_QC, dir_RDS, dir_results)){
        dir.create(path, showWarnings = FALSE)
      }
      
      if(transformationMethod == "Logicle"){ dir.create(trns, recursive = TRUE, showWarnings = FALSE) }
      
      #list the files for import
      files <- list.files(path = dir_raw, pattern = file_pattern)
      
      #create lists for statistics
      preRM = c()
      postRM = c()
      percRM = c()
      binList = c()
      evtsList = c()
      itr2 = 0
      
      
      ##Transformation and/or cleaning
      
      #Making sure minimum number of bins is filled for PeacoQC
      if(cleaningMethod == "PeacoQC"){
        cellNo = c()
        for (i in 1:length(fs)){
          cellNo[i] = keyword(fs[[i]])$"$TOT"
        }
        assign("cellNo", cellNo, envir = .GlobalEnv)
        
        #If minimum number of bins is not filled the amount of events per bin is lowered
        if (ceiling((as.integer(cellNo[which.min(cellNo)])/dEvtsPerBin)) < (minBin-1)){
          dEvtsPerBin = ceiling(as.integer(cellNo[which.min(cellNo)])/minBin)
          cat("\n")
          print(paste0("The number of events per bin has been lowered to ", dEvtsPerBin))
          cat("\n")
        }
      }
      
      for (file in files){
        print(paste0("Current file: ", file))
        itr2 = itr2 + 1
        
        ##import compensation  matrix
        compName = as.character(compMap$compensation_matrix[which(grepl(file, compMap$file_name))])
        compensation_matrix <- read.csv(paste0(cms, "/", compName), check.names = FALSE, row.names = 1) 
        colnames(compensation_matrix) <- sub(" :: .*", "", colnames(compensation_matrix))
        print(paste0("Current compensation matrix: ", compName))
        
        ##Process reference file
        rfile = compMap$reference_file[which(grepl(file, compMap$file_name))]
        print (paste0("Current reference file: ", rfile))
        reference_file <- read.FCS(paste0(rfs, "/", as.character(rfile)), truncate_max_range = FALSE)
        
        
        
        #Check if flow data is conventional, or spectral and removeMargins/compensate reference file
        if(ds == "FC"){
          if(cleaningMethod == "PeacoQC"){
            channels_of_interest <- GetChannels(object = reference_file, markers = markers_of_interest, exact = FALSE)
            ff_m <- PeacoQC::RemoveMargins(reference_file, channels_of_interest)
            ff_c <- flowCore::compensate(ff_m, compensation_matrix)
          }else if(cleaningMethod == "None")  
            ff_c <- flowCore::compensate(reference_file, compensation_matrix)
        }
        
        if(ds == "SFC"){
          if (cleaningMethod == "PeacoQC"){
            warning("Due to often occuring issues with spectral data RemoveMargins is not applied to the reference file")
            channels_of_interest <- GetChannels(object = reference_file, markers = markers_of_interest, exact = FALSE)
          }
          ff_c <- flowCore::compensate(reference_file, compensation_matrix)
        }
        
        if(ds != "FC" & ds != "SFC"){
          stop("Imput for ds in chunk 1 is not valid! Either choose FC or SFC!")
        }
        
        
        #apply Logicle transformation to the reference file
        tryCatchLogicle = tryCatch(estimateLogicle(ff_c, colnames(compensation_matrix)), error = function(e)e)
        if(inherits(tryCatchLogicle, "error")){
          customM = 4.5
          for(customMitr in 1:20){
            tryCatchLogicle = tryCatch(estimateLogicle(ff_c, colnames(compensation_matrix), m = customM), error = function(e)e)
            if(inherits(tryCatchLogicle, "error")){
              customM = customM + 0.1
            }else{
              break
            }
          }
          if(customMitr == 20){
            stop("Couldn't transform data")
          }else{
            message(paste0("Transforming data using m = ", customM))
            translist <- estimateLogicle(ff_c, colnames(compensation_matrix), m = customM)
          }
        }else{
          translist <- estimateLogicle(ff_c, colnames(compensation_matrix))
        }
        
        ff_t <- flowCore::transform(ff_c, translist)
        cat("\n")
        print("Reference file successfully processed")
        
        
        ##process sample files
        #Check if compensation matrix assigned to sample file matches with matrix asigned to reference file
        if (compName != as.character(refCompMap$compensation_matrix[which(grepl(rfile, refCompMap$file_name))])){
          stop (paste0( "The compensation matrix asigned to ", rfile, " does not match the compensation matrix asigned to ", file, ". Either change the mapping of the compensation matrix, or make sure the reference file is correctly mapped"))
        }
        
        #import sample file
        ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE)
        
        #removes margins if cleaning is TRUE
        if(cleaningMethod == "PeacoQC"){
          preRM[itr2] = keyword(ff)$"$TOT"
          ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)
          postRM[itr2] = keyword(ff_m)$"$TOT"
          percRM[itr2] = (((as.integer(preRM[itr2])-as.integer(postRM[itr2]))/as.integer(preRM[itr2]))*100)
          
          #applies compensation
          ff_c <- flowCore::compensate(ff_m, compensation_matrix)
        }else if(cleaningMethod == "None"){
          ff_c <- flowCore::compensate(ff, compensation_matrix)
        }
        
        #applies transformation
        ff_t <- flowCore::transform(ff_c, translist)
        #selected_live <- filter(ff_s, live_gate)
        #ff_l <- ff_s[selected_live@subSet, ]
        ff_l = ff_t
        
        #cleans data with PeacoQC
        if(cleaningMethod == "PeacoQC"){    
          evtsPerBin = dEvtsPerBin
          #checks if minimum number of bins can still be filled after RemoveMargins
          if (ceiling(as.integer(postRM[itr2])/evtsPerBin) < (minBin-1) & itr2 != 1){
            stop(paste0("The RemoveMargins function removed so many events that the minimum number of bins cannot be filled. Please set the dEvtsPerBin value to: ", round(as.integer(postRM[itr2])/minBin, digits = 0), ". Alternatively, lower the value of minBin."))
          } else if (ceiling(as.integer(postRM[itr2])/evtsPerBin) < (minBin-1) & itr2 == 1){
            evtsPerBin = ceiling(postRM[itr2]/minBin)
          }
          
          #Checks if maximum number of bins is not exceeded
          if (as.integer(postRM[itr2])/evtsPerBin > maxBin){
            evtsPerBin = ceiling(as.integer(postRM[itr2])/(maxBin-1))
            cat("\n")
            print(paste0("The maximum amount of bins is exceeded for ", file, ". Increasing amount of events per bin to: ", evtsPerBin))
            cat("\n")
          }       
          
          binList[itr2] = ceiling(as.integer(postRM[itr2])/evtsPerBin)
          evtsList[itr2] = evtsPerBin
          assign("PQC", PeacoQC::PeacoQC(ff = ff_l, channels = channels_of_interest, plot = plotPQC, save_fcs = TRUE, output_directory = dir_QC, MAD = MADValue, IT_limit = ITLimitValue, events_per_bin = evtsPerBin), envir = .GlobalEnv)
          
          #write FCS files
          write.FCS(PQC$FinalFF, file = paste0(dir_prepr, file))
        }else if(cleaningMethod == "None"){
          write.FCS(ff_l, file = paste0(trns, "/", file))
        }
      }
      
      if(cleaningMethod == "PeacoQC"){
        assign("rms", paste0(pqcd, "/PeacoQC metadata.xlsx"), envir = .GlobalEnv)
        write.xlsx(data.frame(index_values = 1:length(files), file_name = files, pre_remove_margins = preRM, post_remove_margins = postRM, percentage_removed = percRM, events_per_bin = evtsList, number_of_bins = binList), rms, row.names = FALSE)
        pdf(file = paste0(dir_results, "/PeacoQC Heatmap.pdf"), width = 10, height = 30)
        PeacoQCHeatmap(paste0(pqcd,"/Preprocessed/QC/PeacoQC_results/PeacoQC_report.txt"))
        dev.off()
        
        write.xlsx(read.delim(paste0(pqcd,"/Preprocessed/QC/PeacoQC_results/PeacoQC_report.txt")), paste0(pqcd,"/Preprocessed/QC/PeacoQC_results/PeacoQC_report.xlsx"), row.names = FALSE)
        
        assign("file_name", c(list.files(paste0(pqcd,"/Preprocessed/QC/PeacoQC_results/fcs_files"))), envir = .GlobalEnv)
        assign("fs_clean", read.flowSet(path=paste0(pqcd,"/Preprocessed/QC/PeacoQC_results/fcs_files"),transformation=FALSE, truncate_max_range=FALSE), envir = .GlobalEnv)
        assign("fs_raw", fs, envir = .GlobalEnv)
        assign("fs", fs_clean, envir = .GlobalEnv)
        assign("markers", c(colnames(fs)), envir = .GlobalEnv)
        assign("internalPQCcheck", TRUE, envir = .GlobalEnv)
        
        cat("\n")
        print("Data successfully cleaned by PeacoQC and transformed using Logicle transformation")
        cat("\n")
        print("Cleaned and transformed FCS files have successfully been loaded into a flowset. Ready for further analysis. The raw data is stored in fs_raw")
        
        to_import = "c"
        assign("to_import", to_import, envir = .GlobalEnv)
      }else if(cleaningMethod == "None"){
        assign("fs_clean", read.flowSet(path=trns,transformation=FALSE, truncate_max_range=FALSE), envir = .GlobalEnv)
        assign("fs", fs_clean, envir = .GlobalEnv)
        assign("markers", c(colnames(fs)), envir = .GlobalEnv)
        assign("internalNONEcheck", TRUE, envir = .GlobalEnv)
        
        cat("\n")
        cat("Data successfully transformed using Logicle transformation. No cleaning was performed.\n\n")
        cat("Transformed FCS files have successfully been loaded into a flowset. Ready for further analysis. The raw data is stored in fs_raw.\n\n")
        
        to_import = "b"
        assign("to_import", to_import, envir = .GlobalEnv)
      }
      
    }
  }
  
  writePanelFile = function(){
    if (writePanel == TRUE){
      tempMarkers = as.character(pData(parameters(fs[[1]]))$desc)
      if(length(which(is.na(as.character(pData(parameters(fs[[1]]))$desc)))) > 0){
        for(i in which(is.na(as.character(pData(parameters(fs[[1]]))$desc)))){
          tempMarkers[i] = pData(parameters(fs[[1]]))$name[i]
        }
      }
      write.xlsx(data.frame(index_values = 1:length(markers), fcs_colname = markers, antigen=tempMarkers, marker_class=vector(length=length(markers))), pnl, row.names = FALSE)
      message ("Panel file template successfully written. You can now indicate which markers should be used for cluster analysis by writing 'type' in the 'marker_class' column. Markers indicated as 'state' will not be taken along for clustering, but the relative expression will still be analysed. Markers indicted as 'none' will not be included in any analyses.\n")
    }
  }
  
  ## Chunk 4: validation panel file and matching internal and external fcs file names.----
  
  panelCheck = function(){
    
    panel = read.xlsx(pnl,1, stringsAsFactors=FALSE) 
    assign("panel", panel, envir = .GlobalEnv)
    
    
    if (length(panel$fcs_colname) == length(markers)){##checks if length of supplied panel matches fcs file panel
      if (all(panel$fcs_colname %in% markers) == TRUE){##checks if markers in supplied channel are all present in fcs file panel
        mismatch_position = c()
        target_position = c()
        for (i in 1:length(markers)){
          if (panel$fcs_colname[i] == markers [i]){##checks if position of supplied channels match position of fcs file channels
            target_position = append(target_position, (i/10))##position is stored
          } else {##if there is a position mismatch between suppled and fcs file channels...
            mismatch_position = append(mismatch_position, i)##...position of mismatch is stored...
            target_position = append(target_position, ((which(markers == panel$fcs_colname[i]))/10))##...and the target position is stored
          }
        }
        if (length(mismatch_position) == 0){
          message("Panel match successfully validated\n")
        } else {##if there are any position mismatches...
          panel$Target = target_position##...the target positions are appended to the dataframe...
          panel = arrange(panel, Target)##...and the dataframe is adjusted so that the positions match with the fcs file...
          assign(panel, subset(panel, select = -c(Target)), envir = .GlobalEnv)
          message("Panel match successfully validated. Panel order has been modified to match fcs file\n")
        }
      } else {##if the supplied channels do not match the fcs file channels...
        mismatch=c(panel$fcs_colname %in% markers)##...a list is made of matches and mismatches
        mismatch_channel = c(which(mismatch == "FALSE"))##positions of mismatches are determined
        cat("supplied panel does not match fcs file panel\n")
        cat(length(mismatch_channel), "conflict(s) found:/n")
        cat("\n")
        conflict_number=1
        cat("channels in FCS file:\n")
        print(markers)
        cat("\n")
        for (i in mismatch_channel){
          cat("conflict", conflict_number, ":\n")
          cat("fcs file panel:", markers[i], "\n")##prints the channel from the fcs file
          cat("supplied panel:", panel$fcs_colname[i], "\n")##prints the incorrect channel from the supplied file
          cat("\n")
          conflict_number = conflict_number + 1
        }
        stop()
      }
    } else {
      stop("Number of supplied channels does not match number of channels in fcs file")
    }
  }
  
  writeMetaDataFile = function(){
    if (writeMeta == TRUE){
      write.xlsx(data.frame(file_name, sample_id=vector(length=length(file_name)), condition=vector(length=length(file_name)), patient_id=vector(length=length(file_name)), sample_control=vector(length=length(file_name)), batch=vector(length=length(file_name))), mtd, row.names = FALSE)
      message("Metadata file successfully written. You can now add metadata to this file. Make sure every file has a unique SampleID! You can also add columns with extra data (e.g. Gender, Age, or Timepoint")
      cat("\n")
    }else{
      message("Metadata file imported.")
    }
  }
  
  
  ## Chunk 5: validation metadata and generating diagnostic plots. ----
  
  specifyChannels = function(){
    markers_CF = c()
    for (i in markers){
      ind = which(grepl(i, panel$fcs_colname))
      if (panel$marker_class[ind] == "type" | panel$marker_class[ind] == "state"){
        markers_CF = append(markers_CF, markers[ind])
      }
    }
    assign("markers_CF", markers_CF, envir = .GlobalEnv)
  }
  
  updateMarkerNames = function(){
    assign("marker_list", panel$antigen, envir = .GlobalEnv)
    assign("channels", c(colnames(fs)), envir = .GlobalEnv)
    names(marker_list) = channels 
    markernames(fs) = marker_list
  }
  
  importData = function(){
    
    if(to_import %in% c("b", "c")){    assign("md", read.xlsx(mtd, 1, stringsAsFactors=FALSE), envir = .GlobalEnv)    }
    if(to_import == "d"){ 
      md <- read.xlsx(file.path(wd, cellType, paste(expID, "Catalyst", cellType, "metadata CN.xlsx")), 1, stringsAsFactors=FALSE)
      md <- dplyr::rename(md, file_name = file_name2)
  
      assign("md", md, envir = .GlobalEnv)
    }
    

    #import metadata
    assign("conList", c(as.character(unique(md$condition))), envir = .GlobalEnv)
    
    
    if(TRUE %in% c(is.na(md))){
      stop("The metadata file contains empty cells")
    }
    
    for (i in seq_along(md$file_name)){
      if (file_name[i] != md$file_name[i]){
        stop("the file names in the metadata file do not match the file names in the file directory and/or flowSet")
      }
    }
    
    #All metadata categories, except for file name and sampleID, are stored in a list to be passed to the prepData function below
    assign("factors_list", c(colnames(md)[-c(1:2)]), envir = .GlobalEnv)
    
    #metadata categories are converted from character to factor
    for (i in factors_list){
    md[, i] = as.factor(md[, i])
    }
    
    ##import panel
    assign("panel", read.xlsx(pnl, 1, stringsAsFactors=FALSE), envir = .GlobalEnv)
    
    
    ##fixes bug where not all selected channels were taken along for analysis (for raw FlowSet)
    for (i in seq_along(fs)){
      keyword(fs[[i]])[["$CYT"]] <- "FACS"
    }
    
    assign("markers", c(colnames(fs)), envir = .GlobalEnv)

  }
  
  TCcheck = function(){
    if(TCs == TRUE){
      if("control" %in% md$sample_control){
        message("Controls detected")
      }else if("Control" %in% md$sample_control){
        TCindexList = which(md$sample_control == "Control")
        for(TCindex in TCindexList){
          md$sample_control[TCindex] = "control"
        }
        message("Changed Technical control label from 'Control' to 'control'")
        cat("\n")
      }else{
        stop("Could not detect Technical Control samples. Please set 'TCs' to FALSE, or make sure you labeled the Technical Controls as çontrol' in the sample_control column of the metadata")
      }  
    }
  }
  
  createSCE = function(){
    sce = prepData(fs, panel = panel, md = md, features=markers_CF, FACS = TRUE, transform = FALSE, md_cols = list(file = "file_name", id = "sample_id", factors = factors_list))
    assayNames(sce)[1] = "exprs" ##when transform = FALSE in prepData the "exprs" assay is not made. The raw data is stored in the "counts" assay. If the input flowSet has already been transformed by FlowVS the "count" assay can be renamed to "exprs". Downstream analyses only accept the "exprs" assay.
    message("Single Cell Experiment successfully generated")
    cat("\n")
    assign("sce", sce, envir = .GlobalEnv)
  }
  
  previousDiagnostics = function(){
    
    # create folder to store them in
    if(to_import == "b"){folder_name = "transformed data"}
    if(to_import == "c"){folder_name = "cleaned and transformed data"}
    if(to_import != "b" & to_import != "c"){stop("Invalid input for variable to_input for this step.")}
    
    folder_path = file.path(paste0(rsl, "/00 Diagnostic plots of ", folder_name))
    assign("folder_path", folder_path, envir = .GlobalEnv)
    
    # when there is already a folder
    if (file.exists(folder_path) & write_files == TRUE) {
      
      continue <- readline("There is already a folder for diagnostics, do you want to overwrite? Note, you will remove the previous files! (y/n): ")
      
      if (tolower(continue) == "n") {
        stop("Please put write_files to FALSE!")
      }
    }
  }
  
  writeDiagnostics = function(){
    
    if (write_files == TRUE){
      
      message("Plots sucessfully created:")
      
      # Create the folder/directory
      if (!file.exists(folder_path)) {dir.create(folder_path, showWarnings = FALSE)}
      
      #writes pdf file containing bar graph of cell counts
      pdf(file = paste0(folder_path, "/Counts ", cellType, ".pdf"), width = 30, height = 10)
          print(plotCounts(sce, color_by = "sample_id") +
                  ggtitle("Count summary of everything in total."))
      
          if(length(conList)>1){
            for (i in 1:length(conList)){
              abc = conList[i]
              print(plotCounts(filterSCE(sce, condition == abc), color_by = "sample_id")+
                      ggtitle(paste0("Plot for: ", abc)))
              }
          }else{
            print(plotCounts(sce, color_by = "sample_id"))
          }  
      dev.off()
      
      cat("> Counts plot.\n")
      
      if(TCs == TRUE){
        
        #writes pdf file containing MDS plot of technical controls
        if(length(which(md$sample_control == "control"))>2){
          pdf(file = paste0(folder_path, "/MDS ", cellType, " TCs.pdf"), width = 10, height = 10)
          print(pbMDS(filterSCE(sce, sample_control == "control"), by = c("sample_id"), fun = "median", features = NULL, assay = "exprs", label_by = NULL, color_by = "sample_id"))
          dev.off()
          cat("> MDS TC plot.\n")
        }
      
        # writes expression plot of technical controls  
        pdf(file = paste0(folder_path, "/Expression ", cellType, " TCs.pdf"), width = 30, height = 10)
        print(plotExprs(filterSCE(sce, sample_control == "control"), color_by = "sample_id"))
        dev.off()
        cat("> Expression TC plot.\n\n")
      }
      
      #writes pdf file containing histogram of transformed expression
      pdf(file = paste0(folder_path, "/Expression ", cellType, ".pdf"), width = 30, height = 10)
      for (i in 1:length(conList)){
        abc = conList[i]
        print(plotExprs(filterSCE(sce, condition == abc), color_by = "sample_id")+
                ggtitle(paste0("Plot for: ", abc)))
      }
      dev.off()
      cat("> Expression plot.\n")
    }
    
    assign("cellNoTot", sum(n_cells(sce)), envir = .GlobalEnv)
    
    if(TCs == TRUE){
      assign("sce_c", (filterSCE(sce, sample_control == "control")), envir = .GlobalEnv)
      assign("sce", filterSCE(sce, sample_control != "control"), envir = .GlobalEnv)
      cat("\nControl group successfully removed from data\n")
      
    }
    
    assign("cellNoNC", sum(n_cells(sce)), envir = .GlobalEnv)
    
    if (write_files == TRUE){
      
      #writes pdf file containing heatmap of relative marker expression
      pdf(file = paste0(folder_path, "/Expression_heatmap ", cellType, ".pdf"), width = 10, height = 10)
      print(plotExprHeatmap(sce, bin_anno = FALSE, row_anno = exprHMRowAnno, row_clust=TRUE, by=c("sample_id"), fun = c("median")))
      dev.off()
      cat("\n> Expression heatmap containing relative marker expressions.\n")
      
      #writes pdf file containing non-redundancy scores, allowing for useful marker selection for downstream analysis
      pdf(file = paste0(folder_path, "/Non-redundancy_score ", cellType, ".pdf"), width = 10, height = 10)
      print(plotNRS(sce, features = NULL, color_by = NRSColorBy, assay = "exprs"))
      dev.off()
      cat("> Heatmap and NRS plots.\n")
      
      #writes pdf file containing non-redundancy scores, allowing for useful marker selection for downstream analysis
      pdf(file = paste0(folder_path, "/MDS ", cellType, ".pdf"), width = 10, height = 10)
      print(pbMDS(sce, by = "sample_id", fun = "median", features = NULL, assay = "exprs", label_by = NULL, color_by = MDSColorBy))
      dev.off()
      cat("> MDS plot.\n")
    }
    
  }
  
  ## Chunk 6: preparing files and directories for cross entropy test. ----
  # functions created based on the code from: https://github.com/AdrianListon/Cross-Entropy-test/tree/main/flow%20analysis, to run CE test more automatically. 
  
  checkInputMarkers = function(){
    
      # remove side/foward scatter variables
      if(sum(channels.of.interest %in% c("SSC", "FSC", "SSC-A", "FSC-A", "SSC-H", "FSC-H", "SSC-W", "FSC-W")) > 0 ) {
        
        message("One of the following markers was detected and will be removed: SSC, FSC, SSC-A, FSC-A, SSC-H, FSC-H, SSC-W, FSC-W.\n")
        channels.of.interest <- markers_of_interest[!markers_of_interest %in% c("SSC", "FSC", "SSC-A", "FSC-A", "SSC-H", "FSC-H", "SSC-W", "FSC-W")]
        
        assign("channels.of.interest", channels.of.interest, envir = .GlobalEnv)
      }
      
        # Display markers of interest
        cat("These markers were noted as interesting before:\n", channels.of.interest, "\n")
        
        # Ask user if they want to continue with these markers

        continue <- readline("Do you want to continue with these markers? (y/n): ")
        
        if (tolower(continue) == "n") {
          stop("Please manually alter the 'markers_of_interest' variable on top of the chunk and run the chunk again.")
        }
    }
  
  prepCE = function (filesCE, md_column) {
    
    if(filesCE == "a"){stop("You did not transform your data or did not import it correctly. Go back to rmd file 1.")}
    if(filesCE == "b"){ imported_data = "trns"}
    if(filesCE == "c"){ imported_data = "pqcd"}
    if(filesCE == "d"){ imported_data = "cn"}
    
    assign("imported_data", imported_data, envir = .GlobalEnv)
    assign("to_compare", md_column, envir = .GlobalEnv)
    
    if (subsettingCE){folder_name_CE = paste0("Files: ", imported_data, " (comparing: ", to_compare,", seed: ", as.character(fcs.seed.base), ", subsetted column: ", subset_column, ", subseted for: ", subset_for, ")")}
    if(!subsettingCE){folder_name_CE = paste0(imported_data, " (", to_compare,", seed of ", as.character(fcs.seed.base), ")")}

    
    # checking folder structure, creating it if necessary
    fcs_folder = paste0(rsl, "/01 Cross_Entropy/", folder_name_CE, "/renamed FCS files")
    ce_folder = paste0(rsl, "/01 Cross_Entropy/", folder_name_CE)
    
    assign("fcs_folder", fcs_folder, envir = .GlobalEnv)
    assign("ce_folder", ce_folder, envir = .GlobalEnv)
    
    
    if ( file.exists(ce_folder)){
      
      continue <- readline("There is already a folder for CE with this input, do you want to overwrite? Note, you will remove the previous files! (y/n): ")
      
      if (tolower(continue) == "n") {
        stop("Please put crossEntropyTest to FALSE!")
      }
    }
    
    unlink(ce_folder, recursive = TRUE) # remove files if there
    
    if ( ! file.exists( fcs_folder )) {
      dir.create( fcs_folder, recursive = TRUE , showWarnings = FALSE)
      }
    
    message("\nFolder structure created.")
    
    # link directories
    source_directory = switch(filesCE,
      "b" = trns,
      "c" = paste0(wd, "/", cellType, "/01 FCS files/02 cleaned and transformed files/Preprocessed/QC/PeacoQC_results/fcs_files"),
      "d" = paste0(wd, "/", cellType, "/01 FCS files/03 normalized files")
    )
    
    destination_directory = fcs_folder
    
    message("Directories correctly linked.")
    
    assign("source_directory", source_directory, envir = .GlobalEnv)
    assign("destination_directory", destination_directory, envir = .GlobalEnv)
    
    # subsetting the data
    if(subsettingCE){
      subsetfiles <- subset(md, md[[subset_column]] == subset_for)
      file_names_to_copy = subsetfiles$file_name

      files_to_copy <- paste(source_directory, file_names_to_copy, sep = "/")
      
    }else{files_to_copy = list.files(source_directory, full.names = TRUE)}
    
    # code to copy files from source to destination directory
      
      
      # actually copying 
      for (file in files_to_copy){
        destination_path <- file.path(destination_directory, basename(file))
        file.copy(file, destination_path)
      }

      message("Files copied to CE folder.")
    
    # rename all the files based on the metadata column specified
      
      # List all files in the destination directory
      files <- list.files(destination_directory)
      
      # get the right metadata file
      md_CE = switch(to_import,
                                "b" = md,
                                "c" = md,
                                "d" = md_2)
      
      assign("md_CE", md_CE, envir = .GlobalEnv)
      
      
      # rename files
      prefixes_CE <- c()  # Initialize an empty vector to store unique prefixes
      
      for (file in files) {
        
        file_info <- md_CE[md_CE$file_name == file, md_column]
        prefix <- paste0(md_column, "_", file_info, "_")
        
        # Check if the prefix is already in prefixes_CE
        if (!(prefix %in% prefixes_CE)) {
          # If not, append it to the vector
          prefixes_CE <- unique(c(prefixes_CE, prefix))
          
        }
        
        new_file_name <- paste0(prefix, file)
        
        # Construct full paths for old and new file names
        old_file_path <- file.path(fcs_folder, file)
        new_file_path <- file.path(fcs_folder, new_file_name)
        
        # Rename the file
        file.rename(old_file_path, new_file_path)
      }
      assign("prefixes_CE", prefixes_CE, envir = .GlobalEnv)
      
      
      
      message("Files successfully renamed.")
      
      
    }
  
  prepVarCE = function() { 

    assign("fcs.dmrd.data.sample.n.per.sample", NULL, envir = .GlobalEnv)
    assign("fcs.dmrd.data.sample.n", NULL, envir = .GlobalEnv)
    
    
    # link all the files and directories
    fcs.data.dir <- destination_directory
    assign("fcs.data.dir", fcs.data.dir, envir = .GlobalEnv)
    files <- list.files(rfs) # Get the list of files in the directory
    
    # save prefixes created before in another variable
    fcs.condition = prefixes_CE
    assign("fcs.condition", fcs.condition, envir = .GlobalEnv)
    
    fcs.condition.label <- setNames(fcs.condition, fcs.condition)
    assign("fcs.condition.label", fcs.condition.label, envir = .GlobalEnv)
    
    flowFrame <- read.FCS(list.files(fcs.data.dir, "\\.fcs$", full.names = TRUE)[1], truncate_max_range = FALSE)
    assign("flowFrame", flowFrame, envir = .GlobalEnv)
    
    channels <- data.frame(name = unname(pData(parameters(flowFrame))$name), 
                           desc = unname(pData(parameters(flowFrame))$desc))
    assign("channels", channels, envir = .GlobalEnv)
    
    for(i in 1:dim(channels)[1] ){
      channels$out1[i] = paste0("\"", channels$name[i], "\", #", channels$desc[i], "\n")
      channels$out2[i] = paste0("\"", channels$name[i], "\" = \"", channels$desc[i], "\",\n") 
      channels$out3[i] = paste0("\"", channels$name[i], "\" = 200, #", channels$desc[i], "\n") 
    }
    
    descs.filtered <- channels$desc[!is.na(channels$desc) & channels$desc != '-'] 
    channels.filtered <- filter(channels, desc %in% descs.filtered)
    assign("descs.filtered", descs.filtered, envir = .GlobalEnv)
    assign("channels.filtered", channels.filtered, envir = .GlobalEnv)
    
    temp <- channels$out1[channels$desc %in% channels.of.interest]
    temp[length(temp)] <- gsub(',', '', temp[length(temp)])
    
    temp <- channels$out2[channels$desc %in% channels.of.interest]
    temp[length(temp)] <- gsub(',', '', temp[length(temp)])
    
    assign("temp", temp, envir = .GlobalEnv)
    
    
    # Initialize an empty list to store the associations
    fcs_channels <- list()
    
    # Loop through each row in the panel data frame
    for (i in 1:nrow(panel)) {
      marker <- panel$antigen[i]
      fcs_colname <- panel$fcs_colname[i]
      
      # Check if the marker is in the markers of interest
      if (marker %in% channels.of.interest) {
        fcs_channels[[marker]] <- fcs_colname
      }
    }
    
    
    
    fcs.channel = unname(unlist(fcs_channels))
    assign("fcs.channel", fcs.channel, envir = .GlobalEnv)
    
    # add the output to fcs.channel.label in the parameter file
    temp <- channels$out2[channels$desc %in% channels.of.interest]
    temp[length(temp)] <- gsub(',', '', temp[length(temp)])
    
    message("\nCopy the output on screen, and add to the fcs.channel.label variable in the next chunk!")
    cat(paste(temp, collapse = ""))
    }
  
  ## Chunk 7: running the cross entropy test. ----
  # functions created based on the code from: https://github.com/AdrianListon/Cross-Entropy-test/tree/main/flow%20analysis, to run CE test more automatically. 
  
  run_CE = function() {
    
    assign("fcs.use.cached.results", FALSE, envir = .GlobalEnv)
    setwd(paste0(rsl, "/"))
    
    # analyze_flow_cytometry_parameter.r ----
    
    library( ConsensusClusterPlus )
    library( digest )
    require( dunn.test )
    library( flowCore )
    library( FlowSOM )
    library( ggplot2 )
    library( ggridges )
    require( RANN )
    library( RColorBrewer )
    library( reshape2 )
    library( Rtsne )
    library( umap )
    
    fcs.condition.n <- length( fcs.condition )
    fcs.channel.n <- length( fcs.channel )
    
    assign("fcs.condition.n", fcs.condition.n, envir = .GlobalEnv)
    assign("fcs.channel.n", fcs.channel.n, envir = .GlobalEnv)
    
    
    fcs.sample.number.width <- 2
    fcs.event.number.width <- 6
    
    assign("fcs.sample.number.width", fcs.sample.number.width, envir = .GlobalEnv)
    assign("fcs.event.number.width", fcs.event.number.width, envir = .GlobalEnv)
    
    folder_name = ce_folder
    
    # graphics parameters
    
    fcs.color.pool <- c( 
      brewer.pal( 8, "Set1" )[ -6 ], 
      brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
      adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
                   red.f = 0.9, green.f = 0.8, blue.f = 0.7 ), 
      adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
                   red.f = 0.9, green.f = 0.8, blue.f = 0.7 ), 
      adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
                   red.f = 0.8, green.f = 0.6, blue.f = 0.5 ), 
      adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
                   red.f = 0.8, green.f = 0.6, blue.f = 0.5 ), 
      adjustcolor( brewer.pal( 8, "Set1" )[ -6 ], 
                   red.f = 0.3, green.f = 0.3, blue.f = 0.3 ), 
      adjustcolor( brewer.pal( 7, "Set2" )[ c( 1, 3, 6 ) ], 
                   red.f = 0.3, green.f = 0.3, blue.f = 0.3 ) )
    fcs.color.pool.n <- length( fcs.color.pool )
    
    assign("fcs.color.pool", fcs.color.pool, envir = .GlobalEnv)
    
    
    fcs.line.type.pool <- 1:6
    fcs.line.type.pool.n <- length( fcs.line.type.pool )
    
    assign("fcs.line.type.pool", fcs.line.type.pool, envir = .GlobalEnv)
    assign("fcs.line.type.pool.n", fcs.line.type.pool.n, envir = .GlobalEnv)
    
    
    fcs.condition.color <- rep( 
      fcs.color.pool, 
      ceiling( fcs.condition.n / fcs.color.pool.n ) 
    )[ 1 : fcs.condition.n ]
    names( fcs.condition.color ) <- fcs.condition
    
    fcs.condition.line.type <- rep( 
      fcs.line.type.pool, 
      ceiling( fcs.condition.n / fcs.line.type.pool.n ) 
    )[ 1 : fcs.condition.n ]
    names( fcs.condition.line.type ) <- fcs.condition
    
    assign("fcs.condition.line.type", fcs.condition.line.type, envir = .GlobalEnv)
    assign("fcs.condition.color", fcs.condition.color, envir = .GlobalEnv)
    
    
    # density parameters
    # Once done, set to NULL and re-run the script on all the data.
    
    fcs.density.data.sample.n <- NULL
    
    assign("fcs.density.data.sample.n", fcs.density.data.sample.n, envir = .GlobalEnv)
    
    fcs.density.partition.all <- "all"
    fcs.density.partition.all.label <- c( "all" = "All" )
    fcs.density.partition.all.color <- c( "all" = "grey" )
    
    assign("fcs.density.partition.all", fcs.density.partition.all, envir = .GlobalEnv)
    assign("fcs.density.partition.all.label", fcs.density.partition.all.label, envir = .GlobalEnv)
    assign("fcs.density.partition.all.color", fcs.density.partition.all.color, envir = .GlobalEnv)
    
    
    fcs.density.font.size <- 4.5
    assign("fcs.density.font.size", fcs.density.font.size, envir = .GlobalEnv)
    
    
    fcs.density.line.size <- 0.2
    fcs.density.line.alpha <- 0.3
    
    assign("fcs.density.line.size", fcs.density.line.size, envir = .GlobalEnv)
    assign("fcs.density.line.alpha", fcs.density.line.alpha, envir = .GlobalEnv)
    
    
    fcs.density.figure.width.base <- 0.4
    fcs.density.figure.height.base <- 0.1
    
    assign("fcs.density.figure.width.base", fcs.density.figure.width.base, envir = .GlobalEnv)
    assign("fcs.density.figure.height.base", fcs.density.figure.height.base, envir = .GlobalEnv)
    
    
    fcs.density.figure.dir <- "./figure_density"
    
    assign("fcs.density.figure.dir", fcs.density.figure.dir, envir = .GlobalEnv)
    
    
    fcs.density.figure.sample <- "density_sample"
    fcs.density.figure.cluster <- "density_cluster"
    
    assign("fcs.density.figure.sample", fcs.density.figure.sample, envir = .GlobalEnv)
    assign("fcs.density.figure.cluster", fcs.density.figure.cluster, envir = .GlobalEnv)
    
    
    # cluster parameters
    # Tip: naming clusters can be done automatically as clusters 1:n via the first command below.
    # To rename your clusters based on marker expression, insert a # before the first command,
    # remove the # before the second command,
    # and rename the clusters appropriately.
    fcs.cluster <- sprintf( "%02d", 1 : fcs.cluster.n )
    
    assign("fcs.cluster", fcs.cluster, envir = .GlobalEnv)
    
    
    #fcs.cluster <- c("B cells", "Naive CD4 T cells", "Activated CD4 T cells", 
    #                 "Naive CD8 T cells", "Activated CD8 T cells", "Macrophages",
    #                 "Monocytes", "Neutrophils", "NK cells", 
    #                 "cDCs", "pDCs", "Eosinophils")
    
    
    fcs.cluster.label <- sprintf( "Cluster-%s", fcs.cluster )
    names( fcs.cluster.label ) <- fcs.cluster
    
    assign("fcs.cluster.label", fcs.cluster.label, envir = .GlobalEnv)
    
    
    # Tip: this controls the grouping of the fcs files that are generated based on the flowSOM clustering.
    # If some of the clusters are related (e.g., CD14+ and CD16+ monocytes), 
    # you might wish to group these into a single output folder.
    # If so, you may use a command such as "a" = c(1, 5, 7),
    # where the numbers are the numbered flowSOM clusters.
    fcs.cluster.group <- as.list(1:fcs.cluster.n)
    names(fcs.cluster.group) <- make.unique(rep(letters, length.out = fcs.cluster.n), sep='')
    
    assign("fcs.cluster.group", fcs.cluster.group, envir = .GlobalEnv)
    
    fcs.cluster.color <- rep( 
      fcs.color.pool, 
      ceiling( fcs.cluster.n / fcs.color.pool.n ) 
    )[ 1 : fcs.cluster.n ]
    names( fcs.cluster.color ) <- fcs.cluster
    
    assign("fcs.cluster.color", fcs.cluster.color, envir = .GlobalEnv)
    
    fcs.cluster.line.type <- rep( 
      fcs.line.type.pool, 
      ceiling( fcs.cluster.n / fcs.line.type.pool.n ) 
    )[ 1 : fcs.cluster.n ]
    names( fcs.cluster.line.type ) <- fcs.cluster
    
    assign("fcs.cluster.line.type", fcs.cluster.line.type, envir = .GlobalEnv)
    
    
    fcs.cluster.table.dir <- "./table_cluster"
    fcs.cluster.data.dir <- "./data_cluster"
    
    assign("fcs.cluster.table.dir", fcs.cluster.table.dir, envir = .GlobalEnv)
    assign("fcs.cluster.data.dir", fcs.cluster.data.dir, envir = .GlobalEnv)
    
    
    fcs.cluster.table.counts <- "cluster_counts"
    
    assign("fcs.cluster.table.counts", fcs.cluster.table.counts, envir = .GlobalEnv)
    
    
    fcs.cluster.data <- "cluster_group"
    
    assign("fcs.cluster.data", fcs.cluster.data, envir = .GlobalEnv)
    
    
    
    # dimensionality reduction parameters
    fcs.dmrd.gradient.color <- c( "black", "blue", "green", "yellow", "red" )
    fcs.dmrd.gradient.palette.n <- 100
    fcs.dmrd.density.palette <- colorRampPalette( fcs.dmrd.gradient.color )( 
      fcs.dmrd.gradient.palette.n )
    
    assign("fcs.dmrd.gradient.color", fcs.dmrd.gradient.color, envir = .GlobalEnv)
    assign("fcs.dmrd.gradient.palette.n", fcs.dmrd.gradient.palette.n, envir = .GlobalEnv)
    assign("fcs.dmrd.density.palette", fcs.dmrd.density.palette, envir = .GlobalEnv)
    
    
    fcs.dmrd.color.alpha <- 0.3
    
    assign("fcs.dmrd.color.alpha", fcs.dmrd.color.alpha, envir = .GlobalEnv)
    
    
    fcs.dmrd.group.title.size <- 8
    
    assign("fcs.dmrd.group.title.size", fcs.dmrd.group.title.size, envir = .GlobalEnv)
    
    
    fcs.dmrd.legend.title.size <- 7
    fcs.dmrd.legend.label.size <- 7
    fcs.dmrd.legend.point.size <- 3
    
    assign("fcs.dmrd.legend.title.size", fcs.dmrd.legend.title.size, envir = .GlobalEnv)
    assign("fcs.dmrd.legend.label.size", fcs.dmrd.legend.label.size, envir = .GlobalEnv)
    assign("fcs.dmrd.legend.point.size", fcs.dmrd.legend.point.size, envir = .GlobalEnv)
    
    
    fcs.dmrd.label.factor.width <- 0.1
    
    assign("fcs.dmrd.label.factor.width", fcs.dmrd.label.factor.width, envir = .GlobalEnv)
    
    
    # Tip: set the number of rows in the output figures here.
    fcs.dmrd.figure.nrow <- 2
    fcs.dmrd.figure.ncol <- ceiling( fcs.condition.n / fcs.dmrd.figure.nrow )
    
    assign("fcs.dmrd.figure.nrow", fcs.dmrd.figure.nrow, envir = .GlobalEnv)
    assign("fcs.dmrd.figure.ncol", fcs.dmrd.figure.ncol, envir = .GlobalEnv)
    
    
    fcs.dmrd.figure.width <- 3
    fcs.dmrd.figure.height <- 3
    
    assign("fcs.dmrd.figure.width", fcs.dmrd.figure.width, envir = .GlobalEnv)
    assign("fcs.dmrd.figure.height", fcs.dmrd.figure.height, envir = .GlobalEnv)
    
    
    
    # umap parameters
    fcs.umap.figure.lims.factor <- 0.8
    fcs.umap.figure.point.size <- 1.2
    
    assign("fcs.umap.figure.lims.factor", fcs.umap.figure.lims.factor, envir = .GlobalEnv)
    assign("fcs.umap.figure.point.size", fcs.umap.figure.point.size, envir = .GlobalEnv)
    
    
    fcs.umap.figure.dir <- "./figure_umap"
    
    fcs.umap.figure.plot <- "umap_plot"
    
    fcs.umap.cache.file.path <- paste0(folder_name, "/umap_cache.dat")
    
    assign("fcs.umap.figure.dir", fcs.umap.figure.dir, envir = .GlobalEnv)
    assign("fcs.umap.figure.plot", fcs.umap.figure.plot, envir = .GlobalEnv)
    assign("fcs.umap.cache.file.path", fcs.umap.cache.file.path, envir = .GlobalEnv)
    
    
    
    # cross-entropy test parameters
    # Tips: Set this depending on your RAM and number of groups.
    # You won't be able to analyze more than about 100k cells unless you have >32GB RAM.
    # The crossentropy test works best with at least 10k cells per group.
    # Multiple hypothesis testing will greatly reduce your ability to distinguish statistical differences.
    fcs.ce.diff.prob.sample.n <- NULL
    
    assign("fcs.ce.diff.prob.sample.n", fcs.ce.diff.prob.sample.n, envir = .GlobalEnv)
    
    
    # Tip: set to "ks" unless you have a good statistical reason for using rank testing.
    # In that case, use "rank" and "median".
    fcs.ce.diff.base.test <- "ks"
    fcs.ce.diff.base.dist <- "ks"
    
    assign("fcs.ce.diff.base.test", fcs.ce.diff.base.test, envir = .GlobalEnv)
    assign("fcs.ce.diff.base.dist", fcs.ce.diff.base.dist, envir = .GlobalEnv)
    
    
    fcs.ce.diff.test.alpha <- 0.05
    
    assign("fcs.ce.diff.test.alpha", fcs.ce.diff.test.alpha, envir = .GlobalEnv)
    
    
    fcs.ce.diff.figure.font.size <- 2
    fcs.ce.diff.figure.line.width <- 3
    
    assign("fcs.ce.diff.figure.font.size", fcs.ce.diff.figure.font.size, envir = .GlobalEnv)
    assign("fcs.ce.diff.figure.line.width", fcs.ce.diff.figure.line.width, envir = .GlobalEnv)
    
    
    fcs.ce.diff.figure.cdf.resolution <- 500
    fcs.ce.diff.figure.cdf.all.color <- "black"
    fcs.ce.diff.figure.cdf.all.label <- "All"
    
    assign("fcs.ce.diff.figure.cdf.resolution", fcs.ce.diff.figure.cdf.resolution, envir = .GlobalEnv)
    assign("fcs.ce.diff.figure.cdf.all.color", fcs.ce.diff.figure.cdf.all.color, envir = .GlobalEnv)
    assign("fcs.ce.diff.figure.cdf.all.label", fcs.ce.diff.figure.cdf.all.label, envir = .GlobalEnv)
    
    
    fcs.ce.diff.figure.dendrogram.weight.condition <- 1 : fcs.condition.n
    names( fcs.ce.diff.figure.dendrogram.weight.condition ) <- fcs.condition
    
    assign("fcs.ce.diff.figure.dendrogram.weight.condition", fcs.ce.diff.figure.dendrogram.weight.condition, envir = .GlobalEnv)
    
    
    fcs.ce.diff.figure.dendrogram.weight.cluster <- 1 : fcs.cluster.n
    names( fcs.ce.diff.figure.dendrogram.weight.cluster ) <- fcs.cluster
    
    assign("fcs.ce.diff.figure.dendrogram.weight.cluster", fcs.ce.diff.figure.dendrogram.weight.cluster, envir = .GlobalEnv)
    
    
    fcs.ce.diff.figure.cdf.width <- 1200
    fcs.ce.diff.figure.cdf.height <- 800
    
    assign("fcs.ce.diff.figure.cdf.width", fcs.ce.diff.figure.cdf.width, envir = .GlobalEnv)
    assign("fcs.ce.diff.figure.cdf.height", fcs.ce.diff.figure.cdf.height, envir = .GlobalEnv)
    
    
    fcs.ce.diff.figure.dendrogram.width <- 2000
    fcs.ce.diff.figure.dendrogram.height <- 800
    
    assign("fcs.ce.diff.figure.dendrogram.width", fcs.ce.diff.figure.dendrogram.width, envir = .GlobalEnv)
    assign("fcs.ce.diff.figure.dendrogram.height", fcs.ce.diff.figure.dendrogram.height, envir = .GlobalEnv)
    
    
    
    # cross-entropy test parameters for umap
    
    fcs.ce.diff.umap.figure.dir <- folder_name 
    
    fcs.ce.diff.umap.figure.cdf <- "umap_ce_diff_cdf"
    fcs.ce.diff.umap.figure.dendrogram <- "umap_ce_diff_dendrogram"
    fcs.ce.diff.umap.result <- "umap_ce_diff_result"
    
    fcs.ce.diff.umap.cache.file.path <- paste0(folder_name, "/umap_ce_diff_cache.dat")
    
    assign("fcs.ce.diff.umap.figure.dir", fcs.ce.diff.umap.figure.dir, envir = .GlobalEnv)
    assign("fcs.ce.diff.umap.figure.cdf", fcs.ce.diff.umap.figure.cdf, envir = .GlobalEnv)
    assign("fcs.ce.diff.umap.figure.dendrogram", fcs.ce.diff.umap.figure.dendrogram, envir = .GlobalEnv)
    assign("fcs.ce.diff.umap.result", fcs.ce.diff.umap.result, envir = .GlobalEnv)
    assign("fcs.ce.diff.umap.cache.file.path", fcs.ce.diff.umap.cache.file.path, envir = .GlobalEnv)
    
    
    # analyze_flow_cytometry.r ----
    # function to set random seed depending on base number and string
    
    set.seed.here <- function( seed.base, seed.char )
    {
      seed.add <- strtoi( substr( digest( seed.char, "xxhash32" ), 2, 8 ), 16 )
      seed.new <- seed.base + seed.add
      set.seed( seed.new )
      invisible( seed.new )
    }
    
    
    
    # check consistency in parameters
    
    stopifnot( names( fcs.channel.label ) == fcs.channel )
    
    stopifnot( names( fcs.condition.label ) == fcs.condition )
    stopifnot( names( fcs.condition.color ) == fcs.condition )
    stopifnot( names( fcs.condition.line.type ) == fcs.condition )
    
    stopifnot( names( fcs.cluster.label ) == fcs.cluster )
    stopifnot( names( fcs.cluster.color ) == fcs.cluster )
    stopifnot( names( fcs.cluster.line.type ) == fcs.cluster )
    
    stopifnot( sum( is.null( fcs.dmrd.data.sample.n ),
                    is.null( fcs.dmrd.data.sample.n.per.condition ),
                    is.null( fcs.dmrd.data.sample.n.per.sample ) ) >= 2 )
    
    stopifnot( unlist( fcs.cluster.group ) >= 1 &
                 unlist( fcs.cluster.group ) <= fcs.cluster.n )
    
    
    # create dirs
    
    figure.dir <- c(
      fcs.ce.diff.umap.figure.dir 
    )
    
    table.dir <- fcs.cluster.table.dir
    
    data.dir <- sapply( names( fcs.cluster.group ), function( fcg.name )
      sprintf( "%s/%s_%s", fcs.cluster.data.dir, fcs.cluster.data, fcg.name ) )
    
    for ( the.dir in c( figure.dir
    ) )
      if ( ! file.exists( the.dir ) )
        dir.create( the.dir, recursive = TRUE , showWarnings = FALSE)
    
    
    # read fcs data
    
    flow.data.filename.all <- list.files( fcs.data.dir, "\\.fcs$" )
    
    flow.data.filename <- grep( paste0( fcs.condition, collapse = "|" ),
                                flow.data.filename.all, value = TRUE )
    
    sample.name.format <- paste0( "%s.%0", fcs.sample.number.width, "d" )
    event.name.format <- paste0( "%s.%0", fcs.event.number.width, "d" )
    
    flow.data.filename.sample <- rep( "", length( flow.data.filename ) )
    names( flow.data.filename.sample ) <- flow.data.filename
    
    sample.idx.next <- rep( 1, fcs.condition.n )
    names( sample.idx.next ) <- fcs.condition
    
    flow.data <- lapply( flow.data.filename, function( flow.data.fn ) {
      #   cat( flow.data.fn, "\n" )
      sample.flow.frame <- read.FCS( file.path( fcs.data.dir, flow.data.fn ),
                                     transformation = NULL, truncate_max_range = FALSE )
      
      condition <- fcs.condition[ sapply( fcs.condition, grepl, flow.data.fn ) ]
      stopifnot( length( condition ) == 1 )
      
      sample.data <- exprs( sample.flow.frame )
      
      if ( ! all( fcs.channel %in% colnames( sample.data ) ) )
      {
        cat( sprintf( "File: %s\n", flow.data.fn ) )
        print( sort( fcs.channel[
          ! fcs.channel %in% colnames( sample.data ) ] ) )
        print( sort( colnames( sample.data ) ) )
        stop( "mismatch in names of fcs channels" )
      }
      
      sample.name <- sprintf( sample.name.format, condition,
                              sample.idx.next[ condition ] )
      
      sample.data <- sample.data[ , fcs.channel, drop = FALSE ]
      
      event.n <- nrow( sample.data )
      if ( event.n > 0 ) {
        event.name <- sprintf( event.name.format, sample.name, 1 : event.n )
        rownames( sample.data ) <- event.name
      }
      
      flow.data.filename.sample[ flow.data.fn ] <<- sample.name
      sample.idx.next[ condition ] <<- sample.idx.next[ condition ] + 1
      
      sample.data
    } )
    
    flow.data <- do.call( rbind, flow.data )
    
    # define samples
    flow.sample <- flow.data.filename.sample
    names( flow.sample ) <- NULL
    
    stopifnot( flow.sample ==
                 unique( sub( "\\.[0-9]+$", "", rownames( flow.data ) ) ) )
    
    flow.sample.n <- length ( flow.sample )
    
    flow.sample.condition <- factor( sub( "\\.[0-9]+$", "", flow.sample ),
                                     levels = fcs.condition )
    names( flow.sample.condition ) <- flow.sample
    
    # reorder samples to follow order of conditions
    flow.sample <- flow.sample[ order( flow.sample.condition ) ]
    flow.sample.condition <- flow.sample.condition[ flow.sample ]
    
    flow.sample.label <- sapply( flow.sample, function( fs ) {
      sample.cond <- sub( "^(.*)\\.[0-9]+$", "\\1", fs )
      sample.num <- sub( "^.*\\.([0-9]+)$", "\\1", fs )
      sprintf( "%s-%s", fcs.condition.label[ sample.cond ], sample.num )
    } )
    
    flow.sample.filename <- sapply( flow.sample, function( fs  )
      names( which( flow.data.filename.sample == fs ) ) )
    
    # define events
    flow.event <- rownames( flow.data )
    flow.event.n <- length( flow.event )
    
    flow.event.sample <- factor( sub( "\\.[0-9]+$", "", flow.event ),
                                 levels = flow.sample )
    names( flow.event.sample ) <- flow.event
    
    flow.event.condition <- factor( sub( "\\.[0-9]+$", "", flow.event.sample ),
                                    levels = fcs.condition )
    names( flow.event.condition ) <- flow.event
    
    # reorder events to follow order of samples
    flow.event.order <- order( flow.event.sample )
    
    flow.data <- flow.data[ flow.event.order, ]
    flow.event <- flow.event[ flow.event.order ]
    flow.event.sample <- flow.event.sample[ flow.event.order ]
    flow.event.condition <- flow.event.condition[ flow.event.order ]
    
    flow.event.sample.n <- as.vector( table( flow.event.sample ) )
    names( flow.event.sample.n ) <- flow.sample
    
    flow.event.condition.n <- as.vector( table( flow.event.condition ) )
    names( flow.event.condition.n ) <- fcs.condition
    
    # flow.data.filename
    
    # table( flow.sample.condition )
    
    # str( flow.data )
    # flow.event.condition.n
    # flow.event.sample.n
    
    
    # define figure parameters for samples
    
    flow.sample.color <- fcs.condition.color[ flow.sample.condition ]
    names( flow.sample.color ) <- flow.sample
    
    flow.sample.color.single <- unlist( lapply( fcs.condition, function( fc ) {
      cond.sample.n <- sum( flow.sample.condition == fc )
      rep(
        fcs.color.pool,
        ceiling( cond.sample.n / fcs.color.pool.n )
      )[ 1 : cond.sample.n ]
    } ) )
    names( flow.sample.color.single ) <- flow.sample
    
    flow.sample.line.type <- fcs.condition.line.type[ flow.sample.condition ]
    names( flow.sample.line.type ) <- flow.sample
    
    flow.sample.line.type.single <- unlist( lapply( fcs.condition, function( fc ) {
      cond.sample.n <- sum( flow.sample.condition == fc )
      rep(
        fcs.line.type.pool,
        ceiling( cond.sample.n / fcs.line.type.pool.n )
      )[ 1 : cond.sample.n ]
    } ) )
    names( flow.sample.line.type.single ) <- flow.sample
    
    flow.ce.diff.figure.dendrogram.weight.sample <-
      fcs.ce.diff.figure.dendrogram.weight.condition[ flow.sample.condition ]
    names( flow.ce.diff.figure.dendrogram.weight.sample ) <- flow.sample
    
    
    # select data for dimensionality reduction
    
    set.seed.here( fcs.seed.base, "select data for dimensionality reduction" )
    
    {
      if ( ! is.null( fcs.dmrd.data.sample.n ) )
      {
        if ( fcs.dmrd.data.sample.n < flow.event.n )
          dmrd.data.idx <- sort( sample( flow.event.n, fcs.dmrd.data.sample.n ) )
        else
          dmrd.data.idx <- 1 : flow.event.n
      }
      else if ( ! is.null( fcs.dmrd.data.sample.n.per.condition ) )
      {
        dmrd.data.idx <- unlist( sapply( fcs.condition, function( fc ) {
          fc.idx <- which( flow.event.condition == fc )
          if ( fcs.dmrd.data.sample.n.per.condition < length( fc.idx ) )
            sort( sample( fc.idx, fcs.dmrd.data.sample.n.per.condition ) )
          else
            fc.idx
        } ) )
        names( dmrd.data.idx ) <- NULL
      }
      else if ( ! is.null( fcs.dmrd.data.sample.n.per.sample ) )
      {
        dmrd.data.idx <- unlist( sapply( flow.sample, function( fs ) {
          fs.idx <- which( flow.event.sample == fs )
          if ( fcs.dmrd.data.sample.n.per.sample < length( fs.idx ) )
            sort( sample( fs.idx, fcs.dmrd.data.sample.n.per.sample ) )
          else
            fs.idx
        } ) )
        names( dmrd.data.idx ) <- NULL
      }
      else
        dmrd.data.idx <- 1 : flow.event.n
    }
    
    dmrd.data <- flow.data[ dmrd.data.idx, ]
    
    dmrd.event.sample <- flow.event.sample[ dmrd.data.idx ]
    dmrd.event.condition <- flow.event.condition[ dmrd.data.idx ]
    
    # str( dmrd.data )
    # table( dmrd.event.condition )
    # table( dmrd.event.sample )
    
    # calculate umap representation
    
    {
      if ( fcs.use.cached.results && file.exists( fcs.umap.cache.file.path ) )
      {
        cat( "Using cached results for umap...\n" )
        
        load( fcs.umap.cache.file.path )
      }
      else
      {
        set.seed.here( fcs.seed.base, "calculate umap representation" )
        
        cat( "Calculating umap...\n" )
        
        umap.config <- umap.defaults
        umap.config$n_epochs <- fcs.umap.iter.n
        umap.config$verbose <- TRUE
        
        umap.result <- umap( dmrd.data, config = umap.config )
        
        save( umap.result, file = fcs.umap.cache.file.path )
      }
    }
    
    umap.data <- umap.result$layout
    dimnames( umap.data ) <- NULL
    
    # str( umap.result )
    
    # calculate cross-entropy test for umap by condition
    
    set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )
    
    ce.diff.test.umap.res <- ce.diff.test.umap( 
      umap.result$knn$distances, umap.result$knn$indexes, 
      umap.data, umap.result$config, 
      dmrd.event.condition, 
      partition.label = fcs.condition.label, 
      partition.color = fcs.condition.color, 
      partition.line.type = fcs.condition.line.type, 
      base.test = fcs.ce.diff.base.test, 
      base.dist = fcs.ce.diff.base.dist, 
      prob.sample.n = fcs.ce.diff.prob.sample.n, 
      dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.condition, 
      result = file.path( fcs.ce.diff.umap.figure.dir, 
                          sprintf( "%s_condition.txt", fcs.ce.diff.umap.result ) ), 
      cdf.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                              sprintf( "%s_condition.png", fcs.ce.diff.umap.figure.cdf ) ), 
      dendrogram.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                                     sprintf( "%s_condition.png", fcs.ce.diff.umap.figure.dendrogram ) ) )
    
    
    # # calculate cross-entropy test for umap by sample
    # 
    # set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )
    # 
    # ce.diff.test.umap.res <- ce.diff.test.umap( 
    #   umap.result$knn$distances, umap.result$knn$indexes, 
    #   umap.data, umap.result$config, 
    #   dmrd.event.sample, 
    #   partition.label = flow.sample.label, 
    #   partition.color = flow.sample.color, 
    #   partition.line.type = flow.sample.line.type.single, 
    #   base.test = fcs.ce.diff.base.test, 
    #   base.dist = fcs.ce.diff.base.dist, 
    #   prob.sample.n = fcs.ce.diff.prob.sample.n, 
    #   dendrogram.order.weight = flow.ce.diff.figure.dendrogram.weight.sample, 
    #   result = file.path( fcs.ce.diff.umap.figure.dir, 
    #                       sprintf( "%s_sample.txt", fcs.ce.diff.umap.result ) ), 
    #   cdf.figure = file.path( fcs.ce.diff.umap.figure.dir, 
    #                           sprintf( "%s_sample.png", fcs.ce.diff.umap.figure.cdf ) ), 
    #   dendrogram.figure = file.path( fcs.ce.diff.umap.figure.dir, 
    #                                  sprintf( "%s_sample.png", fcs.ce.diff.umap.figure.dendrogram ) ) )
    
  }
  
  visualize_CE <- function() {
    # List all files in the specified directory
    files_in_directory <- list.files(
      path = file.path(rfs),
      full.names = TRUE
    )
    
    # Check the number of files in the directory
    num_files <- length(files_in_directory)
    
    if (num_files > 2) {
      CE <- read.table(
        file.path(ce_folder, "umap_ce_diff_result_condition.txt"),
        header = FALSE,
        skip = 2
      )
      
      colnames(CE) <- c("File A", "V2", "File B", "V4", "V5", "D value", "V7", "V8", "P value", "V10", "V11", "Adjusted P value")
      
      CE <- CE %>% select(-"V2", -"V4", -"V5", -"V7", -"V8", -"V10", -"V11")
      
      CE$significant <- ifelse(CE$`Adjusted P value` < as.numeric(threshold_CE), "YES", "NO")
      
      write.csv(CE, file = file.path(ce_folder, "CE_result.csv"), row.names = FALSE)
      View(CE)
      yes_count <- sum(CE$significant == "YES")
      cat(paste0("Amount of significant combinations: ", yes_count, "\n"))
      
    } else if (num_files >= 1) {
      cat("There are 1-2 files or conditions being compared Not running the code.\n")
      
    } else {
      cat("No files found in the directory. Not running the code.\n")
    }
  }
  
  visualize_CE_if_needed <- function() {
    unique_batch_values <- unique(md$batch)
    
    if (length(unique_batch_values) > 2) {
      cat("\nMore than 3 unique batch values detected. An overview will be automatically generated.\n")
      visualize_CE()
    } else if (length(unique_batch_values) == 2) {
      cat("\nOnly 1-2 files or conditions are being compared. Manual verification is required.\n")
    } else {
      cat("Not enough unique batch values for comparison. The visualize_CE function will not run.\n")
    }
  }
  
  # source functions copied from https://github.com/AdrianListon/Cross-Entropy-test/tree/main/flow%20analysis 
  ce.diff.test <- function( cross.entropy, event.partition, partition.label, partition.color, partition.line.type, base.test, base.dist, dendrogram.order.weight, result, cdf.figure, dendrogram.figure ) {
    partition <- levels( event.partition )
    partition.n <- length( partition )
    
    cross.entropy.split <- split( cross.entropy, event.partition )
    
    if ( base.test == "ks" )
    {
      if ( partition.n == 2 )
      {
        ks.test.res <- ks.test( cross.entropy.split[[ 1 ]], 
                                cross.entropy.split[[ 2 ]] )
        
        test.res <- list( ks.single = ks.test.res )
      }
      else if ( partition.n > 2 )
      {
        comparison.n <- partition.n * ( partition.n - 1 ) / 2
        
        ks.pair <- vector( "list", comparison.n )
        comparison <- character( comparison.n )
        D.stat <- numeric( comparison.n )
        p.value <- numeric( comparison.n )
        
        k <- 1
        
        for ( i in 1 : ( partition.n - 1 ) )
          for ( j in (i+1) : partition.n )
          {
            ks.pair[[ k ]] <- ks.test( cross.entropy.split[[ i ]], 
                                       cross.entropy.split[[ j ]] )
            
            comparison[ k ] <- sprintf( "%s - %s", 
                                        partition.label[ partition[ i ] ], 
                                        partition.label[ partition[ j ] ] )
            
            D.stat[ k ] <- ks.pair[[ k ]]$statistic
            
            p.value[ k ] <- ks.pair[[ k ]]$p.value
            
            k <- k + 1
          }
        
        p.value.adj <- p.adjust( p.value, "holm" )
        
        test.res <- list( ks.multiple = list( ks.pair = ks.pair, 
                                              comparison = comparison, D.stat = D.stat, 
                                              p.value = p.value, p.value.adj = p.value.adj ) )
      }
      else
        stop( "no partitions for testing cross-entropy differences" )
    }
    else if ( base.test == "rank" )
    {
      if ( partition.n == 2 )
      {
        wilcox.test.res <- wilcox.test( cross.entropy.split[[ 1 ]], 
                                        cross.entropy.split[[ 2 ]] )
        
        test.res <- list( rank.single = wilcox.test.res )
      }
      else if ( partition.n > 2 )
      {
        kruskal.test.res <- kruskal.test( cross.entropy, 
                                          event.partition )
        
        dunn.test.res <- dunn.test( cross.entropy, event.partition, 
                                    method = "holm", alpha = fcs.ce.diff.test.alpha, altp = TRUE, 
                                    kw = FALSE, table = FALSE, list = TRUE )
        
        test.res <- list( rank.multiple = list( 
          kruskal = kruskal.test.res, dunn = dunn.test.res ) )
      }
      else
        stop( "no partitions for testing cross-entropy differences" )
    }
    else
      stop( "wrong base test for cross-entropy differences" )
    
    if ( ! is.null( result ) )
    {
      result.file <- file( result, "w" )
      sink( result.file )
      
      tr.name <- names( test.res )
      stopifnot( length( tr.name ) == 1 )
      
      if ( tr.name == "ks.single" )
      {
        cat( "\n**  Kolmogorov-Smirnov test\n")
        
        print( test.res$ks.single )
      }
      else if ( tr.name == "ks.multiple" )
      {
        cat( "\n**  Multiple Kolmogorov-Smirnov tests with Holm correction\n\n")
        
        comparison.width <- max( nchar( test.res$ks.multiple$comparison ) )
        
        for ( i in 1 : length( test.res$ks.multiple$comparison ) )
          cat( sprintf( "%-*s\t\tD = %g\t\tpv = %g\t\tadj-pv = %g\n", 
                        comparison.width, 
                        test.res$ks.multiple$comparison[ i ], 
                        test.res$ks.multiple$D.stat[ i ], 
                        test.res$ks.multiple$p.value[ i ], 
                        test.res$ks.multiple$p.value.adj[ i ] ) )
      }
      else if ( tr.name == "rank.single" )
      {
        cat( "\n**  Wilcoxon rank sum test\n")
        
        print( test.res$rank.single )
      }
      else if ( tr.name == "rank.multiple" )
      {
        cat( "\n**  Kruskal-Wallis rank sum test\n")
        
        print( test.res$rank.multiple$kruskal )
        
        cat( "\n**  Dunn post-hoc test with Holm correction\n\n")
        
        comparison.width <- max( nchar( 
          test.res$rank.multiple$dunn$comparison ) )
        
        for ( i in 1 : length( test.res$rank.multiple$dunn$comparisons ) )
          cat( sprintf( "%-*s\t\tZ = %g\t\tpv = %g\t\tadj-pv = %g\n", 
                        comparison.width, 
                        test.res$rank.multiple$dunn$comparisons[ i ], 
                        test.res$rank.multiple$dunn$Z[ i ], 
                        test.res$rank.multiple$dunn$altP[ i ], 
                        test.res$rank.multiple$dunn$altP.adjusted[ i ] ) )
      }
      else
      {
        sink()
        close( result.file )
        stop( "unknown test in ce-diff result" )
      }
      
      sink()
      close( result.file )
    }
    
    if ( ! is.null( cdf.figure ) )
    {
      if ( is.null( partition.label ) )
        partition.label = partition
      
      if ( is.null( partition.color ) )
        partition.color <- rainbow( partition.n )
      
      if ( is.null( partition.line.type ) )
        partition.line.type <- rep( 1, partition.n )
      
      png( filename = cdf.figure, width = fcs.ce.diff.figure.cdf.width, 
           height = fcs.ce.diff.figure.cdf.height )
      
      par( mar = c( 5.5, 6, 2, 1.5 ) )
      
      plot( ecdf( cross.entropy ), ylim = c( 0, 1 ), 
            xlab = "Cross-entropy", ylab = "CDF", main = "", 
            cex.lab = 3, cex.axis = 2.5, 
            col = fcs.ce.diff.figure.cdf.all.color, 
            lwd = fcs.ce.diff.figure.line.width - 1, do.points = FALSE )
      
      for ( pall in partition )
      {
        ces <- cross.entropy.split[[ pall ]]
        ces.n <- length( ces )
        
        if ( ces.n < fcs.ce.diff.figure.cdf.resolution ) {
          plot( ecdf( ces ), col = partition.color[ pall ], 
                lty = partition.line.type[ pall ], 
                lwd = fcs.ce.diff.figure.line.width, 
                do.points = FALSE, add = TRUE )
        }
        else {
          ecdf.x <- sort( ces )
          ecdf.y <- 1 : ces.n / ces.n
          
          lines( ecdf.x, ecdf.y, col = partition.color[ pall ], 
                 lty = partition.line.type[ pall ], 
                 lwd = fcs.ce.diff.figure.line.width )
        }
      }
      
      legend( "bottomright", 
              legend = c( fcs.ce.diff.figure.cdf.all.label, partition.label ), 
              col = c( fcs.ce.diff.figure.cdf.all.color, partition.color ), 
              lty = partition.line.type, 
              lwd = fcs.ce.diff.figure.line.width, 
              cex = fcs.ce.diff.figure.font.size )
      
      dev.off()
    }
    
    if ( ! is.null( dendrogram.figure ) && partition.n > 2 )
    {
      cross.entropy.dist <- matrix( 0, nrow = partition.n, 
                                    ncol = partition.n )
      
      for ( i in 1 : ( partition.n - 1 ) )
        for ( j in (i+1) : partition.n )
        {
          if ( base.dist == "ks" )
            cross.entropy.dist[ i, j ] <- ks.test( 
              cross.entropy.split[[ i ]], 
              cross.entropy.split[[ j ]] 
            )$statistic
          else if ( base.dist == "median" )
            cross.entropy.dist[ i, j ] <- abs( 
              median( cross.entropy.split[[ i ]] ) - 
                median( cross.entropy.split[[ j ]] )
            )
          else
            stop( "wrong base dist for cross-entropy differences" )
          
          cross.entropy.dist[ j, i ] <- cross.entropy.dist[ i, j ]
        }
      
      cross.entropy.hclust <- hclust( as.dist( cross.entropy.dist ) )
      
      if ( ! is.null( dendrogram.order.weight ) )
        cross.entropy.hclust <- as.hclust( reorder( 
          as.dendrogram( cross.entropy.hclust ), 
          dendrogram.order.weight, 
          agglo.FUN = mean
        ) )
      
      if ( is.null( partition.label ) )
        partition.label = partition
      
      png( filename = dendrogram.figure, 
           width = fcs.ce.diff.figure.dendrogram.width, 
           height = fcs.ce.diff.figure.dendrogram.height )
      
      par( mar = c( 5, 5.6, 4, 1.4 ) )
      
      plot( cross.entropy.hclust, 
            labels = partition.label, hang = -1, 
            xlab = "", ylab = "", main = "", sub = "", cex.axis = 3, 
            cex = fcs.ce.diff.figure.font.size )
      
      dev.off()
    }
    
    test.res
  }
  ce.diff.test.umap <- function( orig.dist, orig.knn, umap.data, umap.param, event.partition, partition.label = NULL, partition.color = NULL, partition.line.type = NULL, base.test = "ks", base.dist = "ks", prob.sample.n = NULL, dendrogram.order.weight = NULL, result = NULL, cdf.figure = NULL, dendrogram.figure = NULL ) {
    stopifnot( nrow( orig.dist ) == nrow( umap.data ) && 
                 nrow( orig.dist ) == length( event.partition ) )
    
    data.n <- nrow( orig.dist )
    
    if ( ! is.null( prob.sample.n ) && prob.sample.n < data.n )
      prob.sample.idx <- sample( data.n, prob.sample.n )
    else
      prob.sample.idx <- 1 : data.n
    
    if ( fcs.use.cached.results && 
         file.exists( fcs.ce.diff.umap.cache.file.path ) )
    {
      cat( "Using cached results for probability\n" )
      
      load( fcs.ce.diff.umap.cache.file.path )
    }
    else
    {
      # cat( "Calculating probability...\n" )
      
      orig.umap.prob <- calculate.probability.umap( orig.dist, orig.knn, 
                                                    umap.data, umap.param )
      
      save( orig.umap.prob, file = fcs.ce.diff.umap.cache.file.path )
    }
    
    cross.entropy.all <- calculate.fuzzy.cross.entropy( 
      orig.umap.prob$orig[ prob.sample.idx, ], 
      orig.umap.prob$umap[ prob.sample.idx, ] 
    )
    
    event.partition.all <- event.partition[ prob.sample.idx ]
    
    ce.diff.test( 
      cross.entropy.all, 
      event.partition.all, 
      partition.label, partition.color, partition.line.type, 
      base.test, base.dist, 
      dendrogram.order.weight, 
      result, cdf.figure, dendrogram.figure
    )
  }
  calculate.probability.umap <- function( orig.dis, orig.kn, umap.dat, umap.param ) {
    orig.dat.n <- nrow( orig.dis )
    umap.dat.n <- nrow( umap.dat )
    
    stopifnot( orig.dat.n == umap.dat.n )
    
    # get nearest neighbors in original space and their distances
    
    orig.self.idx <- sapply( 1 : orig.dat.n, function( ri ) {
      ri.idx <- which( orig.kn[ ri, ] == ri )
      ifelse( length( ri.idx ) == 1, ri.idx, NA )
    } )
    
    stopifnot( ! is.na( orig.self.idx ) )
    
    orig.neigh <- t( sapply( 1 : orig.dat.n, function( ri ) 
      orig.kn[ ri, - orig.self.idx[ ri ] ] ) )
    
    orig.dis.reduc <- t( sapply( 1 : orig.dat.n, function( ri ) 
      orig.dis[ ri, - orig.self.idx[ ri ] ] ) )
    
    # calculate probabilities associated to distances in original space
    
    umap.sigma <- apply( orig.dis.reduc, 1, function( dd ) {
      umap.sigma.error <- function( ss, dd ) {
        p <- exp( - pmax( 0, dd - min( dd ) ) / ss )
        sum( p ) - log2( length( p ) )
      }
      
      dd.ascen <- sort( dd )
      dd.ascen <- dd.ascen[ dd.ascen > 0 ]
      ss.lower <- dd.ascen[ 1 ]
      
      dd.descen <- sort( dd, decreasing = TRUE )
      dd.descen <- dd.descen[ ! is.infinite( dd.descen ) ]
      ss.upper <- dd.descen[ 1 ]
      
      while( umap.sigma.error( ss.upper, dd ) < 0 )
      {
        ss.lower <- ss.upper
        ss.upper <- 2 * ss.upper
      }
      
      while( umap.sigma.error( ss.lower, dd ) > 0 )
      {
        ss.upper <- ss.lower
        ss.lower <- ss.lower / 2
      }
      
      uniroot( umap.sigma.error, dd, 
               interval = c( ss.lower, ss.upper ), 
               tol = ( ss.upper - ss.lower ) * .Machine$double.eps^0.25 )$root
    } )
    
    orig.prob <- t( sapply( 1 : orig.dat.n, function( i ) 
      exp( - pmax( 0, orig.dis.reduc[ i, ] - min( orig.dis.reduc[ i, ] ) ) / 
             umap.sigma[ i ] )
    ) )
    
    # symmetrize probabilities in original space
    
    for ( i in 1 : orig.dat.n )
      for ( j2 in 1 : length( orig.neigh[ i, ] ) )
      {
        j <- orig.neigh[ i, j2 ]
        
        i2 <- match( i, orig.neigh[ j, ] )
        
        if ( ! is.na( i2 ) )
        {
          if ( j > i )
          {
            sym.prob <- orig.prob[ i, j2 ] + orig.prob[ j, i2 ] - 
              orig.prob[ i, j2 ] * orig.prob[ j, i2 ]
            orig.prob[ i, j2 ] <- sym.prob
            orig.prob[ j, i2 ] <- sym.prob
          }
        }
      }
    
    # get distances in umap space for closest neighbors in original space
    
    umap.dist2 <- t( sapply( 1 : umap.dat.n, function( i )
      sapply( orig.neigh[ i, ], function( j )
        sum( ( umap.dat[ i, ] - umap.dat[ j, ] )^2 )
      )
    ) )
    
    # calculate probabilities associated to distances in umap representation
    
    umap.a <- umap.param$a
    umap.b <- umap.param$b
    
    umap.prob <- t( apply( umap.dist2, 1, function( dd2 ) 
      p <- 1 / ( 1 + umap.a * dd2 ^ umap.b )
    ) )
    
    list( orig = orig.prob, umap = umap.prob )
  }
  calculate.fuzzy.cross.entropy <- function( prim.prob, secd.prob ) {
    prim.prob.n <- nrow( prim.prob )
    secd.prob.n <- nrow( secd.prob )
    
    prim.prob.m <- ncol( prim.prob )
    secd.prob.m <- ncol( secd.prob )
    
    stopifnot( prim.prob.n == secd.prob.n && prim.prob.m == secd.prob.m )
    
    sapply( 1 : prim.prob.n, function( i ) 
      - sum( prim.prob[ i, ] * log( secd.prob[ i, ] ) + 
               ( 1 - prim.prob[ i, ] ) * log( 1 - secd.prob[ i, ] ) )
    )
  }
  
  ## Chunk 8: check CSV values for CN. ----
  
  testCV_new <- function(fsom, cluster_values = testCVcluster, plot = TRUE, verbose = FALSE, seed = 1) {
    
    nClus <- fsom$map$nNodes
    cluster_labels <- FlowSOM::GetClusters(fsom)
    
    # Determine metacluster labels
    meta_cl <- list()
    for(mc in cluster_values){
      if(verbose) message("Computing ", mc, " metaclusters")
      meta_cl[[as.character(mc)]] <-
        FlowSOM::metaClustering_consensus(fsom$map$codes,
                                          mc,
                                          seed = seed)
    }
    meta_cl[[as.character(nClus)]] <- seq_len(nClus)
    
    # Percentages assigned to each of the clusters per file
    pctgs <- list()
    for(mc in as.character(c(cluster_values, nClus))){
      counts <- matrix(0,
                       nrow = length(unique(fsom$data[,"File"])),
                       ncol = as.numeric(mc),
                       dimnames = list(unique(fsom$data[,"File"]),
                                       as.character(seq_len(as.numeric(mc)))))
      tmp <- table(fsom$data[,"File"],
                   meta_cl[[mc]][cluster_labels])
      counts[rownames(tmp), colnames(tmp)] <- tmp
      pctgs[[mc]] <- t(apply(counts, 1,
                             function(x){ 100 * x/sum(x) }))
    }
    
    # Coefficient of variation for each of the percentages
    cvs <- list()
    for(mc in as.character(c(cluster_values, nClus))){
      cvs[[mc]] <- apply(pctgs[[mc]],
                         2,
                         function(x){ stats::sd(x) / mean(x)})
    }
    res <- named.list(pctgs, cvs, meta_cl)
    if(plot){
      
      # save results in a pdf
      pdf(file = file.path(rsl, "02 CytoNorm CVS values", "CVS_summary_plots.pdf"))
      PlotOverviewCV(fsom, res)
      dev.off()
      
      PlotOverviewCV(fsom, res)
      
      
    }
    
    return(res)
  }
  
  checkVarCN = function(){
    
    if(performCN & file.exists(file.path(wd, cellType, "01 FCS files", "03 normalized files"))){
      continue <- readline("You already performed normalization, do you want to overwrite (and therefore remove) these files? (y/n): ")
      
      if (tolower(continue) == "n") {
        stop("Please put performCN to FALSE!")
      }
    }
    
    if(cvsCheck & file.exists(file.path(rsl, "02 CytoNorm CVS values"))){
      continue <- readline("You already performed a cvs check, do you want to overwrite (and therefore remove) these files? (y/n): ")
      
      if (tolower(continue) == "n") {
        stop("Please put cvsCheck to FALSE!")
      }
    }
    
    
    if(!(nClustersCN %in% testCVcluster)){
      stop("nClustersCN has to be a value that is within the range of testCV cluster, please adjust accordingly.")
    }
    
  }  
  
  checkCVSValues = function(verbose = TRUE){
    
    if(to_import == "c"){
      assign("dirCN", file.path(wd, cellType, "01 FCS files", "02 cleaned and transformed files", "Preprocessed", "QC", "PeacoQC_results", "fcs_files"), envir = .GlobalEnv)
      message("Performing CytoNorm on cleaned files")
    }
    if(to_import == "b"){
      assign("dirCN", trns, envir = .GlobalEnv)
      message("Performing CytoNorm on transformed FCS files")
    } 
    if(to_import == "a"){
      assign("dirCN", fcs, envir = .GlobalEnv)
      message("Performing CytoNorm on files provided in FCS files folder")
    }

    to_import = "d"
    assign("to_import", to_import, envir = .GlobalEnv)
    
    #don't change
    if(verbose) cat("Checking variables\n")
    assign("filesCN", list.files(dirCN), envir = .GlobalEnv)
    if(unique(filesCN == md$file_name) > 1){
      stop("File order of filesCN does not match file order of md")
    }      
    if(verbose) cat("Building dataframe\n")

    dataCN = data.frame(File = filesCN,
                        Path = file.path(dirCN, filesCN),
                        Type = md$sample_control,
                        Batch = md$batch,
                        stringsAsFactors = FALSE)
    
    dataCN$Type = c("control" = "Train", "sample" = "Validation")[dataCN$Type]
    
    assign("dataCN", dataCN, envir = .GlobalEnv)
    assign("train_data", dplyr::filter(dataCN, Type == "Train"), envir = .GlobalEnv)
    assign("validation_data", dplyr::filter(dataCN, Type == "Validation"), envir = .GlobalEnv)
    assign("channelsCN", dplyr::filter(panel)$fcs_colname, envir = .GlobalEnv)
    
    #check parameters
    if(verbose) cat("Preparing FlowSOM\n\n")
    fsom <- prepareFlowSOM(train_data$Path,
                           channelsCN, 
                           nCells = nCellsCN, #total number of cells to use for clustering. Default 1.000.000
                           FlowSOM.params = list(xdim = xdimCN, #default 10
                                                 ydim = ydimCN, #default 10
                                                 nClus = nClustersCN, #default 10
                                                 scale = FALSE), #default FALSE
                           seed = 1, #don't change, for reproducability
                           transformList = NULL,
                           verbose = verbose) #TRUE if you want progress updates 
    
    assign("fsom", fsom, envir = .GlobalEnv)
    
    #don't change
    if(cvsCheck == TRUE) {
      save_dir <- file.path(rsl, "02 CytoNorm CVS values")
      
      if (!file.exists(save_dir)) {
        dir.create(save_dir, showWarnings = FALSE)
      }
      
      cvs <- testCV_new(fsom, cluster_values = testCVcluster)
      cvs$pctgs$`10`
      names(cvs$cvs)
      
      # Initialize max_cvs_values as an empty list
      max_cvs_values <- list()
      
      for(i in names(cvs$cvs)) {
        # Generate the plot
        plot_data <- cvs$cvs[[i]]
        plot_title <- paste("CVS Plot -", i)
        plot(plot_data, main = plot_title)
        
        # Save the plot as a PNG file in the specified directory
        png_filename <- file.path(save_dir, paste0("cvs_plot_", i, ".png"))
        png(png_filename)
        plot(plot_data, main = plot_title)
        dev.off()
        
        # Find and save the maximum value from the current plot
        max_value <- max(cvs$cvs[[i]])
        max_cvs_values[[i]] <- max_value
        
      }
      
      assign("max_cvs_values", max_cvs_values, envir = .GlobalEnv)
      
      # Extracting the names and maximum values from the max_cvs_values list
      plot_names <- names(max_cvs_values)
      max_values <- unlist(max_cvs_values)
      
      # Creating a line graph to display the maximum values
      plot(1:length(max_values), max_values, type = "l", 
           xlab = "Clusters tested", ylab = "Maximum Value", 
           xaxt = "n", ylim = c(0, max(max_values) + 1), 
           main = "Maximum Values from CVS Plots")
      axis(1, at = 1:length(max_values), labels = plot_names)
      
      # Save the line graph as a PNG file in a specific location
      png_filename <- file.path(save_dir, "max_values_plot.png")
      png(png_filename)
      plot(1:length(max_values), max_values, type = "l", 
           xlab = "CVS Plots", ylab = "Maximum Value", 
           xaxt = "n", ylim = c(0, max(max_values) + 1), 
           main = "Maximum Values from CVS Plots")
      axis(1, at = 1:length(max_values), labels = plot_names)
      dev.off()
      
      assign("cvs", cvs, envir = .GlobalEnv)
    } else {
      message("No CVS check performed.")
    }
    
  }
  
  ## Chunk 9: run CytoNorm. ----
  
  performCytoNorm = function(verbose = TRUE){
    
    if(performCN){
      
      if(verbose) message("Starting CytoNorm:")
      if(verbose) cat("> Checking variables\n")
      
      assign("filesCN", list.files(dirCN), envir = .GlobalEnv)
      if(unique(filesCN == md$file_name) > 1){stop("File order of filesCN does not match file order of md")}
      
      if(verbose) cat("> Building dataframe\n")
      #create dataframe with information for cytonorm
      dataCN = data.frame(File = filesCN,
                          Path = file.path(dirCN, filesCN),
                          Type = md$sample_control,
                          Batch = md$batch,
                          stringsAsFactors = FALSE)
      
      dataCN$Type = c("control" = "Train", "sample" = "Validation")[dataCN$Type]
      
      assign("dataCN", dataCN, envir = .GlobalEnv)
      assign("train_data", dplyr::filter(dataCN, Type == "Train"), envir = .GlobalEnv)
      assign("validation_data", dplyr::filter(dataCN, Type == "Validation"), envir = .GlobalEnv)
      assign("channelsCN", dplyr::filter(panel)$fcs_colname, envir = .GlobalEnv)      
      if(verbose) cat("> Preparing FlowSOM\n\n")
      fsom <- prepareFlowSOM(train_data$Path,
                             channelsCN, 
                             nCells = nCellsCN, #total number of cells to use for clustering. Default 1.000.000
                             FlowSOM.params = list(xdim = xdimCN, #default 10
                                                   ydim = ydimCN, #default 10
                                                   nClus = nClustersCN, #default 10
                                                   scale = FALSE), #default FALSE
                             seed = 1, #don't change, for reproducability
                             transformList = NULL,
                             verbose = TRUE) #TRUE if you want progress updates 
      
      assign("fsom", fsom, envir = .GlobalEnv)
      
      train_data$Batch = 1:length(list.files(rfs))
      
      
      if(verbose) message("\nBuilding model")
      modelCN = CytoNorm.train(files = train_data$Path,
                               labels = train_data$Batch,
                               channels = channelsCN,
                               transformList = NULL,
                               FlowSOM.params = list(nCells = nCellsCN, #total number of cells to use for clustering. Default 1000000
                                                     xdim = xdimCN, #default 10
                                                     ydim = ydimCN, #default 10
                                                     nClus = nClustersCN, #default 10
                                                     scale = FALSE), #default FALSE
                               
                               normMethod.train = QuantileNorm.train, #normalization method to use for each cluster. Default = QuantileNorm.train
                               
                               normParams = list(nQ = nQuantilesCN, #number of quantiles
                                                 goal = goalDistributionCN, #goal distribution. Default = mean, can also be nQ numeric values or one of the batch labels
                                                 verbose = TRUE, #if TRUE progress updates are printed
                                                 limit = normParamLimit,
                                                 plot = FALSE, #if TRUE, plot is generated showing all quantiles
                                                 plotTitle = "Quantiles"), #title to use in the plot
                               
                               seed = 1, #don't change, for reproducability
                               outputDir = rsl,
                               verbose = TRUE) #TRUE if you want progress updates
      
      assign("modelCN" , modelCN, envir = .GlobalEnv)
      
      
      if(verbose) cat("\n> Applying model to data\n")
      
      CytoNorm.normalize(model = modelCN,
                         files = dataCN$Path,
                         labels = dataCN$Batch,
                         transformList = NULL,
                         transformList.reverse = NULL,
                         normMethod.normalize = QuantileNorm.normalize,
                         outputDir = paste0(wd, "/", cellType ,"/01 FCS files/03 normalized files"),
                         prefix = "Norm_",
                         clean = TRUE, #if TRUE in between files are removed
                         verbose = FALSE) #if TRUE progress updates are printed
      
      total.list = list.files(paste0(wd, "/", cellType ,"/01 FCS files/03 normalized files"))
      to.keep = which(grepl("Norm", list.files(paste0(wd, "/", cellType ,"/01 FCS files/03 normalized files"))))
      to.remove = total.list[-to.keep]
      setwd(paste0(wd, "/", cellType ,"/01 FCS files/03 normalized files"))
      file.remove(to.remove)
      
      file_name2 = c(list.files(paste0(wd, "/", cellType ,"/01 FCS files/03 normalized files")))
      fsNorm = read.flowSet(path=paste0(wd, "/", cellType ,"/01 FCS files/03 normalized files"),transformation=FALSE, truncate_max_range=FALSE)
      
      assign("file_name2", file_name2, envir = .GlobalEnv)
      assign("fsNorm", fsNorm, envir = .GlobalEnv)
      
      if(verbose) cat("> Preparing metadata file\n")
      
      md_2 = read.xlsx(mtd, 1, stringsAsFactors=FALSE)
      md_2$file_name <- paste("Norm_", md_2$file_name, sep = "")

      assign("md_2", md_2, envir = .GlobalEnv)
      
      # save CN metadata as an excel file
      write.xlsx(md_2, file = file.path(wd, cellType, paste(expID, "Catalyst", cellType, "metadata CN.xlsx")), row.names = FALSE)
      
      if(verbose) cat("> Creating single cell experiment\n")
      
      sce_CytoNorm = prepData(fsNorm, #flowset
                              panel = panel, #panel file
                              md = md_2, #metadata file
                              features=markers_CF, #vector of channel names
                              FACS = TRUE, 
                              transform = FALSE,
                              md_cols = list(file = "file_name2", id = "sample_id", factors = factors_list))
      
      assayNames(sce_CytoNorm)[1] = "exprs" #this line seems to refuse to run within the function
      
      assign("sce_CytoNorm", sce_CytoNorm, envir = .GlobalEnv)
      assign("sce", sce_CytoNorm, envir = .GlobalEnv)
      assign("fs", fsNorm, envir = .GlobalEnv)
      
      
      if(verbose) cat("\nDone with performing normalization with CytoNorm.")
      
      # rename the settings file 
      file.rename(from = file.path(rsl, "CytoNorm_FlowSOM.RDS"),
                  to = file.path(rsl, "/03 CytoNorm diagnostic plots/00 Settings_CytoNorm_FlowSOM.RDS"))
      

      
    } else {
      message("No CytoNorm performed, no values or variables changed.")
    }

  if(to_import == "d"){
    # import new md file
    md <- readxl::read_excel(file.path(wd, cellType, paste(expID, "Catalyst", cellType, "metadata CN.xlsx")))
    md <- dplyr::rename(md, file_name = file_name2)
    md_2 <- dplyr::rename(md_2, file_name = file_name2)
    
    assign("md", md, envir = .GlobalEnv)
    
    cat("\n\nMetadata file of cytonorm is now imported.")  }
  }
  
  writeDiagnosticsCN = function(){
    
    if (write_files == TRUE){
      
      if(verbose) cat("\nWriting diagnostic plots.\n\n")
      
      setwd(rsl)
      dir.create("03 CytoNorm diagnostic plots", showWarnings = FALSE)
      assign("res_CN", paste0(rsl, "/03 CytoNorm diagnostic plots"), envir = .GlobalEnv)
      
      if(verbose) cat("Plots sucessfully created:")
      
      # Create the folder/directory
      if (!file.exists(res_CN)) {dir.create(res_CN, showWarnings = FALSE)}
      
      #writes pdf file containing bar graph of cell counts
      pdf(file = paste0(res_CN, "/Counts ", cellType, ".pdf"), width = 30, height = 10)
      print(plotCounts(sce, color_by = "sample_id") +
              ggtitle("Count summary of everything in total."))
      
      if(length(conList)>1){
        for (i in 1:length(conList)){
          abc = conList[i]
          print(plotCounts(filterSCE(sce, condition == abc), color_by = "sample_id")+
                  ggtitle(paste0("Plot for: ", abc)))
        }
      }else{
        print(plotCounts(sce, color_by = "sample_id"))
      }  
      dev.off()
      
      if(verbose) cat("\n> Counts plot.\n")
      
      if(TCs == TRUE){
        
        #writes pdf file containing MDS plot of technical controls
        if(length(which(md$sample_control == "control"))>2){
          pdf(file = paste0(res_CN, "/MDS ", cellType, " TCs.pdf"), width = 10, height = 10)
          print(pbMDS(filterSCE(sce, sample_control == "control"), by = c("sample_id"), fun = "median", features = NULL, assay = "exprs", label_by = NULL, color_by = "sample_id"))
          dev.off()
          if(verbose) cat("> MDS TC plot.\n")
        }
        
        # writes expression plot of technical controls  
        pdf(file = paste0(res_CN, "/Expression ", cellType, " TCs.pdf"), width = 30, height = 10)
        print(plotExprs(filterSCE(sce, sample_control == "control"), color_by = "sample_id"))
        dev.off()
        if(verbose) cat("> Expression TC plot.\n")
      }
      
      #writes pdf file containing histogram of transformed expression
      pdf(file = paste0(res_CN, "/Expression ", cellType, ".pdf"), width = 30, height = 10)
      for (i in 1:length(conList)){
        abc = conList[i]
        print(plotExprs(filterSCE(sce, condition == abc), color_by = "sample_id")+
                ggtitle(paste0("Plot for: ", abc)))
      }
      dev.off()
      if(verbose) cat("> Expression plot.\n\n")
    }
    
    assign("cellNoTot", sum(n_cells(sce)), envir = .GlobalEnv)
    
    if(TCs == TRUE){
      assign("sce_c", (filterSCE(sce, sample_control == "control")), envir = .GlobalEnv)
      assign("sce", filterSCE(sce, sample_control != "control"), envir = .GlobalEnv)
      cat("\n\nControl group successfully removed from data\n")
      
    }
    
    assign("cellNoNC", sum(n_cells(sce)), envir = .GlobalEnv)
    
    if (write_files == TRUE){
      
      #writes pdf file containing heatmap of relative marker expression
      pdf(file = paste0(res_CN, "/Expression_heatmap ", cellType, ".pdf"), width = 10, height = 10)
      print(plotExprHeatmap(sce, bin_anno = FALSE, row_anno = exprHMRowAnno, row_clust=TRUE, by=c("sample_id"), fun = c("median")))
      dev.off()
      if(verbose) cat("\n> Expression heatmap containing relative marker expressions.\n")
      
      #writes pdf file containing non-redundancy scores, allowing for useful marker selection for downstream analysis
      pdf(file = paste0(res_CN, "/Non-redundancy_score ", cellType, ".pdf"), width = 10, height = 10)
      print(plotNRS(sce, features = NULL, color_by = NRSColorBy, assay = "exprs"))
      dev.off()
      if(verbose) cat("> Heatmap and NRS plots.\n")
      
      #writes pdf file containing non-redundancy scores, allowing for useful marker selection for downstream analysis
      pdf(file = paste0(res_CN, "/MDS ", cellType, ".pdf"), width = 10, height = 10)
      print(pbMDS(sce, by = "sample_id", fun = "median", features = NULL, assay = "exprs", label_by = NULL, color_by = MDSColorBy))
      dev.off()
      if(verbose) cat("> MDS plot.\n")
    }
    
  }
  
  

# II. FlowSOM clustering ----

  ## Chunk 2: loading data into the working environment, importing md and panel files and creating sce object. ---- 
  
  load_data = function(){
      
      specifyChannels()
      createSCE()
      
      for (i in seq_along(fs)){keyword(fs[[i]])[["$CYT"]] <- "FACS"}
      assign("cellNoNC", sum(n_cells(sce)), envir = .GlobalEnv)
      assign("cellNoTot", sum(n_cells(sce)), envir = .GlobalEnv)
    
  }
  
  ## Chunk 3: down sampling, clustering and dimension reduction.----
  
  downSample = function(){
    
    assign("conList", c(as.character(unique(md$condition))), envir = .GlobalEnv)
    
    # Downsampling
    message("Initiating downsampling")
    message("\nBefore downsampling: \n")
    print(sce)

    # tabulate number of cells in each sample
    n_cells <- table(sce$sample_id)
    
    # exclude samples with with less than 'min_cells'
    assign("rmv", n_cells > min_cells, envir = .GlobalEnv)
    assign("sce", filterSCE(sce, sample_id %in% names(which(rmv))), envir = .GlobalEnv)
    
    # downsample samples with more than 'max_cells'
    cells_by_sample <- split(seq(ncol(sce)), sce$sample_id)
    cells_keep <- unlist(lapply(cells_by_sample, 
                                function(.) sample(.,min(length(.), max_cells))))
    assign("sce", sce[, cells_keep], envir = .GlobalEnv)
    
    message("\nAfter downsampling:")
    print(sce)

    #writes pdf file containing bar graph of cell counts
      if (length(conList) > 1) {
        pdf(file = paste0(rsl, "/Counts ", cellType, " post downsampling.pdf"), width = 30, height = 10)
        print(plotCounts(sce, color_by = "sample_id"))
        for (i in 1:length(conList)) {
          abc = conList[i]
          print(plotCounts(filterSCE(sce, condition == abc), color_by = "sample_id"))
        }
        dev.off()
        cat("Plot Counts successfully written\n")
      } else {
        pdf(file = paste0(rsl, "/Counts ", cellType, " post downsampling.pdf"), width = 30, height = 10)
        print(plotCounts(sce, color_by = "sample_id"))
        dev.off()
        cat("Plot Counts successfully written\n")
      }
    
    
    assign("exclCellNo", cellNoNC - sum(n_cells(sce)), envir = .GlobalEnv)
    assign("exclCellPer", paste0((exclCellNo/cellNoNC)*100, "%"), envir = .GlobalEnv)
    assign("exclFileNo", length(grep("FALSE", c(data.frame(rmv)[,1]))), envir = .GlobalEnv)

  }
  
  removeMarkers = function(){
    markers_CL = markers_CF
    for (i in excluded_markers){
      markers_CL = markers_CL[-grep(panel$fcs_colname[grep(i, panel$antigen)], markers_CL)]
      sce = sce[-grep(i, rownames(sce))]
    }
    assign("markers_CL", markers_CL, envir = .GlobalEnv)
    assign("sce", sce, envir = .GlobalEnv)
    assign("tested_antigens", c(rownames(sce)), envir = .GlobalEnv)
  }
  
  performClustering = function(){
    assign("meta", gsub(" ", "", paste("meta", metClust)), envir = .GlobalEnv)##creates variable equal to amount of metaclusters for k=meta
    sce <- cluster(sce, features = type_markers(sce), xdim = FSxdim, ydim=FSydim, maxK=metClust, seed=seedFS) ##xdim * ydim = total SOM clusters; maxK = total metaclusters

    if(length(table(cluster_ids(sce)))<(FSxdim*FSydim)){
      stop(paste0("Only ", length(table(cluster_ids(sce)))," of ", (FSxdim*FSydim), " SOMnodes are filled. To fix this first execute the createSCE() function in the console. Then increase the amount of cells in the FlowSOM run (min_cells and max_cells) and/or reduce the amount of SOMnodes (FSxdim and FSydim"))
    }
    
    assign("sce" , sce, envir = .GlobalEnv)
  }
  
  performDimRed = function(){
    set.seed(clusterSeed)
    if ("UMAP" %in% anType){
      sce <- runDR(sce, dr=c("UMAP"), cells = UMAP_cells, features=type_markers(sce))
      assign("sce" , sce, envir = .GlobalEnv)
    }
    if ("T-SNE" %in% anType){
      sce <- runDR(sce, dr=c("TSNE"), cells = UMAP_cells, features=type_markers(sce))
      assign("sce" , sce, envir = .GlobalEnv)
    }
  }
  
  
  
  ## Chunk 4: determination requirement of cluster merging ----
  
  plot_sequences_DR <- function() {
    
    if (!plot_sequences) {stop("continue to next chunk")}
    
    message("Printing FlowSOM expression heatmaps and Dimension Reduction Sequences for Expression Levels:")
    
    # Create main folders
    main_folders <- c("00 Expression levels per marker in DR", "01 Cluster sequences", "02 Expression heatmap per cluster")
    for(folder in main_folders){ dir.create(path  = file.path(rsl, "00 Cluster number determination", folder), recursive = T, showWarnings = F) }
    
    #print flowSOM expression heatmaps
    flowsom_folder = file.path(rsl, "00 Cluster number determination", "02 Expression heatmap per cluster")
    
    for (nr in start_nr:end_nr){
      
      pdf(file = paste0(flowsom_folder, "/Expression_Heatmap_", nr, ".pdf"), width = 12, height = 10, title = "Expression_Heatmap_merged")  
      print(plotExprHeatmap(sce, k = paste0("meta",nr), bin_anno = TRUE, row_anno = TRUE, row_clust=TRUE, col_clust = FALSE, bars = TRUE, perc = TRUE, by=c("cluster_id"), fun = c("median")))
      dev.off()
      
    }
    
    
    if ("UMAP" %in% anType) {
      umap_folder1 = file.path(rsl, "00 Cluster number determination", "00 Expression levels per marker in DR")
      
      # Print UMAPs for expression levels of all markers
      for (i in rownames(sce)) {
        png(file.path(umap_folder1, paste0("UMAP_", i, ".png")), width = 3600, height = 3000, res = 300)
        print(plotDR(sce, "UMAP", color_by = i, facet_by = NULL) + theme_classic())
        dev.off()
      }
      
      # Print UMAPs for clusters
      umap_folder2 = file.path(rsl, "00 Cluster number determination", "01 Cluster sequences")
      dir.create(umap_folder2, recursive = TRUE, showWarnings = FALSE)
      
      for (nr in start_nr:end_nr){
        png(file.path(umap_folder2, paste0("UMAP_MetaClusters_", nr, ".png")), width = 10, height = 10, units = "in", res = 300, title = paste0("UMAP of ", nr, " metaclusters"))
        print(plotDR_new(sce, "UMAP", color_by = paste0("meta", nr)))
        dev.off()
        
        pdf(file.path(flowsom_folder, paste0("Expression_Heatmap_", nr, ".pdf")), width = 12, height = 10, title = "Expression_Heatmap_merged")  
        print(plotExprHeatmap(sce, k = paste0("meta",nr), bin_anno = TRUE, row_anno = TRUE, row_clust=TRUE, col_clust = FALSE, bars = TRUE, perc = TRUE, by=c("cluster_id"), fun = c("median")))
        dev.off()
      }
    }
    
    
    if ("T-SNE" %in% anType) {
      
      tsne_folder_expr = file.path(rsl, "00 Cluster number determination", "00 Expression levels per marker in DR")
      
      # Print T-SNEs for expression levels of all markers
      for (i in rownames(sce)) {
        png(file.path(tsne_folder_expr, paste0("TSNE_", i, ".png")), width = 3600, height = 3000, res = 300)
        print(plotDR(sce, "TSNE", color_by = i, facet_by = NULL) + theme_classic())
        dev.off()
      }
      
      tsne_folder_clusters = file.path(rsl, "00 Cluster number determination", "01 Cluster sequences")
      
      # Print T-SNEs for clusters
      for (nr in start_nr:end_nr){
        png(file.path(tsne_folder_clusters, paste0("TSNE_MetaClusters_", nr, ".png")), width = 10, height = 10, units = "in", res = 300, title = paste0("TSNE of ", nr, " metaclusters"))
        print(plotDR_new(sce, "TSNE", color_by = paste0("meta", nr)))
        dev.off()
      }
    }
    
    
    cat("> Done printing DR sequences.\n")
  }
  
  ## Chunk 5: cluster merging. ----
  
  createCMfolder <- function() {
    if (performClusterMerging) {
      
      if(to_import == "b"){ imported_data = "trns"}
      if(to_import == "c"){ imported_data = "pqcd"}
      if(to_import == "d"){ imported_data = "cn"}
      
      folder_name <- paste0('01 Cluster merging of ', imported_data) 
      cmr <- file.path(rsl, folder_name)
      
      if (!dir.exists(cmr)) {
        dir.create(cmr)
      } 
      
      assign("cmr", cmr, envir = .GlobalEnv)
    } else {
      cat("No cluster merging will be performed")
    }
  }
  
  clusterMergingFile = function(metClust, writeCM = FALSE) {
    if (!performClusterMerging) return
    
    cmf <- paste0(cmr, "/merging_table1.xlsx")
    
    if (writeCM) {
      data <- data.frame(
        original_cluster = 1:metClust,
        new_cluster = rep(NA, metClust)
      )
      
      write.xlsx(data, cmf, sheetName = "merging_table1", row.names = FALSE)
      cat("Excel file for merging clusters successfully created.\n")
      stop("Adapt your merging file! And put writeCM to FALSE!")
    } else if (file.exists(cmf)) {
      file_data <- read.xlsx(cmf, sheetName = "merging_table1")
      
      if (anyNA(file_data)) {
        stop("The file merging_table1.xlsx already exists, but it contains NA values. You should update the file with non-NA values or set writeCM to TRUE to overwrite it.\n")
      }
    } else {
      cat("No file named merging_table1.xlsx found...\n")
      cat("> You should either set writeCM to TRUE to create the file or upload your own Excel sheet matching the requirements.\n")
    }
    
    if (!writeCM) {
      merging_table1 <- read_excel(cmf, col_names = TRUE)
      merging_table1$new_cluster <- as.character(merging_table1$new_cluster)
      merging_table1$original_cluster <- as.character(merging_table1$original_cluster)
      
      assign("merging_table1", merging_table1, envir = .GlobalEnv)
    }
    
    assign("cmf", cmf, envir = .GlobalEnv)
  }
  
  plot_MC_before = function() {
    
    if (plotDiagnosticsBefore == TRUE) {
      
      pdf(file = file.path(cmr, "00_FlowSOM_heatmap_cluster_before.pdf"), width = 10, height = 10)  
      
      print(plotExprHeatmap(sce, 
                            k = paste0("meta",customMeta), 
                            bin_anno = TRUE, 
                            row_anno = TRUE, 
                            row_clust= TRUE, 
                            col_clust = FALSE, 
                            bars = TRUE, 
                            perc = TRUE, 
                            by=c("cluster_id"), 
                            fun = c("median")))
      dev.off()
      
      if ("UMAP" %in% anType) {
        png(file = file.path(cmr, "00_UMAP_cluster_before.png"), width = 10, height = 10, units = "in", res = 300)
        print(plotDR_new(sce, "UMAP", color_by = paste0("meta", customMeta)) +
                theme_classic())
        dev.off()
      }
      
      if ("T-SNE" %in% anType) {
        png(file = file.path(cmr, "00_TSNE_cluster_before.png"), width = 10, height = 10, units = "in", res = 300)
        print(plotDR_new(sce, "TSNE", color_by = paste0("meta", customMeta)) +
                theme_classic())
        dev.off()
      }
      
      message("FlowSOM expression heatmap and specified plots are saved.")
      
    } 
    
  }
  
  plot_MC_results = function(){
    
    if (plotDiagnosticsCM == TRUE) {
      
      pdf(file = file.path(cmr, "01_FlowSOM_heatmap_cluster_merged.pdf"), width = 10, height = 10)  
      
      print(plotExprHeatmap(sce_merge, 
                            k = "merging1", 
                            bin_anno = TRUE, 
                            row_anno = TRUE, 
                            row_clust= TRUE, 
                            col_clust = FALSE, 
                            bars = TRUE, 
                            perc = TRUE, 
                            by=c("cluster_id"), 
                            fun = c("median")))
      dev.off()
      
      if ("UMAP" %in% anType) {
        png(file = file.path(cmr, "01_UMAP_cluster_merged.png"), width = 10, height = 10, units = "in", res = 300)
        print(plotDR_new(sce_merge, "UMAP", color_by = "merging1") +
                theme_classic())
        dev.off()
      }
      
      if ("T-SNE" %in% anType) {
        png(file = file.path(cmr, "01_TSNE_cluster_merged.png"), width = 10, height = 10, units = "in", res = 300)
        print(plotDR_new(sce_merge, "TSNE", color_by = "merging1") +
                theme_classic())
        dev.off()
      }
      
      message("FlowSOM expression heatmap and specified plots are saved.")
      
    } 
  }
  
  save_MC_freq = function(){
    
    
    if ( ! file.exists( file.path(cmr) )) {dir.create( file.path(cmr), recursive = TRUE )}
    
    
    # Define the file path
    file_path <- file.path(cmr, "freq_som_merging1.xlsx")
    
    res <- table(cluster_id = cluster_ids(sce_merge, "merging1"), 
                 sample_id = sample_ids(sce_merge))
    
    fq <- prop.table(res, 2) * 100 
    
    # Save the data to Excel file
    write.xlsx(as.data.frame(fq), file = file_path)
    cat("Results saved successfully.\n\n")
    
  }  
  
  cluster_merging = function(){
    
    if (performClusterMerging){ 
      
      # create folder for cluster merging results in results folder
      createCMfolder()
      # create excel sheet for cluster merging
      clusterMergingFile(customMeta, writeCM)
      # plot clusters before merging
      plot_MC_before()
      # apply manual merging
      sce_merge <- mergeClusters(sce, k = paste0("meta", customMeta), table = merging_table1, id = "merging1") 
      assign("sce_merge", sce_merge, envir = .GlobalEnv)
      # plot UMAP and flowsom heatmap after cluster merging
      plot_MC_results() 
      # save frequencies
      save_MC_freq()
      
      message("cluster merging performed")
      
    } else{message("No cluster merging performed.")}
    
  }

  createSettingsFile <- function(prefix) {
    
    if(to_import == "b"){ imported_data = "trns"}
    if(to_import == "c"){ imported_data = "pqcd"}
    if(to_import == "d"){ imported_data = "cn"}
    
    # Create the file path
    filePath <- file.path(rsl, paste0(prefix, expID, "_", cellType, ".txt"))
    
    # Open the file for writing
    fileConn <- file(filePath)
    
    # Define the settings information
    settingsInfo <- c(
      paste0("Files loaded into directory: ", imported_data),
      "",
      paste0("Total number of files: ", nrow(md)),
      paste0("Total number of cells: ", cellNoTot),
      "",
      "DOWNSAMPLING",
      paste0("Markers excluded from analysis: ", as.character(excluded_markers)),
      paste0("Cell threshold: ", as.character(min_cells)),
      paste0("Amount of cells excluded: ", as.character(exclCellNo)),
      paste0("Percentage of cells excluded: ", exclCellPer),
      paste0("Amount of files excluded: ", as.character(exclFileNo)),
      "",
      "FLOWSOM CLUSTERING",
      paste0("Amount of files in FlowSOM: ", as.character((nrow(md) - exclFileNo))),
      paste0("Total number of cells in FlowSOM: ", sum(n_cells(sce))),
      paste0("Max number of cells per file in FlowSOM: ", as.character(max_cells)),
      paste0("FlowSOM seed: ", as.character(seedFS)),
      paste0("FlowSOM metaclusters: ", as.character(metClust)),
      paste0("Amount of clusters chosen to continue: ", customMeta),
      paste0("Cluster merging performed: ", as.character(performClusterMerging)),
      "",
      paste0("Types of analysis done: ", as.character(anType)),
      paste0("UMAP number of cells per file: ", as.character(UMAP_cells)),
      paste0("T-SNE number of cells per file: ", as.character(TSNE_cells)),
      paste0("UMAP/T-SNE seed: ", as.character(clusterSeed)),
      ""
    )
    
    # Write the settings information to the file
    writeLines(settingsInfo, fileConn)
    
    # Close the file connection
    close(fileConn)
    
  }
  
  save_variables_rmd2 = function() {
    # Create a list to store common variables
    variables_to_save <- list(
      sce = sce,
      to_import = to_import,
      excluded_markers = excluded_markers, 
      min_cells = min_cells, 
      max_cells = max_cells, 
      metClust = metClust, 
      seedFS = seedFS, 
      FSxdim = FSxdim, 
      FSydim = FSydim, 
      anType = anType, 
      clusterSeed = clusterSeed,
      UMAP_cells = UMAP_cells, 
      TSNE_cells = TSNE_cells, 
      customMeta = customMeta, 
      performClusterMerging = performClusterMerging,
      exclCellNo = exclCellNo
    )
    
    # Add additional variables if performClusterMerging is TRUE
    if (performClusterMerging) {
      variables_to_save$sce_merge <- sce_merge
      variables_to_save$exclCellPer <- exclCellPer
      variables_to_save$exclFileNo <- exclFileNo
    }
    
    # Serialize and save to a file
    saveRDS(variables_to_save, file = file.path(rsl, "saved_variables.rds"))
  }
  
  
  ## Chunk 6: cluster evaluation ----
  
  save_CM_plots <- function() {
    if (!save_plots) return  # Do nothing if save_plots is FALSE
    
    plot_folder <- file.path(rsl, "02 Frequency boxplots")
    
    if (!dir.exists(plot_folder)) dir.create(plot_folder, recursive = TRUE)
    
    file_path <- file.path(plot_folder, paste0(name_CM_result, ".pdf"))
    
    if (file.exists(file_path)) {
      overwrite <- tolower(readline(prompt = "The file already exists. Overwrite? (y/n): "))
      if (overwrite == "y") {
        ggsave(file = file_path, width = 10, height = 10)
        cat("File overwritten.\n")
      } else {
        cat("File not overwritten.\n")
      }
    } else {
      ggsave(file = file_path, width = 10, height = 10)
    }
  }
  
  
  
  
# III. Downstream analysis ----
  
  ## Chunk 1: loading previous working environment ----
  
  load_variables_rmd2 <- function() {
    # Define the file path where variables were saved
    saved_file_path <- file.path(wd, cellType, "03 Results",  "02 FlowSOM_clustering", "saved_variables.rds")
    
    # Check if the file exists
    if (file.exists(saved_file_path)) {
      # Deserialize and load variables from the file
      loaded_variables <- readRDS(saved_file_path)
      
      # Access the variables
      assign("to_import", loaded_variables$to_import, envir = .GlobalEnv)
      assign("excluded_markers", loaded_variables$excluded_markers, envir = .GlobalEnv)
      assign("min_cells", loaded_variables$min_cells, envir = .GlobalEnv)
      assign("max_cells", loaded_variables$max_cells, envir = .GlobalEnv)
      assign("metClust", loaded_variables$metClust, envir = .GlobalEnv)
      assign("seedFS", loaded_variables$seedFS, envir = .GlobalEnv)
      assign("FSxdim", loaded_variables$FSxdim, envir = .GlobalEnv)
      assign("FSydim", loaded_variables$FSydim, envir = .GlobalEnv)
      assign("anType", loaded_variables$anType, envir = .GlobalEnv)
      assign("clusterSeed", loaded_variables$clusterSeed, envir = .GlobalEnv)
      assign("UMAP_cells", loaded_variables$UMAP_cells, envir = .GlobalEnv)
      assign("TSNE_cells", loaded_variables$TSNE_cells, envir = .GlobalEnv)
      assign("customMeta", loaded_variables$customMeta, envir = .GlobalEnv)
      assign("sce", loaded_variables$sce, envir = .GlobalEnv)
      assign("sce_merge", loaded_variables$sce_merge, envir = .GlobalEnv)
      
      assign("performClusterMerging", loaded_variables$performClusterMerging, envir = .GlobalEnv)
      assign("exclCellNo", loaded_variables$exclCellNo, envir = .GlobalEnv)
      assign("exclCellPer", loaded_variables$exclCellPer, envir = .GlobalEnv)
      assign("exclFileNo", loaded_variables$exclFileNo, envir = .GlobalEnv)
      
      
    } else {
      # Handle the case where the file does not exist
      stop("Saved variables file not found.\n")  }
  }
  
  checkFilesPresent <- function() {
    
    fcs_dir <- file.path(wd, cellType, "01 FCS files")
    
    raw_files_dir <- file.path(fcs_dir, "00 raw files")
    raw_files <- dir(raw_files_dir)
    
    cat("The following fcs files were detected:\n")
    
    check_and_print_files(raw_files_dir, "raw files")
    check_and_print_files(file.path(fcs_dir, "01 transformed files"), "transformed files")
    check_and_print_files(file.path(fcs_dir, "02 cleaned and transformed files/Preprocessed/QC/PeacoQC_results/fcs_files"), "cleaned and transformed files")
    check_and_print_files(file.path(fcs_dir, "03 normalized files"), "normalized files")
    
    cat("\nFiles loaded in the current environment: ")
    
    if(to_import == "a" | to_import == "b"){
      cat("logicle transformed files")
    }
    
    if(to_import == "c"){
      cat("cleaned and logicle transformed files")}
    
    if(to_import == "d"){
      cat("normalized files")
      
      
    }
    if(performClusterMerging) {cat("\n\nManual cluster merging has been performed")}
  }
  
  
  ## Chunk 2: create downstream analysis folder. ----
  
  createSubsetList = function(){
    subsetList = c()
    
    for (i in 1:length(subsetCat)){
      values = as.character(unique(sce_trunc_subs[[as.character(subsetCat[i])]]))
      
      if (length(values) > 1){
        formattedValues = paste(values, collapse = ", ")
        subsetList = append(subsetList, paste0(subsetCat[i], " = ", formattedValues))
      } else {
        subsetList = append(subsetList, paste0(subsetCat[i], " = ", values))
      }
    }
    
    assign("subsetList", subsetList, envir = .GlobalEnv)
  }
  
  createIdentifier = function(){
    
      vc = as.character(Sys.time()) #to create a unique identifier for every set of files the system date and time is appended to the pdf file name
      substr(vc, 14, 14) = "h"
      substr(vc, 17, 17) = "m"
      substr(vc, 20, 20) = "s"
      vc = substring(vc, 1, nchar(vc) - 6)
      
      customSuffix = if (performClusterMerging) "sce_merge" else "sce"
      
      vc = paste0(vc, "s ", customSuffix, " ", customMeta, "MCs")

      assign("vc", vc, envir = .GlobalEnv)
      assign("customSuffix", customSuffix, envir = .GlobalEnv)
      
      
    }
  
  writeFiles = function(){
    
    if (write_files == TRUE){
      
      name_DA = "00 PDFs of the results"
      dir.create(file.path(rsl, vc, name_DA), recursive = T, showWarnings = F)
      
      
      if(performClusterMerging){ k_val = "merging1" 
      } else {k_val = paste0("meta", customMeta)}
      
      # Set file paths
      pdf_file_path <- paste0(rsl, "/", vc, "/", name_DA, "/")
      
      # Function to generate and save plots
      generate_and_save_plot <- function(file_name, plot_function, width, height, title) {
        pdf(file = paste0(pdf_file_path, file_name), width = width, height = height, title = title)
        print(plot_function)
        dev.off()
      }
      
      # Plot 1: Expression Heatmap
      generate_and_save_plot("Expression_Heatmap_merged.pdf", 
                             plotExprHeatmap(sce_trunc_subs, k = k_val, bin_anno = TRUE, 
                                             row_anno = TRUE, row_clust = TRUE, col_clust = FALSE, 
                                             bars = TRUE, perc = TRUE, by = c("cluster_id"), fun = c("median")),
                             12, 10, "Expression Heatmap Merged")
      
      # Plot 2: SOM Code Heatmap
      pMHcheck <- if ("state" %in% panel$marker_class) "state" else "abundances"
      generate_and_save_plot("SOMcode_heatmap_v1.pdf", 
                             plotMultiHeatmap(sce_trunc_subs, hm1 = "type", hm2 = pMHcheck, 
                                              k = k_val, row_clust = TRUE, perc = TRUE, bars = TRUE),
                             15, 10, "SOM Code Heatmap V1")
      
      # Plot 3: t-SNE and PCA of SOM Codes
      generate_and_save_plot("T-SNE_and_PCA_SOMcodes.pdf", 
                             plotCodes(sce_trunc_subs, k = k_val),
                             10, 10, "t-SNE and PCA of SOM Codes")
      
      # Plot 4: Another SOM Code Heatmap
      generate_and_save_plot("SOMcode_heatmap_v2.pdf", 
                             plotMultiHeatmap(sce_trunc_subs, hm1 = "type", hm2 = pMHcheck, 
                                              k = paste0("som", FSxdim * FSydim), m = k_val, row_clust = TRUE),
                             10, 20, "SOM Code Heatmap V2")
      
      # Plot 5: Metacluster Relative Abundance
      generate_and_save_plot("Metacluster_relative_abundance.pdf", 
                             plotAbundances_new(sce_trunc_subs, k = k_val, 
                                                group_by = groupVar, by = "cluster_id", shape = shapeVar),
                             15, 15, "Metacluster Relative Abundance")
      
      # Plot 6: Delta Area
      generate_and_save_plot("Delta_area.pdf", 
                             delta_area(sce_trunc_subs),
                             10, 10, "Delta Area")
      
      # Display success message
      message("FlowSOM plots successfully written\n")
      
      
      # write a settings file
      output_directory = file.path(rsl, vc)
      file_path <- file.path(output_directory, paste0("Settings ", expID, " ", cellType, " ", vc, ".txt"))
      
      fileConn <- file(file_path)
      writeLines(c(
        paste0("Total files: ", nrow(md) - exclFileNo),
        paste0("Total cells in FlowSOM: ", sum(n_cells(sce))),
        paste0("Max cells per file in FlowSOM: ", as.character(max_cells)),
        paste0("FlowSOM seed: ", as.character(seedFS)),
        paste0("FlowSOM metaclusters: ", as.character(metClust)),
        paste0("Subsetted for: ", subsetList),
        paste0("Grouping variable: ", as.character(groupVar)),
        paste0("Shape variable: ", as.character(shapeVar)),
        paste0("Analysis types: ", as.character(anType)),
        paste0("UMAP cells per file: ", as.character(UMAP_cells)),
        paste0("T-SNE cells per file: ", as.character(TSNE_cells)),
        paste0("UMAP/T-SNE seed: ", as.character(clusterSeed)),
        paste0("Facet variable for UMAP and T-SNE: ", as.character(facetVar)),
        paste0("Dataset: ", as.character(customSuffix)),
        paste0("customMeta: ", as.character(customMeta))
      ), fileConn)
      close(fileConn)
      
      cat("Settings successfully written\n")
      
      
      write.xlsx(flowsomdf, paste0(rsl, "/", vc, "/FlowSOM_data_", customMeta, "MCs.xlsx"))
      
      if("UMAP" %in% anType){
        png(file = paste0(rsl, "/", vc, "/UMAP_MetaClusters.png"), width = 10, height = 10, units = "in", res = 300, title = "UMAP of metaclusters")
        print(plotDR_new(sce_trunc_subs, "UMAP", color_by = meta, facet_by = facetVar))
        dev.off()
        
        png(file = paste0(rsl, "/", vc, "/UMAP_Markers.png"), width = 10, height = 10, units = "in", res = 300, title = "UMAP of markers")
        print(plotDR_new(sce_trunc_subs, "UMAP", color_by = rownames(sce), facet_by = NULL))
        dev.off()
        
        png(file = paste0(rsl, "/", vc, "/UMAP_Batch.png"), width = 10, height = 10, units = "in", res = 300, title = "UMAP of batch")
        print(plotDR_new(sce_trunc_subs, "UMAP", color_by = "batch", facet_by = NULL))
        dev.off()
        
        cat("UMAP plots successfully written\n")

        }
      
      if("T-SNE" %in% anType){
        pdf(file = paste0(rsl, "/", vc, "/", expID, " ", cellType, " T-SNE.pdf"), width = 10, height = 10, title = "T-SNE of metaclusters")
        print(plotDR_new(sce_trunc_subs, "TSNE", color_by = meta, facet_by = facetVar))
        dev.off()
        cat("T-SNE plots successfully written\n")
      }
    }
  }
  
  plotGraphs = function(){
    if (plot_graph == TRUE){
      
      print(plotExprHeatmap(sce_trunc_subs, k = meta, bin_anno = TRUE, row_anno = TRUE, row_clust=TRUE, col_clust = FALSE, bars = TRUE, perc = TRUE, by=c("cluster_id"), fun = c("median")))
      print(plotCodes(sce_trunc_subs, k=meta))
      print(plotClusterHeatmap_new(sce_trunc_subs, hm2 = NULL, k = paste0("som", FSxdim*FSydim), m = meta, cluster_anno = TRUE, draw_freqs = TRUE))
      print(plotAbundances_new(sce_trunc_subs, k = meta, group_by = groupVar, by = "cluster_id", shape = shapeVar))

      if("UMAP" %in% anType){
        print(plotDR_new(sce_trunc_subs, "UMAP", color_by = meta, facet_by = facetVar))
        print(plotDR_new(sce_trunc_subs, "UMAP", color_by = rownames(sce), facet_by = NULL))
        print(plotDR_new(sce_trunc_subs, "UMAP", color_by = "batch", facet_by = NULL))
      }
      
      if("T-SNE" %in% anType){
        print(plotDR_new(sce_trunc_subs, "TSNE", color_by = meta, facet_by = facetVar))
      }
    }
  }
  
  writeFCS = function(){
    
    if(write_FCS == TRUE){
      dir.create(file.path(rsl, vc, "01 FCS files"), recursive = T, showWarnings = F)
      sce_trunc_subs[[meta]] = cluster_ids(sce_trunc_subs, meta)
      redDim = reducedDim(sce_trunc_subs, "UMAP") #create list of reduced dim values
      sceExp = sce_trunc_subs[, !is.na(redDim[, 1])] #remove events containing no reduced dim values
      fsExport = sce2fcs(sceExp, split_by = "sample_id", assay = "exprs", keep_cd = TRUE, keep_dr = TRUE)
      write.flowSet(fsExport, outdir = paste0(rsl,"/", vc,"/01 FCS files"))
      file.remove(paste0(rsl,"/", vc, "/01 FCS files/annotation.txt"))
      concat = concatenate_fcs_files(files.list = paste0(rsl,"/", vc,"/01 FCS files/", list.files(paste0(rsl,"/", vc,"/01 FCS files"))))
      flowCore::write.FCS(concat, filename = paste0(rsl,"/", vc,"/01 FCS files/concatenated.fcs"))
      
    }
  }
  
  ## Chunk 3: statistical analysis. ----
  
  runStat = function(){
    if(runStatistics){
  setwd(paste0(rsl, "/", vc))
  dir.create("02 Statistics")
  
  #import and prepare dataframe containing metacluster abundances
  dfStat = read.xlsx(paste0(rsl, "/", vc, "/FlowSOM_data_", customMeta, "MCs.xlsx"), 1, stringsAsFactors = FALSE) #import data for statistics
  dfStat$cluster_id = factor(dfStat$cluster_id, levels = c(1:length(unique(dfStat$cluster_id)))) #change order of cluster_IDs for faceting
  dfStat[[groupVar]] = factor(dfStat[[groupVar]], levels = groupOrder) #change order of x-axis group
  patientIDlist = c()
  for(i in unique(dfStat$sample_id)){ #append patient_id's to data fame
    sampleIDindex = which(md$sample_id == i)
    patientIDlist = append(patientIDlist, rep(md$patient_id[sampleIDindex], customMeta))
  }
  dfStat$patient_id = patientIDlist
  
  statGroup = colnames(dfStat[ncol(dfStat)-1]) #get the x-axis group
  
  #if samples are paired, remove all donors that do not occur in all groups specified by groupVar
  if(pairedData == TRUE){
    patientIDlist = c(unique(dfStat$patient_id))
    for(i in unique(dfStat[[statGroup]])){
      tempPatientIDlist = dfStat[dfStat[[statGroup]] == i,]$patient_id
      dfStat = dfStat[which(dfStat$patient_id %in% tempPatientIDlist),]
    }
  }
  
  
  dfStatHold = dfStat
  plotList = list()
  aov_df = data.frame()
  pwc_df = data.frame()
  
  #running analysis of variance test
  for(i in unique(dfStat$cluster_id)){
    dfStat = subset(dfStatHold, cluster_id == i)
    if(AOV == "kruskal.test" & pairedData == FALSE){
      res.aov = dfStat %>% kruskal_test(as.formula(paste0("Freq ~ ",statGroup)))
    }else if(AOV != "anova" & pairedData == TRUE){
      res.aov = dfStat %>% friedman_test(as.formula(paste0("Freq ~ ",statGroup, "| patient_id")))
    }else if(AOV == "anova" & pairedData == FALSE){
      res.aov = dfStat %>% anova_test(as.formula(paste0("Freq ~ ",statGroup)))
    }else if(AOV == "anova" & pairedData == TRUE){
      res.aov = dfStat %>% anova_test(as.formula(paste0("Freq ~ ",statGroup, "+ Error(patient_id/", statGroup, ")")))
    }else{
      stop('Invalid enrty in AOV, or pairedData variable. Valid entries are "kruskal.test" and "anova". Set pairedData to TRUE or FALSE')
    }
    
    temp_aov_df = as.data.frame(res.aov)
    temp_aov_df$cluster_id = c(i)
    aov_df = rbind(aov_df, temp_aov_df)
    
    
    if(PWC == "wilcox.test" & pairedData == FALSE){
      pwc = dfStat %>% wilcox_test(as.formula(paste0("Freq ~ ", statGroup)), paired = FALSE, p.adjust.method = multCompAdjust)
    }else if(PWC == "wilcox.test" & pairedData == TRUE){
      pwc = dfStat %>% wilcox_test(as.formula(paste0("Freq ~ ", statGroup)), paired = TRUE, p.adjust.method = multCompAdjust)
    }else if(PWC == "t.test" & pairedData == FALSE){
      pwc = dfStat %>% t_test(as.formula(paste0("Freq ~ ", statGroup)), paired = FALSE, p.adjust.method = multCompAdjust)
    }else if(PWC == "t.test" & pairedData == TRUE){
      pwc = dfStat %>% t_test(as.formula(paste0("Freq ~ ", statGroup)), paired = TRUE, p.adjust.method = multCompAdjust)
    }else{
      stop('Invalid entry in PWC, or pairedData variable. Valid entries for PWC are "wilcox.test" and "t.test". Set pairedData to TRUE or FALSE')
    }
    
    temp_pwc_df = as.data.frame(pwc)
    temp_pwc_df$cluster_id = c(rep(i, nrow(temp_pwc_df)))
    pwc_df = rbind(pwc_df, temp_pwc_df)
    
    pwc <- pwc %>% add_xy_position(x = statGroup)
    
    
    statPlot =  ggplot(dfStat, aes_string(y="Freq", x=statGroup))+
      geom_boxplot(aes_string(color = statGroup, fill = statGroup),
                   position = position_dodge(), alpha = 0.2, 
                   outlier.color = NA, show.legend = FALSE) +
      geom_point(aes_string(color = statGroup),
                 position = position_jitter(width = 0.2)) +
      facet_wrap("cluster_id", scales = "free_y", ncol = 4)+
      theme_classic()+
      stat_pvalue_manual(pwc, hide.ns = TRUE) +
      labs(
        subtitle = get_test_label(res.aov,  detailed = TRUE),
        caption = get_pwc_label(pwc)
      )
    
    plotList[[i]] = statPlot
    
    
    
  }
  
  if(write_files == TRUE){
    message("Writing statistics plots")
    write.xlsx(aov_df, paste0(rsl, "/", vc, "/02 Statistics/Analysis of variance.xlsx"))
    write.xlsx(pwc_df, paste0(rsl, "/", vc, "/02 Statistics/Pairwise comparisons.xlsx"))
    for(i in 1:length(plotList)){
      png(file = paste0(rsl, "/", vc, "/02 Statistics/MetaCluster_", i, ".png"), width = 15, height = 15, units = "in", res = 300, title = "Stat")
      print(plotList[[i]])
      dev.off()
    }
  }
  
  print(aov_df)
  print(pwc_df)
  
  if(plot_graph == TRUE){
    lapply(plotList,print)
  }
}}







  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# IV. WIP  ----
##define PlotAbundances function
  
.check_sce_new <- function(x, y = FALSE) {
  stopifnot(
    is(x, "SingleCellExperiment"), 
    !is.null(x$sample_id))
  if (y) 
    stopifnot(
      !is.null(x$cluster_id),
      !is.null(metadata(x)$cluster_codes))
}

.check_k_new <- function(x, k) {
  kids <- names(cluster_codes(x))
  if (is.null(k)) return(kids[1])
  stopifnot(length(k) == 1, is.character(k))
  if (!k %in% kids)
    stop("Clustering ", dQuote(k), 
         " doesnt't exist; valid are",
         " 'names(cluster_codes(x))'.")
  return(k)
}

.check_cd_factor_new <- function(x, y, n = 1) {
  if (is.null(y))
    return(TRUE)
  if (!is.null(n))
    stopifnot(length(y) == n)
  stopifnot(
    is.character(y), 
    all(y %in% names(colData(x))),
    !vapply(colData(x)[y], is.numeric, logical(1)))
  return(TRUE)
}


.check_pal_new <- function(x, n = 2) {
  if (is.null(x)) 
    return(TRUE)
  stopifnot(
    length(x) >= n,
    is.character(x))
  foo <- tryCatch(col2rgb(x),
                  error = function(e) {})
  if (is.null(foo)) {
    arg_nm <- deparse(substitute(x))
    stop(sprintf("'%s' is invalid.", arg_nm))
  }
  return(TRUE)
}

.get_shapes_new <- function(x, shape_by) {
  if (is.null(shape_by))
    return(NULL)
  # default shapes
  shapes <- c(16, 17, 15, 3, 7, 8) 
  n <- nlevels(x[[shape_by]])
  if (n > 18) {
    message(
      "At most 17 shapes are currently supported but ",
      n, " are required; setting 'shape_by' to NULL.")
    return(NULL)
  } else if (n > 6) {
    more_shapes <- setdiff(c(seq_len(16)-1, 18), shapes)
    shapes <- c(shapes, more_shapes[seq_len(n-length(shapes))])
  } else shapes <- shapes[seq_len(n)]
  return(shapes)
}

.scale_exprs_new <- function(x, margin = 1, q = 0.01) {
  if (!is(x, "matrix")) x <- as.matrix(x)
  qs <- c(rowQuantiles, colQuantiles)[[margin]]
  qs <- qs(x, probs = c(q, 1-q))
  qs <- matrix(qs, ncol = 2)
  x <- switch(margin,
              "1" = (x - qs[, 1]) / (qs[, 2] - qs[, 1]),
              "2" = t((t(x) - qs[, 1]) / (qs[, 2] - qs[, 1])))
  x[x < 0 | is.na(x)] <- 0
  x[x > 1] <- 1
  return(x)
}

.check_assay_new <- function(x, y) {
  stopifnot(
    length(y) == 1, 
    is.character(y),
    sum(y == assayNames(x)) == 1)
  return(TRUE)
}





##relative cluster abundance plots
plotAbundances_new <- function(x, k = "meta20", 
                               by = c("sample_id", "cluster_id"), 
                               group_by = "condition", shape_by = NULL,
                               col_clust = TRUE, 
                               distance = c(
                                 "euclidean", "maximum", "manhattan", 
                                 "canberra", "binary", "minkowski"), 
                               linkage = c(
                                 "average", "ward.D", "single", "complete", 
                                 "mcquitty", "median", "centroid", "ward.D2"),
                               k_pal = CATALYST:::.cluster_cols) {
  
  # check validity of input arguments
  by <- match.arg(by)
  .check_sce_new(x, TRUE)
  k <- .check_k_new(x, k)
  .check_cd_factor_new(x, group_by)
  .check_cd_factor_new(x, shape_by)
  .check_pal_new(k_pal)
  linkage <- match.arg(linkage)
  distance <- match.arg(distance)
  stopifnot(is.logical(col_clust), length(col_clust) == 1)
  
  shapes <- .get_shapes_new(x, shape_by)
  if (is.null(shapes)) shape_by <- NULL
  
  # ramp cluster color palette
  if (by == "sample_id") {
    nk <- nlevels(cluster_ids(x, k))
    if (length(k_pal) < nk)
      k_pal <- colorRampPalette(k_pal)(nk)
  }
  
  # get frequencies by cluster & sample
  ns <- table(
    cluster_id = cluster_ids(x, k), 
    sample_id = sample_ids(x))
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  
  # add relevant cell metadata
  m <- match(df$sample_id, x$sample_id)
  for (i in c(shape_by, group_by))
    df[[i]] <- x[[i]][m]
  
  if (by == "sample_id" && col_clust 
      && length(unique(df$sample_id)) > 1) {
    d <- dist(t(fq), distance)
    h <- hclust(d, linkage)
    o <- colnames(fq)[h$order]
    df$sample_id <- factor(df$sample_id, o)
  }
  
  assign("flowsomdf", df, envir = .GlobalEnv)
  
  if(any(groupOrder == md[[groupVar]])){
    df[[groupVar]] = factor(df[[groupVar]], levels = groupOrder) #change order of x-axis group
  }
  # specify shared aesthetics
  p <- ggplot(df, aes_string(y = "Freq")) +
    labs(x = NULL, y = "Proportion [%]") + 
    theme_bw() + theme(
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = NA, color = NA),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.key.height  =  unit(0.8, "lines"))
  
  switch(by,
         sample_id = p + (if (!is.null(group_by)) 
           facet_wrap(group_by, scales = "free_x")) +
           geom_bar(
             aes_string(x = "sample_id", fill = "cluster_id"), 
             position = "fill", stat = "identity") +
           scale_fill_manual("cluster_id", values = k_pal) +
           scale_x_discrete(expand = c(0, 0)) +
           scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
           theme(
             panel.border = element_blank(),
             panel.spacing.x = unit(1, "lines"))
         ,
         cluster_id = {
           p <- p + scale_shape_manual(values = shapes) + guides(
             col = guide_legend(order = 1, override.aes = list(size = 3)),
             shape = guide_legend(override.aes = list(size = 3)))
           if (is.null(group_by)) {
             p + geom_boxplot(aes_string(x = "cluster_id"), alpha = 0.2,
                              position = position_dodge(), outlier.color = NA) + 
               geom_point(aes_string("cluster_id", shape = shape_by),
                          position = position_jitter(width = 0.2))
           } else {
             if ("violin_plot" %in% rcaPlotType == TRUE){
               violinq = p + geom_violin(aes_string(x = group_by, 
                                                    color = group_by, fill = group_by),
                                         position = position_dodge(), alpha = 0.2, 
                                         outlier.color = NA, show.legend = FALSE) +
                 
                 geom_point(aes_string(x = group_by, 
                                       col = group_by, shape = shape_by),
                            position = position_jitter(width = 0.2)) +
                 
                 stat_summary(aes_string(x = group_by), fun.y = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.5, fatten = 1.5) +
                 
                 facet_wrap("cluster_id", scales = "free_y", ncol = 4)
               
               print(violinq)
               
             } 
             if ("box_plot" %in% rcaPlotType == TRUE){
               boxq = p + geom_boxplot(aes_string(x = group_by, 
                                                  color = group_by, fill = group_by),
                                       position = position_dodge(), alpha = 0.2, 
                                       outlier.color = NA, show.legend = FALSE) +
                 
                 geom_point(aes_string(x = group_by, 
                                       col = group_by, shape = shape_by),
                            position = position_jitter(width = 0.2)) +
                 
                 
                 facet_wrap("cluster_id", scales = "free_y", ncol = 4)
               
               boxq
             } 
             
           }
         }
  )
}


plotDR_new <- function(x, dr = NULL, 
                       color_by = "condition", facet_by = NULL, ncol = NULL,
                       assay = "exprs", scale = TRUE, q = 0.01, dims = c(1, 2),
                       k_pal = CATALYST:::.cluster_cols, 
                       a_pal = hcl.colors(10, "viridis")) {
  
  # check validity of input arguments
  stopifnot(
    is(x, "SingleCellExperiment"),
    .check_assay_new(x, assay),
    length(reducedDims(x)) != 0,
    is.logical(scale), length(scale) == 1,
    is.numeric(q), length(q) == 1, q >= 0, q < 0.5)
  .check_pal_new(a_pal)
  .check_cd_factor_new(x, facet_by)
  
  if (!is.null(ncol)) 
    stopifnot(is.numeric(ncol), length(ncol) == 1, ncol %% 1 == 0)
  
  if (is.null(dr)) {
    dr <- reducedDimNames(x)[1]
  } else {
    stopifnot(
      is.character(dr), length(dr) == 1, 
      dr %in% reducedDimNames(x))
  }
  stopifnot(is.numeric(dims), length(dims) == 2, 
            dims %in% seq_len(ncol(reducedDim(x, dr))))
  
  if (!all(color_by %in% rownames(x))) {
    stopifnot(length(color_by) == 1)
    if (!color_by %in% names(colData(x))) {
      .check_sce_new(x, TRUE)
      .check_pal_new(k_pal)
      .check_k_new(x, color_by)
      kids <- cluster_ids(x, color_by)
      nk <- nlevels(kids)
      if (length(k_pal) < nk)
        k_pal <- colorRampPalette(k_pal)(nk)
    } else kids <- NULL
  }
  
  # construct data.frame of reduced dimensions & relevant cell metadata
  xy <- reducedDim(x, dr)[, dims]
  colnames(xy) <- c("x", "y")
  df <- data.frame(colData(x), xy, check.names = FALSE)
  if (all(color_by %in% rownames(x))) {
    es <- as.matrix(assay(x, assay))
    es <- es[color_by, , drop = FALSE]
    if (scale) 
      es <- .scale_exprs_new(es, 1, q)
    df <- melt(
      cbind(df, t(es)), 
      id.vars = colnames(df))
    l <- switch(assay, exprs = "expression", assay)
    l <- paste0("scaled\n"[scale], l)
    scale <- scale_colour_gradientn(l, colors = a_pal)
    thm <- guide <- NULL
    color_by <- "value"
    facet <- facet_wrap("variable", ncol = ncol)
  } else if (is.numeric(df[[color_by]])) {
    if (scale) {
      vs <- as.matrix(df[[color_by]])
      df[[color_by]] <- .scale_exprs_new(vs, 2, q)
    }
    l <- paste0("scaled\n"[scale], color_by)
    scale <- scale_colour_gradientn(l, colors = a_pal)
    color_by <- sprintf("`%s`", color_by)
    facet <- thm <- guide <- NULL
  } else {
    facet <- NULL
    if (!is.null(kids)) {
      df[[color_by]] <- kids
      scale <- scale_color_manual(values = k_pal)
    } else scale <- NULL
    n <- nlevels(droplevels(factor(df[[color_by]])))
    guide <- guides(col = guide_legend(
      ncol = ifelse(n > 12, 2, 1),
      override.aes = list(alpha = 1, size = 3))) 
    thm <- theme(legend.key.height = unit(0.8, "lines"))
  }
  
  # set axes equal for linear dimension reductions
  if (dr %in% c("PCA", "MDS")) {
    asp <- coord_equal()
  } else asp <- NULL
  
  # get axes labels
  if (dr == "PCA") {
    labs <- paste0("PC", dims)
  } else labs <- paste(dr, "dim.", dims)
  
  # remove cells for which no reduced dimensions are available
  df <- df[!(is.na(df$x) | is.na(df$y)), ]
  
  #assign("umapdf", df, envir = .GlobalEnv)
  
  #q <- ggplot(df, aes_string("x", "y")) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour = "white") +
  #xlim(min(df$x)-1.3, (max(df$x)+1.3))+
  #ylim(min(df$y)-1.3, (max(df$y)+1.3))+
  ##coord_cartesian(xlim=c((min(df$x)-1.3), (max(df$x)+1.3)), ylim=c((min(df$y)-1.3), (max(df$y)+1.3))) +
  
  #scale_fill_distiller(palette = "RdYlBu") +
  #labs(x = labs[1], y = labs[2]) +
  #facet + scale + guide + asp + 
  
  #theme_minimal() + thm + theme(
  # panel.grid.minor = element_blank(),
  #strip.text = element_text(face = "bold"),
  #axis.text = element_text(color = "black"),
  #aspect.ratio = if (is.null(asp)) 1 else NULL)
  
  
  p <- ggplot(df, aes_string("x", "y", col = color_by)) +
    geom_point(size = 0.4, alpha = 0.8) + 
    labs(x = labs[1], y = labs[2]) +
    xlim(min(df$x)-1.3, (max(df$x)+1.3)) +
    ylim(min(df$y)-1.3, (max(df$y)+1.3)) +
    #coord_cartesian(xlim=c((min(df$x)-1.3), (max(df$x)+1.3)), ylim=c((min(df$y)-1.3), (max(df$y)+1.3))) +
    facet + scale + guide + asp + 
    theme_minimal() + thm + theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      aspect.ratio = if (is.null(asp)) 1 else NULL)
  
  
  
  if (is.null(facet_by)) 
    return(p)
  
  if (is.null(facet)) {
    p + facet_wrap(facet_by)
    #q + facet_wrap(facet_by)
  } else {
    if (nlevels(df$variable) == 1) {
      p + facet_wrap(facet_by, ncol = ncol) + 
        ggtitle(levels(df$variable))
      #print(q + facet_wrap(facet_by, ncol = ncol) + 
      #ggtitle(levels(df$variable)))
    } else {
      fs <- c("variable", facet_by)
      ns <- vapply(df[fs], nlevels, numeric(1))
      if (ns[2] > ns[1]) fs <- rev(fs)
      p + facet_grid(reformulate(fs[1], fs[2]))
      #print(q + facet_grid(reformulate(fs[1], fs[2])))
    }
  }
} 


plotClusterHeatmap_new <- function(x, hm2 = NULL, 
                                   k = "meta20", m = NULL, fun = c("median", "mean"), 
                                   cluster_anno = TRUE, split_by = NULL, scale = TRUE, 
                                   draw_dend = TRUE, draw_freqs = FALSE, 
                                   palette = rev(brewer.pal(11, "RdYlBu"))) {
  
  .Deprecated(
    new = "plotMultiHeatmap",
    old = "plotClusterHeatmap",
    msg = paste(sep = "\n",
                "'plotClusterHeatmap' is deprecated; instead, please use",
                " o 'plotExprHeatmap' for aggregated expression heatmaps",
                " o 'plotFreqHeatmap' for cluster frequency heatmaps",
                " o 'plotMultiHeatmap' to combine multiple heatmaps"))
  
  if (is.null(hm2)) {
    plotExprHeatmap(x, features = TSMarkers, 
                    by = "cluster_id", k = k, m = m,
                    assay = "exprs", fun = match.arg(fun),
                    scale = "first", q = 0.01,
                    row_anno = cluster_anno, col_anno = FALSE,
                    row_clust = TRUE, col_clust = FALSE,
                    row_dend = TRUE, col_dend = FALSE,
                    bars = draw_freqs, perc = draw_freqs, 
                    hm_pal = palette)
  } else {
    plotMultiHeatmap(x, 
                     hm1 = "type", hm2 = hm2, 
                     k = k, m = m, 
                     assay = "exprs", fun = match.arg(fun), 
                     scale = ifelse(scale, "first", "never"), 
                     q = 0.01, normalize = FALSE,
                     row_anno = cluster_anno, col_anno = FALSE, 
                     row_clust = TRUE, col_clust = FALSE, 
                     row_dend = TRUE, col_dend = FALSE, 
                     bars = draw_freqs, perc = draw_freqs, 
                     hm1_pal = palette)
  }
}


plotExprs_new <- function(x, features = NULL, 
                          color_by = "condition", assay = "exprs") {
  # check validity of input arguments
  .check_sce_new(x)
  .check_assay_new(x, assay)
  .check_cd_factor_new(x, color_by)
  
  # subset features to use
  features <- .get_features(x, features)
  y <- assay(x, assay)[features, ]
  
  # construct data.frame include cell metadata
  df <- data.frame(t(y), colData(x), check.names = FALSE)
  value <- ifelse(assay == "exprs", "expression", assay)
  
  gg_df <- melt(df, 
                value.name = value,
                variable.name = "antigen", 
                id.vars = names(colData(x)))
  
  ggplot(gg_df, fill = NULL, 
         aes_string(
           x = value, y = "..ndensity..",
           col = color_by, group = "sample_id")) + 
    facet_wrap(~ antigen, scales = "free_x") +
    geom_density() + 
    ylab("normalized density") +
    theme_classic() + theme(
      panel.grid = element_blank(), 
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black"), 
      axis.title = element_text(color = "black"))
}

# x: SingleCellExperiment
# y: logical; should cluster() have been run?
#' @importFrom methods is
#' @importFrom S4Vectors metadata
.check_sce_new <- function(x, y = FALSE) {
  stopifnot(
    is(x, "SingleCellExperiment"), 
    !is.null(x$sample_id))
  if (y) 
    stopifnot(
      !is.null(x$cluster_id),
      !is.null(metadata(x)$cluster_codes))
}


runStat_old = function(){
  setwd(paste0(rsl, "/", vc))
  dir.create("02 Statistics")
  dfStat = read.xlsx(paste0(rsl, "/", vc, "/FlowSOM_data_", customMeta, "MCs.xlsx"), 1, stringsAsFactors = FALSE) #import data for statistics
  dfStat$cluster_id = factor(dfStat$cluster_id, levels = c(1:length(unique(dfStat$cluster_id)))) #change order of cluster_IDs for faceting
  dfStat[[groupVar]] = factor(dfStat[[groupVar]], levels = groupOrder) #change order of x-axis group
  patientIDlist = c()
  for(i in unique(dfStat$sample_id)){ #append patient_id's
    sampleIDindex = which(md$sample_id == i)
    patientIDlist = append(patientIDlist, rep(md$patient_id[sampleIDindex], customMeta))
  }
  dfStat$patient_id = patientIDlist
  
  statGroup = colnames(dfStat[ncol(dfStat)-1]) #get the x-axis group
  statPlot = ggplot(dfStat, aes_string(y="Freq", x=statGroup))+
    geom_boxplot(aes_string(color = statGroup, fill = statGroup),
                 position = position_dodge(), alpha = 0.2, 
                 outlier.color = NA, show.legend = FALSE) +
    geom_point(aes_string(color = statGroup),
               position = position_jitter(width = 0.2)) +
    facet_wrap("cluster_id", scales = "free_y", ncol = 4)+
    theme_classic()
  
  if(length(unique(dfStat[[statGroup]]))<3){ #statistics on 2 groups
    if(pairedData == TRUE){
      cat("\n")
      message(paste0("Performing statistics on two paired groups using ", PWC,". Showing only significant results in plots"))
    }else{
      cat("\n")
      message(paste0("Performing statistics on two unpaired groups using ", PWC,". Showing only significant results in plots"))
    }
    
    
    for(i in unique(dfStat$cluster_id)){
      
    }
    
    
    statPlot = statPlot + stat_compare_means(paired = pairedData, method = PWC, label = "p.signif")
    
    if(write_files == TRUE){
      message("Writing statistics plots")
      png(file = paste0(rsl, "/", vc, "/02 Statistics/Statistics.png"), width = 15, height = 15, units = "in", res = 300, title = "Stat")
      print(statPlot)
      dev.off()
    }
    
    if(plot_graph == TRUE){
      print(statPlot)
    }
    
  }else if(pairedData == FALSE & length(unique(dfStat[[statGroup]]))>2){ #statistics on unpaired >2 groups
    cat("\n")
    message(paste0("Performing statistics on multiple unpaired groups using ", AOV,". Showing only significant results in plots"))
    
    #create list with index numbers of groups to be compared
    multCompList = vector(mode = "list")
    tempVec = c()
    holdList = c(1:length(unique(dfStat[[statGroup]])))
    statItr = 0
    for(i in 1:length(unique(dfStat[[statGroup]]))){
      holdList = holdList[holdList != i]
      
      for(j in holdList){
        statItr = statItr + 1
        tempVec[1] = i
        tempVec[2] = j
        multCompList[[statItr]] = tempVec
      }
    }
    
    statPlot = statPlot + stat_compare_means(method = AOV, paired = FALSE, label.y.npc = 0.5, comparisons = multCompList)
    
    if(write_files == TRUE){
      message("Writing statistics plots")
      png(file = paste0(rsl, "/", vc, "/Statistics/Statistics.png"), width = 15, height = 15, units = "in", res = 300, title = "Stat")
      print(statPlot)
      dev.off()
    }
    
    if(plot_graph == TRUE){
      print(statPlot)
    }
    
  }else if(pairedData == TRUE & length(unique(dfStat[[statGroup]]))>2){ #statistics on paired >2 groups
    cat("\n")
    message(paste0("Performing statistics on multiple paired groups using Friedman test with wilcoxon for multiple comparisons and ", multCompAdjust, " correction"))
    message("Showing only significant results in plots")
    message(paste0("Removing all samples that are not present in the ", statGroup, " group"))
    pID.table = table(dfStat$patient_id) #tabulate the occurences of patien_id's
    pID.max.occurence = max(pID.table) #find which donors occur the most times
    #pID.max.occurrence = length(unique(dfStat[[statGroup]]))
    dfStat = dfStat[which(dfStat$patient_id %in% c(names(pID.table)[pID.table == pID.max.occurence])),] #remove samples that occur less than others from data frame
    dfStat$patient_id = factor(dfStat$patient_id)
    dfStat[[statGroup]] = factor(dfStat[[statGroup]])
    
    dfStatHold = dfStat
    plotList = list()
    friedman_df = data.frame()
    pwc_df = data.frame()
    for(i in unique(dfStat$cluster_id)){
      dfStat = subset(dfStatHold, cluster_id == i)
      res.friedman = dfStat %>% friedman_test(as.formula(paste0("Freq ~ ",statGroup, "| patient_id")))
      temp_friedman_df = as.data.frame(res.friedman)
      temp_friedman_df$cluster_id = c(i)
      friedman_df = rbind(friedman_df, temp_friedman_df)
      
      pwc = dfStat %>%
        wilcox_test(as.formula(paste0("Freq ~ ", statGroup)), paired = TRUE, p.adjust.method = multCompAdjust)
      temp_pwc_df = as.data.frame(pwc)
      temp_pwc_df$cluster_id = c(rep(i, nrow(temp_pwc_df)))
      pwc_df = rbind(pwc_df, temp_pwc_df)
      
      pwc <- pwc %>% add_xy_position(x = statGroup)
      
      statPlot =  ggplot(dfStat, aes_string(y="Freq", x=statGroup))+
        geom_boxplot(aes_string(color = statGroup, fill = statGroup),
                     position = position_dodge(), alpha = 0.2, 
                     outlier.color = NA, show.legend = FALSE) +
        geom_point(aes_string(color = statGroup),
                   position = position_jitter(width = 0.2)) +
        facet_wrap("cluster_id", scales = "free_y", ncol = 4)+
        theme_classic()+
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.friedman,  detailed = TRUE),
          caption = get_pwc_label(pwc)
        )
      
      plotList[[i]] = statPlot
      
      
      
    }
    
    if(write_files == TRUE){
      message("Writing statistics plots")
      write.xlsx(friedman_df, paste0(rsl, "/", vc, "/Statistics/Friedman_test.xlsx"))
      write.xlsx(pwc_df, paste0(rsl, "/", vc, "/Statistics/Wilcoxon_signed_rank_test.xlsx"))
      for(i in 1:length(plotList)){
        png(file = paste0(rsl, "/", vc, "/Statistics/MetaCluster_", i, ".png"), width = 15, height = 15, units = "in", res = 300, title = "Stat")
        print(plotList[[i]])
        dev.off()
      }
    }
    
    print(friedman_df)
    print(pwc_df)
    
    if(plot_graph == TRUE){
      print(friedman_df)
      print(pwc_df)
      lapply(plotList,print)
    }
    
  }
}















