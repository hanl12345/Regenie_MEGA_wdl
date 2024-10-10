version 1.0

## Version 0719_2024
## This WDL workflow runs Regenie.
## This workflow assumes users have thoroughly read the Regenie docs for caveats and details.
## Regenie's documentation: https://rgcgithub.github.io/regenie/options/
#Most important setting
#Analyzed population list:EUR,AFR
#Analyzed chromosomes list
#1:23 for EUR and AFR
#Input file list,use EUR as an example, 
#Two kinds of genetic file:
#gs://fc-secure-52fe426d-2d5f-48f5-b819-de205040508b/GWASResult/MEGA_v4

#MEGA whole-genome dataset
#AFR
#gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/AFR/AFR_HighQuality_wg.bgen
#EUR
#gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/EUR/EUR_HighQuality_wg.bgen

#MEGA AFR per-chromosome imputed bgen files 
#gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/AFR/R203_mac10/AFR-chr*-mac10.bgen
#gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/AFR/R203_mac200/AFR-chr*-mac200.bgen

#MEGA EUR per-chromosome imputed bgen files
#gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/EUR/R203_mac10/EUR-chr*-mac10.bgen
#gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/EUR/R203_mac200/EUR-chr*-mac200.bgen

#Phecode file
#gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/megaex_v4_phecode_noExclusion_2024_primary.txt
#Covariate file
#gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/megaex_v4_cov_2024_primary.txt

## Cromwell version support - Successfully tested on v77
##
## Distributed under terms of the MIT License


workflow Regenie {
    # Refer to Regenie's documentation for the descriptions to most of these parameters.
    input { 
        String SuperPop="AFR" 
        Int macFile=200
        String SourceDir="gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409"
        
         String geneticFileFormat="bgen"  

         String outfile_prefix="gs://fc-secure-52fe426d-2d5f-48f5-b819-de205040508b/GWASResult/MEGA_v4"
         #--bsize 1000 --cv 10 --bt --spa --pThresh 0.05;--bsize 1000 --bt --spa --pThresh 0.05
         #--bsize 1000 --cv 10 --qt;--bsize 1000 --qt
         
         String step1_Options ="--bsize 1000 --cv 10 --bt --firth --approx --pThresh 0.05 --firth-se "
         String step2_Options ="--bsize 1000 --bt --firth --approx --pThresh 0.05 --firth-se "
         
         
         File phenoFile ="gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/megaex_v4_phecode_noExclusion_2024_primary.txt"
 
         File covarFile ="gs://working-set-redeposit/megaex/megaex_regenie/megaex_regenie_v4_202409/megaex_v4_cov_2024_primary.txt"
 
         String? covarColList  #added

        Int? minCaseCount=500 #added on 0901_2023
        String? test   #added on 0901_2023
        Float? minMAC=200   #added on 0901_2023
                 
        Int threads=64 #Added 0n 09/01/2023
               
        Array[String] chr_list=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"] # List of chromosomes for Step2.
       #Array[Int] chr_list=[21,22]
       Array[String]+ phenoColList # Phenotypes you want to analyze. (Column names).
                        
        String regenie_docker = "ghcr.io/rgcgithub/regenie/regenie:v3.6.gz" #revise 1003_2024
        
        #String regenie_docker = "briansha/regenie:v3.0_boost" # Compiled with Boost IOSTREAM: https://github.com/rgcgithub/regenie/wiki/Using-docker
        
        #Int preemptible = 1 #added on 0901_2023
        Int maxRetries = 0 #added on 0901_2023
      }
       

    call fitModel {
        input:
            SourceDir=SourceDir,
            #geneticFilePrefix_t1=geneticFilePrefix_Step1,
            SuperPop=SuperPop,
            geneticFileFormat=geneticFileFormat,#Add it on 0912_2023
            
            outfile_prefix = outfile_prefix, #Added on 0902_2023
            covarFile = covarFile,                            
            phenoFile = phenoFile, 
            docker = regenie_docker,
            step1_Options=step1_Options,
            
            covarColList=covarColList,    #added        
            phenoColList=phenoColList, #Added it on 09032023             
            minCaseCount=minCaseCount, #added on 0901_2023
            test =test,   #added on 0901_2023
            minMAC = minMAC,  #added on 0901_2023,remved on 12_15_2023, no useful for step1
            threads=threads,  #added on 0901_2023
             #preemptible =preemptible,  #added on 0901_2023
            maxRetries =maxRetries  #added on 0901_2023            
            
    }

    scatter (chromosome in chr_list) {
      call chrAssociationTest {
          input:
          	  SourceDir=SourceDir,
              #geneticFilePrefix_t2 = ChrFileDir_Step2,
              SuperPop=SuperPop,
              macFile=macFile,
              step2_Options=step2_Options,
              geneticFileFormat=geneticFileFormat, #Added it on 0912_2023
              outfile_prefix = outfile_prefix, #Added on 0902_2023
              chr = chromosome, # May need to change to chrList = chr_list later, as chr is Int only...not accounting for X, Y, etc.
              covarFile = covarFile,
              phenoFile = phenoFile,
              pred = fitModel.fit_bin_out,
              #use_null_firth=fitModel.firth_list,  #Add it on 12_12_2023
              docker = regenie_docker,
              output_locos_gz = fitModel.output_locos_gz,
              covarColList=covarColList,    #added
             
              phenoColList=phenoColList,#   removed it on 09022023          
              
              minCaseCount=minCaseCount, #added on 0901_2023
              test =test,   #added on 0901_2023
              minMAC = minMAC,  #added on 0901_2023
            
              threads=threads,  #added on 0901_202
              
             # preemptible =preemptible,  #added on 0901_2023
              maxRetries =maxRetries  #added on 0901_2023
      }
    }

    call joinGWASResults {
      input:
        output_files = chrAssociationTest.test_bin_out_firth, # Refers implicitly to the entire array of files that were scattered.
        SuperPop=SuperPop,
        outfile_prefix = outfile_prefix,#Added on 0902_2023
        chr_list = chr_list,       
        phenoColList = phenoColList,        
        maxRetries =maxRetries  #added on 0903_2023
    }

    call drawManhattanQQPlots {
      input:
        phenoColList=phenoColList, #Added it on 09032023
        SuperPop=SuperPop,
        outfile_prefix = outfile_prefix, #Added on 0902_2023
        chr_list = chr_list,
        file_input = joinGWASResults.outputs,
        maxRetries =maxRetries  #added on 0903_2023
    }
    
    #added it on 0909_2023
    call moveFiles2TargetDir {
         input: 
          outfile_prefix = outfile_prefix, #Added on 0902_2023 
          SuperPop=SuperPop,
          phenoColList=phenoColList, #Added it on 09032023
          FullRegenieFile = drawManhattanQQPlots.output_regenie, 
          output_qq_plots = drawManhattanQQPlots.output_qq_plots,
          output_manhattan_plots = drawManhattanQQPlots.output_manhattan_plots,
          #Run time

          maxRetries =maxRetries,  #added on 0903_2023
       }
       #revised it on 0909_2023
    output {
           File allOutfilesList=moveFiles2TargetDir.bucket_outprefix          
       }
    meta {
    	author : "Lide Han"
        email : "lide.han@vumc.org"
        description : "This workflow runs Regenie"
    }
}

# In Step 1, the whole genome regression model is fit to the traits
# and a set of genomic predictions are produced as output.
# Refer to Regenie's documentation for the descriptions to most of these parameters.
task fitModel { 
    input {     
        String outfile_prefix #Added it on 0902_2023
        String SuperPop
        String geneticFileFormat
        String step1_Options

        String SourceDir                 
        File? phenoFile
        File? covarFile
        String? covarColList
        
        Int? minCaseCount

        String? test
        Float? minMAC
        
        Int threads  #added on 0901_2023      
             
        Array[String] phenoColList #Added it on 0903_2023   

        # Runtime
        String docker
        Int maxRetries       
    }
    
     String geneticFilePrefix_t1=SourceDir + "/"+ SuperPop + "/" + SuperPop + "_HighQuality_wg"
     String outfile_prefix_1=outfile_prefix + SuperPop
    #Int cpu=threads
     Int cpu=threads/2
    Int phecodelen = length(phenoColList) #
     #Revised it 0912_2023
    File familyFile=geneticFilePrefix_t1+ (if geneticFileFormat == "bed" then ".fam" else if geneticFileFormat == "bgen" then ".bgen" else ".psam") #revised it on 0912_2023
 
    File bedFile=geneticFilePrefix_t1+ (if geneticFileFormat == "bed" then ".bed" else if geneticFileFormat == "bgen" then ".bgen" else ".pgen")
    
    File bimFile=geneticFilePrefix_t1 + (if geneticFileFormat == "bed" then ".bim" else if geneticFileFormat == "bgen" then ".bgen" else ".pvar")
     
    Float genotype_size = size(bedFile, "GiB")
    
    #File extractFile
    
    Int disk =ceil((1+phecodelen) * genotype_size)#1.4+0.20
    Float memory =ceil((1.4+phecodelen*0.15) * genotype_size) #edit it on 11_24_2023,revised it on 0728_2024 from 1 to 2
     
    String OutfilePrefix= basename(outfile_prefix_1) #Added it on 0902_2023

    command <<<
        set -euo pipefail  
        echo "disk ~{disk} gb"
         echo "memory ~{memory} gb"
         echo "cpu: ~{cpu} "
         
        Step1flags='--threads=~{threads} --gz '
        Step1flags+='--phenoFile=~{phenoFile} '
        Step1flags+='~{step1_Options} '
        Step1flags+='--phenoColList=~{sep=',' phenoColList} '
        Step1flags+='~{if defined(covarFile) then "--covarFile=~{covarFile} " else " "} '
        Step1flags+='~{if defined(covarColList) then "--covarColList=~{covarColList} " else " "} '
        Step1flags+='~{if defined(test) then "--test=~{test} " else " "} '      
       
       if [[ ~{geneticFileFormat} == "bed" ]]; then       
        regenie --step 1 \
        --bed="~{sub(bedFile,'\\.bed$','')}" \
        ${Step1flags} \
        --out=~{OutfilePrefix}
        elif [[~{geneticFileFormat} == "pgen" ]]; then 
        regenie --step 1 \
        --pgen="~{sub(bedFile,'\\.pgen$','')}" \
        ${Step1flags} \
        --out=~{OutfilePrefix}
        else
        regenie --step 1 \
        --bgen="~{bedFile}" \
        ${Step1flags} \
        --out=~{OutfilePrefix}
        fi               
    >>>

    output {
        File fit_bin_out = "${OutfilePrefix}_pred.list" # Refers to the list of loco files written. Lists the files as being located in the current working directory.      
        Array[File] output_locos_gz=glob("*.loco.gz")
   }

    runtime {
        docker:docker # Compiled with Boost IOSTREAM: https://github.com/rgcgithub/regenie/wiki/Using-docker
        cpu: cpu
        memory: memory + " GiB"
	    disks: "local-disk " + disk + " HDD"
        #preemptible: preemptible
        maxRetries: maxRetries
    }
}

# In Step 2, a set of imputed SNPs are tested for association using a Firth logistic regression model.
# Refer to Regenie's documentation for the descriptions to most of these parameters.
task chrAssociationTest {
    
    input {
    
         String SourceDir 
         #String  geneticFilePrefix_t2
         String SuperPop
         String geneticFileFormat
         Int macFile
         
         String step2_Options
        
         Array[String] phenoColList #Added it on 0903_2023
         
        String outfile_prefix
        Array[File] output_locos_gz

        # Basic options
        File? phenoFile

        File? covarFile
        String? covarColList
        File? pred     # File containing predictions from Step 1

        # Options         
        Int? minCaseCount
        String? test
        Float? minMAC
        String chr  #delete on 0903_2023        
        Int threads  #added on 0901_2023 #   removed it on 10032024 

        # Runtime
        String docker
        Int maxRetries
    }
    
     String  geneticFilePrefix_t2= SourceDir + "/"+ SuperPop + "/R203_mac" + macFile + "/"
    
    String outfile_prefix_1=outfile_prefix + SuperPop
    String EachChrBaseFilePrefix=geneticFilePrefix_t2+SuperPop+"-chr"+chr+"-mac"+macFile  
    File familyFile=EachChrBaseFilePrefix + (if geneticFileFormat == "bed" then ".fam" else if geneticFileFormat == "bgen" then ".sample" else ".psam") #revised it on 0912_2023
 
    File bedFile=EachChrBaseFilePrefix + (if geneticFileFormat == "bed" then ".bed" else if geneticFileFormat == "bgen" then ".bgen" else ".pgen")
    
    File bimFile=EachChrBaseFilePrefix + (if geneticFileFormat == "bed" then ".bim" else if geneticFileFormat == "bgen" then ".bgen" else ".pvar")
   
     Int phecodelen = length(phenoColList) #
     Float dosage_size = size(bedFile, "GiB") 
   
    Int disk = ceil((1.0+phecodelen) * dosage_size) #edit it on 11_24_2023  
    
     Int cpu=2
      Float memory= ceil(1.1 * dosage_size)
     String OutfilePrefix=basename(outfile_prefix_1) + ".Chr" + chr #Added it on 0901_2023
     
      
# Loco files are moved to the current working directory due to the list of predictions (pred) expecting them to be there.
    command <<<
           
        echo "disk:~{disk} gb"
        echo "membory: ~{memory} gb"
        set -euo pipefail
        
        for file in ~{sep=' ' output_locos_gz}; do \
              mv $file .; \
         done
                      
         Step2flags='--threads=~{threads} '
         Step2flags+='~{step2_Options} '
         Step2flags+='--phenoFile=~{phenoFile} '
         Step2flags+='--phenoColList=~{sep=',' phenoColList} '
         Step2flags+='--chr=~{chr} '
         Step2flags+='~{if defined(covarFile) then "--covarFile=~{covarFile} " else " "} '
         Step2flags+='~{if defined(covarColList) then "--covarColList=~{covarColList} " else " "} '
         Step2flags+='~{if defined(pred) then "--pred=~{pred} " else " "} '
         Step2flags+='~{if defined(test) then "--test=~{test} " else " "} '
         Step2flags+='~{if defined(minMAC) then "--minMAC=~{minMAC} " else " "} '
        
        if [[ ~{geneticFileFormat} == "bed" ]]; then
        regenie --step 2 \
        --bed=~{sub(bedFile,'\\.bed$','')} \
        ${Step2flags} \
        --out=~{OutfilePrefix}
        elif [[~{geneticFileFormat} == "pgen" ]]; then 
        regenie --step 2 \
        --pgen=~{sub(bedFile,'\\.pgen$','')} \
        ${Step2flags} \
        --out=~{OutfilePrefix}
        else
        regenie --step 2 \
        --bgen=~{bedFile} \
        ${Step2flags} \
        --out=~{OutfilePrefix}
        fi   
    >>>

    output {
        Array[File] test_bin_out_firth = glob("*.regenie")
     }

    runtime {
        docker:docker # Compiled with Boost IOSTREAM: https://github.com/rgcgithub/regenie/wiki/Using-docker
        memory: memory + " GiB"
	    disks: "local-disk " + disk + " HDD"
        cpu: cpu
        #preemptible: preemptible
        maxRetries: maxRetries
    }
}


# Join together all output relating to the scattered tasks in Step 2.
# For each phenotype, one file will be made containing all analysis from Step 1 on each chromosome used in testing.
task joinGWASResults {
  input {
      Array[Array[File]] output_files # All output from Step 2.
      String SuperPop
      Array[String] phenoColList
      String outfile_prefix #added it on 0903_2023
      Array[String] chr_list
      #String docker

      Int maxRetries = 0
  }
  
  String outfile_prefix_1=outfile_prefix + SuperPop
  
  Array[File] all_output_files = flatten(output_files)
  String OutfileName= basename(outfile_prefix_1, ".bgen") #Added it on 0903_2023
  Float regenie_files_size = size(all_output_files, "GiB")
  
  Int phecodelen = length(phenoColList)
  Int disk = ceil(30.0 + (1.0+phecodelen) * regenie_files_size)
  Float memory= ceil(1.5 * regenie_files_size+2)  #Added it on 11_21_2023


  command <<<
        set -euo pipefail
        echo "disk ~{disk} gb"
        echo "memory ~{memory} gb"
        echo "OutfileName: ~{OutfileName}"
        for array in ~{sep=' ' all_output_files}; do 
          for file in $array; do 
            mv $file .; 
        done done
        
        for pheno in ~{sep=' ' phenoColList}; do 
        for chr in ~{sep= ' ' chr_list}; do 
            cat ~{OutfileName}.Chr${chr}_${pheno}.regenie >> ~{OutfileName}_${pheno}.regenie; 
            rm ~{OutfileName}.Chr${chr}_${pheno}.regenie; 
        done done
        
        for pheno in ~{sep=' ' phenoColList}; do 
        gzip -fc ~{OutfileName}_${pheno}.regenie > ~{OutfileName}_${pheno}.regenie.gz; 
        done; 

  >>>

  output {
        Array[File] outputs=glob("*.regenie.gz")
   }

  runtime {
        docker: "briansha/regenie_r_base:4.1.0"
        memory: memory + " GiB"
	    disks: "local-disk " + disk + " HDD"
       # cpu: cpu
       # preemptible: preemptible
        maxRetries: maxRetries
  }
}


# QQ and Manhattan plots
task drawManhattanQQPlots {
  input {
  
     String outfile_prefix #Added it on 0903_2023
     String SuperPop
     Array[String] phenoColList #Added it on 0903_2023    
    #Array[String] phenoColList_regenie # Format: "<phenotype_name>.regenie". Names of the files produced in joinGWASResults.
    Array[String] chr_list
    Array[File] file_input
    #String docker

    Int maxRetries = 0
  }
   
   String outfile_prefix_1=outfile_prefix + SuperPop
   
   String OutfileName= basename(outfile_prefix_1, ".bgen") + "_" #Added it on 0903_2023
   
   Array[String] sourcephenoColList=prefix(OutfileName, phenoColList) #Added it on 0903_2023
   
   String finalOutfileName= sub(outfile_prefix_1,".bgen","") + "_" #Added it on 0903_2023   
   Array[String] full_phenotypes=prefix(finalOutfileName, phenoColList) #Added it on 0903_2023
   Int totalNumber = length(phenoColList)
   Int phecodelen = length(phenoColList)
   Boolean gz = true #Added it on 0910_2023 
     
  Float regenie_files_size = size(file_input, "GiB")
  #68
   Int disk =  ceil(20.0 +(3.0+phecodelen)*regenie_files_size)
  Float memory= ceil(2.0* regenie_files_size+20.0)  #Added it on 11_21_2023,change into 20, 68
  
  command <<<
    set -euo pipefail
     echo "disk ~{disk} gb"
     echo "memory ~{memory} gb"
     
    file_input_Array=(~{sep=' ' file_input});
    sourcephenotype_Array=(~{sep=' ' sourcephenoColList}); 
    for((c = 0; c < ~{totalNumber}; c++)); do 
         mv ${file_input_Array[$c]} .; 
    done 
    R --no-save --args ~{sep=' ' sourcephenoColList} ~{gz} <<RSCRIPT
    library(data.table)
    library(qqman)
    args <- commandArgs(trailingOnly = TRUE)
    NumberofArgs=length(args)
    gzStatus=args[NumberofArgs]
    filelist=head(args,-1)
    for (file in filelist) {   
        if(gzStatus)
        {
           file1=paste(file,".regenie.gz",sep="")
         }
         else
         {
           file1=paste(file,".regenie",sep="")
         }      
        regenie_output <- fread(file1)
        regenie_ADD_subset <-subset.data.frame(regenie_output, TEST=="ADD")
        regenie_ADD_subset=as.data.frame(regenie_ADD_subset);#data.frame conversion
        print(colnames(regenie_ADD_subset));
        regenie_ADD_subset[,"CHROM"] <-as.numeric(unlist(regenie_ADD_subset[,"CHROM"]))
        regenie_ADD_subset[,"LOG10P"] <-as.numeric(unlist(regenie_ADD_subset[,"LOG10P"]))
        regenie_ADD_subset[,"GENPOS"] <-as.numeric(unlist(regenie_ADD_subset[,"GENPOS"]))
         # Check for non-finite values,Ensure your data does not contain non-finite values
         if(any(!is.finite(regenie_ADD_subset[,"LOG10P"]))) {
          regenie_ADD_subset <- regenie_ADD_subset[is.finite(regenie_ADD_subset[,"LOG10P"]), ];  # Remove non-finite values
         }
        qq_plot = paste(file, "qqplot.png", sep="_")
        png(qq_plot)
        p = 10 ^ (-1 * (as.numeric(unlist(regenie_ADD_subset[,"LOG10P"]))))
        qq(p,main = "Q-Q plot of GWAS p-values",pch = 18, col = "blue4", cex = 1.5, las = 1);
        dev.off()       
        manhattan_plot = paste(file, "manhattan.png", sep="_")
        png(manhattan_plot)
        manhattan(regenie_ADD_subset, chr="CHROM", bp="GENPOS", snp="ID",p="LOG10P", logp=FALSE, genomewideline = -log10(5e-08),suggestiveline=FALSE)
        dev.off()               
      }
    RSCRIPT
  >>>

  output {
         Array[File] output_qq_plots = glob("*qqplot.png")
         Array[File] output_manhattan_plots = glob("*manhattan.png")
         Array[File] output_regenie = glob("*.regenie.gz")
  }


  runtime {
        docker: "briansha/regenie_r_base:4.1.0"  # Ubuntu 18.04, R 4.1.0, and a few Ubuntu and R packages.
        memory: memory + " GiB"
	    disks: "local-disk " + disk + " HDD"
        #preemptible: preemptible
        maxRetries: maxRetries
  }
}

#Array[File]+ PlotsFile
task moveFiles2TargetDir {
     input {
            Array[String]+ phenoColList 
            String outfile_prefix
            String SuperPop
            Array[File]+  FullRegenieFile
            Array[File]+ output_qq_plots   #Revised it on 0905_2023
            Array[File]+ output_manhattan_plots  #Revised it on 0905_2023
            #run time
            #String docker
            Int maxRetries
         }
         
         Float memory=10
         String outfile_prefix_1=outfile_prefix + SuperPop
         
         String OutfileName= basename(outfile_prefix_1, ".bgen") + "_" #Added it on 0903_2023   
         Array[String] base_phenotypes=prefix(OutfileName, phenoColList) #Added it on 0903_2023
         #Remove the last string    
         String finalOutfileName= sub(outfile_prefix_1,".bgen","") + "_" #Added it on 0903_2023   
         Array[String] full_phenotypes=prefix(finalOutfileName, phenoColList) #Added it on 0903_2023
         Int totalNumberofPhenos = length(phenoColList)
         Float regenie_files_size = size(FullRegenieFile, "GiB")
         Int disk =ceil(10.0 + regenie_files_size) 
        # Int disk = select_first([disk_size_override, ceil(10.0 + regenie_files_size)])   
      
         command <<<   
          set -euo pipefail
            echo "OutfilePrefix ~{outfile_prefix_1}"
            regenie_input_Array=(~{sep=' ' FullRegenieFile}); 
            qqplots_input_Array=(~{sep=' ' output_qq_plots});
            manhattanplots_input_Array=(~{sep=' ' output_manhattan_plots}); 
            targetphenotype_Array=(~{sep=' ' full_phenotypes});
            for((c = 0; c < ~{totalNumberofPhenos}; c++)); do 
           # out_prefix= ${targetphenotype_Array[$c]}
            gsutil mv ${qqplots_input_Array[$c]} ${targetphenotype_Array[$c]}_qqplot.png; 
            gsutil mv ${manhattanplots_input_Array[$c]} ${targetphenotype_Array[$c]}_manhattan.png; 
            echo "qqplot file: ${targetphenotype_Array[$c]}_qqplot.png"; 
            echo "manhattan file: ${targetphenotype_Array[$c]}_manhattan.png"; 
            done 
            
            for((c = 0; c < ~{totalNumberofPhenos}; c++)); do 
              out_prefix=~{outfile_prefix_1}${targetphenotype_Array[$c]}
              gsutil mv ${regenie_input_Array[$c]} ${targetphenotype_Array[$c]}.regenie.gz;
              echo "regenie output file (step 2): ${targetphenotype_Array[$c]}.regenie.gz";
            done 
           >>>
           
       runtime {
                  docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
                  memory: memory + " GiB"
	              disks: "local-disk " + disk + " HDD"
                 # preemptible: preemptible
                  maxRetries: maxRetries
           }
        output {
               File bucket_outprefix = stdout()
           }          
}