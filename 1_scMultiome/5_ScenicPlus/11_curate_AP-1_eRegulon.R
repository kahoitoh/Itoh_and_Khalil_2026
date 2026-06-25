#!/usr/bin/env Rscript

# ==============================================================================
# Extract AP-1 eRegulon from SCENIC+ result
# ------------------------------------------------------------------------------
# Purpose:
#   This script extracts AP-1-related eRegulon from the SCENIC+ output.
#
# Inputs:
#   - Manually curated eRegulon names associated with AP-1 motif activity.
#   - Candidate eRegulons are selected by inspecting motif logos in ctx_results.html.
#     This manual curation is required because SCENIC+ eRegulon names do not
#     always directly correspond to the displayed enriched motif. Therefore,
#     eRegulons with non-AP-1 TF names can still represent AP-1-like motif
#     enrichment.
#   - Only +/+ eRegulons are selected.
#
# Notes:
#   - Edit the CONFIG section before running.
#   - For CA1, "extended" result was used instead of "direct" result.
#   - For mPFC, Pvalb and Sst AP-1 eRegulons are concatenated and treated as
#     Inhibitory neuron AP-1 eRegulons.
#
# Outputs:
#   - AP1_eRegulon_<REGION>.csv
# ==============================================================================

library(dplyr)
library(purrr)

# config -----------------------------------------------------------------------
REGION <- c("BLA", "Hippo", "mPFC")[3]

SCENICPLUS_RES_DIR_LIST <- if(REGION == "BLA"){
  
  list("BA" = "path/to/ScenicPlus/res/BA",
       "LA_Chst9" = "path/to/ScenicPlus/res/LA_Chst9",
       "Sst" = "path/to/ScenicPlus/res/Sst")
  
}else if(REGION == "Hippo"){
  
  list("DG" = "path/to/ScenicPlus/res/DG",
       "CA1" = "path/to/ScenicPlus/res/CA1",
       "CA3" = "path/to/ScenicPlus/res/CA3")
  
}else if(REGION == "mPFC"){
  
  list("Exc" = "path/to/ScenicPlus/res/Exc", # result using c("L2_3_IT", "L4_5_IT", "L5_ET", "L5_NP", "L6_CT", "L6_IT")
       "Sst" = "path/to/ScenicPlus/res/Sst",
       "Pvalb" = "path/to/ScenicPlus/res/Pvalb")
  
}

AP1_eRegulon_NAME_LIST <- if(REGION == "BLA"){
  
  list("BA" = "Fosl2_direct_+/+_(265g)",
       "LA_Chst9" = c("Bach1_direct_+/+_(337g)", "Fosb_direct_+/+_(59g)", "Jund_direct_+/+_(487g)"),
       "Sst" = "Maf_direct_+/+_(187g)")
  
}else if(REGION == "Hippo"){
  
  list("DG" = c("Jdp2_direct_+/+_(38g)", "Jun_direct_+/+_(172g)", "Smad3_direct_+/+_(203g)"),
       "CA1" = c("Bcl11a_extended_+/+_(511g)"),
       "CA3" = c("Fosb_direct_+/+_(72g)", "Jund_direct_+/+_(319g)", "Maf_direct_+/+_(60g)", "Smad3_direct_+/+_(174g)"))
  
}else if(REGION == "mPFC"){
  
  list("Exc" = c("Tcf4_direct_+/+_(440g)", "Tcf12_direct_+/+_(399g)"),
       "Sst" = "Atf2_direct_+/+_(34g)",
       "Pvalb" = "Bach2_direct_+/+_(264g)")
  
}

FINAL_RES_DIR <- "path/to/your/final/res"


# import GRN result ------------------------------------------------------------
res_list <- lapply(SCENICPLUS_RES_DIR_LIST, 
                   function(p){

                     read.csv(paste0(p, "/df_direct_e_regulon_uns.csv")) 
                     
                   })

if(REGION == "Hippo"){
  
  res_list$CA1 <- read.csv(paste0(SCENICPLUS_RES_DIR_LIST$CA1, "/df_extended_e_regulon_uns.csv")) 
  
}

AP1_res_list <- map2(res_list, AP1_eRegulon_NAME_LIST,
                     function(res, AP1eReg){
                       
                       res[which(res$Gene_signature_name %in% AP1eReg),]
                       
                     })

AP1_res_list <- map2(AP1_res_list, names(AP1_res_list),
                     function(res, n){
                       
                       res$CellType <- n
                       res
                       
                     })

AP1_res <- rbind(AP1_res_list[[1]], AP1_res_list[[2]], AP1_res_list[[3]])

if(REGION == "mPFC"){
  
  AP1_res$CellType <- sub("Sst|Pvalb", "Inh", AP1_res$CellType)
  
}

write.table(AP1_res, paste0(FINAL_RES_DIR, "/AP1_eRegulon_", REGION, ".csv"), 
            sep = ",", quote = F, row.names = F)

