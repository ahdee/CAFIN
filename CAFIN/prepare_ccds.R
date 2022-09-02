library ( dplyr )
bp= 20 # what length you want 


file =  "./CCDS/hg38/CCDS.GRCh38.p12.txt"

hgi = read.table( file , header=TRUE,sep="\t",stringsAsFactors = FALSE,
                  na.strings=".", quote = "", fill = TRUE, comment.char = '!')

hgi$X.chromosome = paste0("chr",hgi$X.chromosome)

table ( hgi$ccds_status)

hgi = hgi [ !grepl ( "withdraw", hgi$ccds_status,  ignore.case = T), ]

hgi$cds_to = as.numeric ( hgi$cds_to  )
hgi$cds_from = as.numeric (hgi$cds_from  ) 


hgi = hgi [!is.na ( hgi$cds_from) | !is.na ( hgi$cds_to) , ]
unique ( hgi$ccds_status)

hgi$length = hgi$cds_to - hgi$cds_from

# reorder by "Public", "Reviewed, update pending", "Under review, update" and decending length 
# then remove duplicates 

hgi = hgi[ order ( hgi$ccds_status, hgi$gene, -hgi$length),]
hgi = hgi[!duplicated (hgi$gene), ]

# cleanup ranges and split them into two columns
hgi$cds_locations = gsub ( "\\[|\\]| ", "", hgi$cds_locations  )

hgi = hgi %>% 
  mutate(exons = strsplit(as.character(cds_locations), ",")) %>% 
  unnest(exons)

hgi$match_type = NULL 
hgi$cds_locations = NULL 




# split into two columns 
hgi <- hgi %>% tidyr::separate (exons, c("LeftBorder", "RightBorder"), "-")  %>%  data.frame ( stringsAsFactors = F)

hgi$LeftBorder = as.numeric ( hgi$LeftBorder) + 1
hgi$RightBorder = as.numeric ( hgi$RightBorder) + 1



hgi$fusion_5_prime = ifelse ( hgi$cds_strand == "+", 
                             paste0 (hgi$X.chromosome,":", hgi$RightBorder - bp , "-", hgi$RightBorder),
                            paste0 (hgi$X.chromosome,":", hgi$LeftBorder  , "-",  hgi$LeftBorder + bp )
                            
                            )

hgi[hgi$gene=="EWSR1", ] 
hgi[hgi$gene=="PRIM1", ] 
