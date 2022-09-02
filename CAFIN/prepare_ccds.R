library ( dplyr )

hg19 = read.table("./CCDS/hg19/CCDS.GRCh37.p13.txt", header=TRUE,sep="\t",stringsAsFactors = FALSE,
                  na.strings=".", quote = "", fill = TRUE, comment.char = '!')
hg19$X.chromosome = paste0("chr",hg19$X.chromosome)

table ( hg19$ccds_status)

hg19 = hg19 [ !grepl ( "withdraw", hg19$ccds_status,  ignore.case = T), ]

hg19$cds_to = as.numeric ( hg19$cds_to  )
hg19$cds_from = as.numeric (hg19$cds_from  ) 


hg19 = hg19 [!is.na ( hg19$cds_from) | !is.na ( hg19$cds_to) , ]
unique ( hg19$ccds_status)

hg19$length = hg19$cds_to - hg19$cds_from

# reorder by "Public", "Reviewed, update pending", "Under review, update" and decending length 
# then remove duplicates 

hg19 = hg19[ order ( hg19$ccds_status, hg19$gene, -hg19$length),]
hg19 = hg19[!duplicated (hg19$gene), ]

# cleanup ranges and split them into two columns
hg19$cds_locations = gsub ( "\\[|\\]| ", "", hg19$cds_locations  )

hg19 = hg19 %>% 
  mutate(exons = strsplit(as.character(cds_locations), ",")) %>% 
  unnest(exons)

hg19$match_type = NULL 
hg19$cds_locations = NULL 



hg19[hg19$gene=="EWSR1", ] %>% tidyr::separate (exons, c("Left", "Right"), "-")
