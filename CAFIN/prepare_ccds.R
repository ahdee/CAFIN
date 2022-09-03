library ( dplyr )
library ( glue )
library( stringr )
bp= 20 # what length you want 

samtools = '/ehome/bin/samtools_rstudio/samtools-1.15.1/bin/bin/samtools'
seqtk = '/ehome/bin/seqtk/seqtk'
ihg8 = '/ehome/resource/genecode/p13/GRCh38.p13.genome.fa'

chg8_pos = glue ( samtools, " faidx ", ihg8, ' range;')
chg8_neg = glue ( samtools, " faidx ", ihg8, ' range|', seqtk, " seq -r;")

# 
out_dir = "./tmp/"
dir.create(out_dir)

# this function will return sequence of range, if strand is negative 
# it will return  a reverse complementary 

getseq<- function(range, strand, build_pos, build_neg){
  
  build = build_pos
  
  if ( strand == "-"){
    build = build_neg
  }
  
  build = gsub ( "range", range, build )
  
  seq = system  (build, intern = TRUE )
  seq = as.character ( stringr::str_split(seq, " ")[2] )
  
  return ( as.character ( seq)  )
}


# takes a ccds formated file and return the coordinates for 5' and 3' partners 
# bp is how many extra base pair to get 
# this will only return coordinates and the actual seq for that use getseq

get_db = function ( file , bp ){
  
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
  
  hgi$fusion_3_prime = ifelse ( hgi$cds_strand == "+", 
                                paste0 (hgi$X.chromosome,":", hgi$LeftBorder , "-",  hgi$LeftBorder + bp ),
                                paste0 (hgi$X.chromosome,":", hgi$RightBorder - bp , "-",   hgi$RightBorder  )
                                
  )
  
  
  return ( hgi )
  
}

file =  "./CCDS/hg38/CCDS.GRCh38.p12.txt"

hg38 = get_db ( file, bp=25)



# see script test_hg38.R for testing. 

# -/-
hg38[hg38$gene=="EWSR1", ] 
hg38[hg38$gene=="FLI1", ] 
# +/+
hg38[hg38$gene=="PRIM1", ] 
hg38[hg38$gene=="NACA", ] 



# get sequences 
# depending on the strand and whether its for 5' or 3' the commands vary. 
# thus will split into 4 different tables 

pos = hg38[hg38$cds_strand=="+", ]
neg = hg38[hg38$cds_strand=="-", ]


# get positive strand 5' 


posin5 = paste0(out_dir, 'fusion5p_pos', ".IN.tmp")
posout5 =  paste0(out_dir, 'fusion5p_pos', ".OUT.tmp")
write( pos$fusion_5_prime , posin5 , append = F)
# ihg8 is the index 
cmd = glue ( samtools, " faidx ", ihg8, ' -r ', posin5 , " -o ",  posout5 )
system ( cmd , intern = T )

# get negative strand 5' 

negin5 = paste0(out_dir, 'fusion5p_neg', ".IN.tmp")
negout5 =  paste0(out_dir, 'fusion5p_neg', ".OUT.tmp")
write( neg$fusion_5_prime , negin5 , append = F)
# ihg8 is the index 
cmd = glue ( samtools, " faidx ", ihg8, ' --reverse-complement -r ', negin5 , " -o ",  negout5 )
system ( cmd , intern = T )

### 3' positive
posin3 = paste0(out_dir, 'fusion3p_pos', ".IN.tmp")
posout3 =  paste0(out_dir, 'fusion3p_pos', ".OUT.tmp")

write( pos$fusion_3_prime , posin3 , append = F)


cmd = glue ( samtools, " faidx ", ihg8, ' -r ', posin3 , " -o ",  posout3 )
system ( cmd , intern = T )


### 3' neg
negin3 = paste0(out_dir, 'fusion3p_neg', ".IN.tmp")
negout3 =  paste0(out_dir, 'fusion3p_neg', ".OUT.tmp")

write( neg$fusion_3_prime , negin3 , append = F)


cmd = glue ( samtools, " faidx ", ihg8, ' --reverse-complement -r ', negin3 , " -o ",  negout3 )
system ( cmd , intern = T )





library("Biostrings")

fastaFile <- readDNAStringSet(  paste0(out_dir, 'fusion5p_pos', ".OUT_tmp") )
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)






hg38$`3'_seq` = as.character ( mapply(getseq, 
                                         range=hg38$fusion_3_prime, 
                                         strand=hg38$cds_strand, 
                                         build_pos =chg8_pos, 
                                         build_neg=chg8_neg )
)











