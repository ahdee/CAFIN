library ( dplyr )
library ( glue )
library( stringr )
library("Biostrings")
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


# takes a ccds formated file and return the coordinates if exons is a 5' and 3' partners
# two things need to be consider, strand and whether the gene is on the 5' or  3' fusion. 
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



# slurp fastq  into dataframe
fastqToDF <- function ( file, cname=NA ){

  
  
  fastaFile <- readDNAStringSet(  file  )
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df <- data.frame(match = seq_name, cname=sequence)
  
  colnames ( df ) = c( 'match', cname)
  
  return ( df )
  
}

getSeqRange <- function ( use_index, use_genome ) {
  # get sequences 
  # depending on the strand and whether its for 5' or 3' the commands vary. 
  # thus will split into 4 different tables 
  
  # for testing only 
  ## use_index = ihg8
  ## use_genome = hg38 
  
  # split into + and - strands since sequences for negative requires
  # complimentary reverse 
  
  pos = use_genome[use_genome$cds_strand=="+", ]
  neg = use_genome[use_genome$cds_strand=="-", ]
  
  ## there are going to be total of 2 columns ( 5' and 3' )
  ### however we need to run this twice, one for + and then - 
  ### samtools used to pull out sequence of each region 
  
  ### 5' positive strand 
  
  
  posin5 = paste0(out_dir, 'fusion5p_pos', ".IN.tmp")
  posout5 =  paste0(out_dir, 'fusion5p_pos', ".OUT.tmp")
  write( pos$fusion_5_prime , posin5 , append = F)
  cmd = glue ( samtools, " faidx ", use_index, ' -r ', posin5 , " -o ",  posout5 )
  system ( cmd , intern = T )
  
  
  
  ### 5' negative strand
  
  negin5 = paste0(out_dir, 'fusion5p_neg', ".IN.tmp")
  negout5 =  paste0(out_dir, 'fusion5p_neg', ".OUT.tmp")
  write( neg$fusion_5_prime , negin5 , append = F)
  cmd = glue ( samtools, " faidx ", use_index, ' --reverse-complement -r ', negin5 , " -o ",  negout5 )
  system ( cmd , intern = T )
  
  
  ### 3' positive strand
  posin3 = paste0(out_dir, 'fusion3p_pos', ".IN.tmp")
  posout3 =  paste0(out_dir, 'fusion3p_pos', ".OUT.tmp")
  write( pos$fusion_3_prime , posin3 , append = F)
  cmd = glue ( samtools, " faidx ", use_index, ' -r ', posin3 , " -o ",  posout3 )
  system ( cmd , intern = T )
  
  
  ### 3' negative strand
  negin3 = paste0(out_dir, 'fusion3p_neg', ".IN.tmp")
  negout3 =  paste0(out_dir, 'fusion3p_neg', ".OUT.tmp")
  write( neg$fusion_3_prime , negin3 , append = F)
  cmd = glue ( samtools, " faidx ", use_index, ' --reverse-complement -r ', negin3 , " -o ",  negout3 )
  system ( cmd , intern = T )
  
  
  # conver fastq to df and merge 
  
  pos5 = fastqToDF ( posout5 , cname = "5'_seq")
  pos3 = fastqToDF ( posout3 , cname = "3'_seq")
  
  # don't merge 
  pos = cbind ( pos, pos5 )
  p1 = all.equal(pos$fusion_5_prime, pos$match) # this need to be true to preserve order
  pos$match = NULL 
  
  pos = cbind ( pos, pos3 )
  p2 = all.equal(pos$fusion_3_prime, pos$match)
  pos$match = NULL 
  
  ## negative strand 
  neg5 = fastqToDF ( negout5 , cname = "5'_seq")
  neg3 = fastqToDF ( negout3 , cname = "3'_seq")
  # remove /rc trailing
  neg5$match = gsub ( "\\/rc", "", neg5$match )
  neg3$match = gsub ( "\\/rc", "", neg3$match )
  
  neg = cbind ( neg, neg5 )
  p3 = all.equal(neg$fusion_5_prime, neg$match) # this need to be true to preserve order
  neg$match = NULL 
  
  neg = cbind ( neg, neg3 )
  p4 = all.equal(neg$fusion_3_prime, neg$match)
  neg$match = NULL 
  
  if (  all ( p1,p2,p3,p4 ) ) {
    final_df = rbind ( pos, neg )
  } else {
    return ( "ERROR: something went wrong and unable to find sequences")
  }
  
  return ( final_df )
  
  
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

### 
hg38s = getSeqRange( ihg8, hg38)

hg38s[hg38s$gene=="EWSR1", ] 
hg38s[hg38s$gene=="FLI1", ] 

hg38s[hg38s$gene=="PRIM1", ] 
hg38s[hg38s$gene=="NACA", ] 


combodf = hg38s 

can = c("EWSR1-FLI1", "PRIM1-NACA")

for ( fc in can ){
  
  temp = unlist ( stringr::str_split( fc, "-") )
  gene5 = temp[1]
  gene3 = temp[2]
  
  j5 = combodf[combodf$gene == gene5, ]$`5'_seq`
  j3 = combodf[combodf$gene == gene3, ]$`3'_seq`
  
  testj = c()
  
  # By adding the --ignore-case option the command becomes 57x times slower 
  # so will not lowercase 
  
  for ( j5s in j5 ){
    for ( j3s in j3){
      seq = paste0(j5s, "", j3s )
      testj= c( testj, seq)
    }
  }
  
  # for EWSR1-FLI set to: "GCTACGGGCAGCAGAacccttcttatgactca"
  testj[grepl ( 'GCTACGGGCAGCAGAacccttcttatgactca' , testj, ignore.case = T) ]
  
  
  
}



