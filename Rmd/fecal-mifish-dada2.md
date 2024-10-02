NFS_bogo_fecal_mifish_dada2
================
dsb
2024-09-24

Fecal samples metabarcoded with MiFish and sequenced on an Illumina
MiSeq.

Lab work by Jamie Musbach

Primers were trimmed using Cutadapt prior to processing here.

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
# file location
path <- "/genetics/edna/workdir/northernfurseals/bogo_fecal/trimmed/"

path
```

    ## [1] "/genetics/edna/workdir/northernfurseals/bogo_fecal/trimmed/"

``` r
list.files(path)
```

    ##   [1] "e03719-plate1-a_S1_R1.fastq.gz"    "e03719-plate1-a_S1_R2.fastq.gz"   
    ##   [3] "e03719-plate1-b_S49_R1.fastq.gz"   "e03719-plate1-b_S49_R2.fastq.gz"  
    ##   [5] "e03719-plate2-a_S97_R1.fastq.gz"   "e03719-plate2-a_S97_R2.fastq.gz"  
    ##   [7] "e03719-plate2-b_S145_R1.fastq.gz"  "e03719-plate2-b_S145_R2.fastq.gz" 
    ##   [9] "e03719-plate3-a_S193_R1.fastq.gz"  "e03719-plate3-a_S193_R2.fastq.gz" 
    ##  [11] "e03719-plate3-b_S241_R1.fastq.gz"  "e03719-plate3-b_S241_R2.fastq.gz" 
    ##  [13] "e03720-plate1-a_S2_R1.fastq.gz"    "e03720-plate1-a_S2_R2.fastq.gz"   
    ##  [15] "e03720-plate1-b_S50_R1.fastq.gz"   "e03720-plate1-b_S50_R2.fastq.gz"  
    ##  [17] "e03720-plate2-a_S98_R1.fastq.gz"   "e03720-plate2-a_S98_R2.fastq.gz"  
    ##  [19] "e03720-plate2-b_S146_R1.fastq.gz"  "e03720-plate2-b_S146_R2.fastq.gz" 
    ##  [21] "e03720-plate3-a_S194_R1.fastq.gz"  "e03720-plate3-a_S194_R2.fastq.gz" 
    ##  [23] "e03720-plate3-b_S242_R1.fastq.gz"  "e03720-plate3-b_S242_R2.fastq.gz" 
    ##  [25] "e03721-plate1-a_S3_R1.fastq.gz"    "e03721-plate1-a_S3_R2.fastq.gz"   
    ##  [27] "e03721-plate1-b_S51_R1.fastq.gz"   "e03721-plate1-b_S51_R2.fastq.gz"  
    ##  [29] "e03721-plate2-a_S99_R1.fastq.gz"   "e03721-plate2-a_S99_R2.fastq.gz"  
    ##  [31] "e03721-plate2-b_S147_R1.fastq.gz"  "e03721-plate2-b_S147_R2.fastq.gz" 
    ##  [33] "e03721-plate3-a_S195_R1.fastq.gz"  "e03721-plate3-a_S195_R2.fastq.gz" 
    ##  [35] "e03721-plate3-b_S243_R1.fastq.gz"  "e03721-plate3-b_S243_R2.fastq.gz" 
    ##  [37] "e03722-plate1-a_S4_R1.fastq.gz"    "e03722-plate1-a_S4_R2.fastq.gz"   
    ##  [39] "e03722-plate1-b_S52_R1.fastq.gz"   "e03722-plate1-b_S52_R2.fastq.gz"  
    ##  [41] "e03722-plate2-a_S100_R1.fastq.gz"  "e03722-plate2-a_S100_R2.fastq.gz" 
    ##  [43] "e03722-plate2-b_S148_R1.fastq.gz"  "e03722-plate2-b_S148_R2.fastq.gz" 
    ##  [45] "e03722-plate3-a_S196_R1.fastq.gz"  "e03722-plate3-a_S196_R2.fastq.gz" 
    ##  [47] "e03722-plate3-b_S244_R1.fastq.gz"  "e03722-plate3-b_S244_R2.fastq.gz" 
    ##  [49] "e03723-plate1-a_S5_R1.fastq.gz"    "e03723-plate1-a_S5_R2.fastq.gz"   
    ##  [51] "e03723-plate1-b_S53_R1.fastq.gz"   "e03723-plate1-b_S53_R2.fastq.gz"  
    ##  [53] "e03723-plate2-a_S101_R1.fastq.gz"  "e03723-plate2-a_S101_R2.fastq.gz" 
    ##  [55] "e03723-plate2-b_S149_R1.fastq.gz"  "e03723-plate2-b_S149_R2.fastq.gz" 
    ##  [57] "e03723-plate3-a_S197_R1.fastq.gz"  "e03723-plate3-a_S197_R2.fastq.gz" 
    ##  [59] "e03723-plate3-b_S245_R1.fastq.gz"  "e03723-plate3-b_S245_R2.fastq.gz" 
    ##  [61] "e03724-plate1-a_S6_R1.fastq.gz"    "e03724-plate1-a_S6_R2.fastq.gz"   
    ##  [63] "e03724-plate1-b_S54_R1.fastq.gz"   "e03724-plate1-b_S54_R2.fastq.gz"  
    ##  [65] "e03724-plate2-a_S102_R1.fastq.gz"  "e03724-plate2-a_S102_R2.fastq.gz" 
    ##  [67] "e03724-plate2-b_S150_R1.fastq.gz"  "e03724-plate2-b_S150_R2.fastq.gz" 
    ##  [69] "e03724-plate3-a_S198_R1.fastq.gz"  "e03724-plate3-a_S198_R2.fastq.gz" 
    ##  [71] "e03724-plate3-b_S246_R1.fastq.gz"  "e03724-plate3-b_S246_R2.fastq.gz" 
    ##  [73] "e03725-plate1-a_S7_R1.fastq.gz"    "e03725-plate1-a_S7_R2.fastq.gz"   
    ##  [75] "e03725-plate1-b_S55_R1.fastq.gz"   "e03725-plate1-b_S55_R2.fastq.gz"  
    ##  [77] "e03725-plate2-a_S103_R1.fastq.gz"  "e03725-plate2-a_S103_R2.fastq.gz" 
    ##  [79] "e03725-plate2-b_S151_R1.fastq.gz"  "e03725-plate2-b_S151_R2.fastq.gz" 
    ##  [81] "e03725-plate3-a_S199_R1.fastq.gz"  "e03725-plate3-a_S199_R2.fastq.gz" 
    ##  [83] "e03725-plate3-b_S247_R1.fastq.gz"  "e03725-plate3-b_S247_R2.fastq.gz" 
    ##  [85] "e03726-plate1-a_S8_R1.fastq.gz"    "e03726-plate1-a_S8_R2.fastq.gz"   
    ##  [87] "e03726-plate1-b_S56_R1.fastq.gz"   "e03726-plate1-b_S56_R2.fastq.gz"  
    ##  [89] "e03726-plate2-a_S104_R1.fastq.gz"  "e03726-plate2-a_S104_R2.fastq.gz" 
    ##  [91] "e03726-plate2-b_S152_R1.fastq.gz"  "e03726-plate2-b_S152_R2.fastq.gz" 
    ##  [93] "e03726-plate3-a_S200_R1.fastq.gz"  "e03726-plate3-a_S200_R2.fastq.gz" 
    ##  [95] "e03726-plate3-b_S248_R1.fastq.gz"  "e03726-plate3-b_S248_R2.fastq.gz" 
    ##  [97] "e03727-plate1-a_S9_R1.fastq.gz"    "e03727-plate1-a_S9_R2.fastq.gz"   
    ##  [99] "e03727-plate1-b_S57_R1.fastq.gz"   "e03727-plate1-b_S57_R2.fastq.gz"  
    ## [101] "e03727-plate2-a_S105_R1.fastq.gz"  "e03727-plate2-a_S105_R2.fastq.gz" 
    ## [103] "e03727-plate2-b_S153_R1.fastq.gz"  "e03727-plate2-b_S153_R2.fastq.gz" 
    ## [105] "e03727-plate3-a_S201_R1.fastq.gz"  "e03727-plate3-a_S201_R2.fastq.gz" 
    ## [107] "e03727-plate3-b_S249_R1.fastq.gz"  "e03727-plate3-b_S249_R2.fastq.gz" 
    ## [109] "e03728-plate1-a_S10_R1.fastq.gz"   "e03728-plate1-a_S10_R2.fastq.gz"  
    ## [111] "e03728-plate1-b_S59_R1.fastq.gz"   "e03728-plate1-b_S59_R2.fastq.gz"  
    ## [113] "e03728-plate2-a_S106_R1.fastq.gz"  "e03728-plate2-a_S106_R2.fastq.gz" 
    ## [115] "e03728-plate2-b_S155_R1.fastq.gz"  "e03728-plate2-b_S155_R2.fastq.gz" 
    ## [117] "e03728-plate3-a_S202_R1.fastq.gz"  "e03728-plate3-a_S202_R2.fastq.gz" 
    ## [119] "e03728-plate3-b_S251_R1.fastq.gz"  "e03728-plate3-b_S251_R2.fastq.gz" 
    ## [121] "e03729-plate1-a_S11_R1.fastq.gz"   "e03729-plate1-a_S11_R2.fastq.gz"  
    ## [123] "e03729-plate1-b_S60_R1.fastq.gz"   "e03729-plate1-b_S60_R2.fastq.gz"  
    ## [125] "e03729-plate2-a_S107_R1.fastq.gz"  "e03729-plate2-a_S107_R2.fastq.gz" 
    ## [127] "e03729-plate2-b_S156_R1.fastq.gz"  "e03729-plate2-b_S156_R2.fastq.gz" 
    ## [129] "e03729-plate3-a_S203_R1.fastq.gz"  "e03729-plate3-a_S203_R2.fastq.gz" 
    ## [131] "e03729-plate3-b_S252_R1.fastq.gz"  "e03729-plate3-b_S252_R2.fastq.gz" 
    ## [133] "e03730-plate1-a_S13_R1.fastq.gz"   "e03730-plate1-a_S13_R2.fastq.gz"  
    ## [135] "e03730-plate1-b_S61_R1.fastq.gz"   "e03730-plate1-b_S61_R2.fastq.gz"  
    ## [137] "e03730-plate2-a_S109_R1.fastq.gz"  "e03730-plate2-a_S109_R2.fastq.gz" 
    ## [139] "e03730-plate2-b_S157_R1.fastq.gz"  "e03730-plate2-b_S157_R2.fastq.gz" 
    ## [141] "e03730-plate3-a_S205_R1.fastq.gz"  "e03730-plate3-a_S205_R2.fastq.gz" 
    ## [143] "e03730-plate3-b_S253_R1.fastq.gz"  "e03730-plate3-b_S253_R2.fastq.gz" 
    ## [145] "e03731-plate1-a_S14_R1.fastq.gz"   "e03731-plate1-a_S14_R2.fastq.gz"  
    ## [147] "e03731-plate1-b_S62_R1.fastq.gz"   "e03731-plate1-b_S62_R2.fastq.gz"  
    ## [149] "e03731-plate2-a_S110_R1.fastq.gz"  "e03731-plate2-a_S110_R2.fastq.gz" 
    ## [151] "e03731-plate2-b_S158_R1.fastq.gz"  "e03731-plate2-b_S158_R2.fastq.gz" 
    ## [153] "e03731-plate3-a_S206_R1.fastq.gz"  "e03731-plate3-a_S206_R2.fastq.gz" 
    ## [155] "e03731-plate3-b_S254_R1.fastq.gz"  "e03731-plate3-b_S254_R2.fastq.gz" 
    ## [157] "e03732-plate1-a_S15_R1.fastq.gz"   "e03732-plate1-a_S15_R2.fastq.gz"  
    ## [159] "e03732-plate1-b_S63_R1.fastq.gz"   "e03732-plate1-b_S63_R2.fastq.gz"  
    ## [161] "e03732-plate2-a_S111_R1.fastq.gz"  "e03732-plate2-a_S111_R2.fastq.gz" 
    ## [163] "e03732-plate2-b_S159_R1.fastq.gz"  "e03732-plate2-b_S159_R2.fastq.gz" 
    ## [165] "e03732-plate3-a_S207_R1.fastq.gz"  "e03732-plate3-a_S207_R2.fastq.gz" 
    ## [167] "e03732-plate3-b_S255_R1.fastq.gz"  "e03732-plate3-b_S255_R2.fastq.gz" 
    ## [169] "e03733-plate1-a_S16_R1.fastq.gz"   "e03733-plate1-a_S16_R2.fastq.gz"  
    ## [171] "e03733-plate1-b_S64_R1.fastq.gz"   "e03733-plate1-b_S64_R2.fastq.gz"  
    ## [173] "e03733-plate2-a_S112_R1.fastq.gz"  "e03733-plate2-a_S112_R2.fastq.gz" 
    ## [175] "e03733-plate2-b_S160_R1.fastq.gz"  "e03733-plate2-b_S160_R2.fastq.gz" 
    ## [177] "e03733-plate3-a_S208_R1.fastq.gz"  "e03733-plate3-a_S208_R2.fastq.gz" 
    ## [179] "e03733-plate3-b_S256_R1.fastq.gz"  "e03733-plate3-b_S256_R2.fastq.gz" 
    ## [181] "e03734-plate1-a_S17_R1.fastq.gz"   "e03734-plate1-a_S17_R2.fastq.gz"  
    ## [183] "e03734-plate1-b_S65_R1.fastq.gz"   "e03734-plate1-b_S65_R2.fastq.gz"  
    ## [185] "e03734-plate2-a_S113_R1.fastq.gz"  "e03734-plate2-a_S113_R2.fastq.gz" 
    ## [187] "e03734-plate2-b_S161_R1.fastq.gz"  "e03734-plate2-b_S161_R2.fastq.gz" 
    ## [189] "e03734-plate3-a_S209_R1.fastq.gz"  "e03734-plate3-a_S209_R2.fastq.gz" 
    ## [191] "e03734-plate3-b_S257_R1.fastq.gz"  "e03734-plate3-b_S257_R2.fastq.gz" 
    ## [193] "e03735-plate1-a_S18_R1.fastq.gz"   "e03735-plate1-a_S18_R2.fastq.gz"  
    ## [195] "e03735-plate1-b_S66_R1.fastq.gz"   "e03735-plate1-b_S66_R2.fastq.gz"  
    ## [197] "e03735-plate2-a_S114_R1.fastq.gz"  "e03735-plate2-a_S114_R2.fastq.gz" 
    ## [199] "e03735-plate2-b_S162_R1.fastq.gz"  "e03735-plate2-b_S162_R2.fastq.gz" 
    ## [201] "e03735-plate3-a_S210_R1.fastq.gz"  "e03735-plate3-a_S210_R2.fastq.gz" 
    ## [203] "e03735-plate3-b_S258_R1.fastq.gz"  "e03735-plate3-b_S258_R2.fastq.gz" 
    ## [205] "e03736-plate1-a_S19_R1.fastq.gz"   "e03736-plate1-a_S19_R2.fastq.gz"  
    ## [207] "e03736-plate1-b_S67_R1.fastq.gz"   "e03736-plate1-b_S67_R2.fastq.gz"  
    ## [209] "e03736-plate2-a_S115_R1.fastq.gz"  "e03736-plate2-a_S115_R2.fastq.gz" 
    ## [211] "e03736-plate2-b_S163_R1.fastq.gz"  "e03736-plate2-b_S163_R2.fastq.gz" 
    ## [213] "e03736-plate3-a_S211_R1.fastq.gz"  "e03736-plate3-a_S211_R2.fastq.gz" 
    ## [215] "e03736-plate3-b_S259_R1.fastq.gz"  "e03736-plate3-b_S259_R2.fastq.gz" 
    ## [217] "e03737-plate1-a_S20_R1.fastq.gz"   "e03737-plate1-a_S20_R2.fastq.gz"  
    ## [219] "e03737-plate1-b_S68_R1.fastq.gz"   "e03737-plate1-b_S68_R2.fastq.gz"  
    ## [221] "e03737-plate2-a_S116_R1.fastq.gz"  "e03737-plate2-a_S116_R2.fastq.gz" 
    ## [223] "e03737-plate2-b_S164_R1.fastq.gz"  "e03737-plate2-b_S164_R2.fastq.gz" 
    ## [225] "e03737-plate3-a_S212_R1.fastq.gz"  "e03737-plate3-a_S212_R2.fastq.gz" 
    ## [227] "e03737-plate3-b_S260_R1.fastq.gz"  "e03737-plate3-b_S260_R2.fastq.gz" 
    ## [229] "e03738-plate1-a_S21_R1.fastq.gz"   "e03738-plate1-a_S21_R2.fastq.gz"  
    ## [231] "e03738-plate1-b_S69_R1.fastq.gz"   "e03738-plate1-b_S69_R2.fastq.gz"  
    ## [233] "e03738-plate2-a_S117_R1.fastq.gz"  "e03738-plate2-a_S117_R2.fastq.gz" 
    ## [235] "e03738-plate2-b_S165_R1.fastq.gz"  "e03738-plate2-b_S165_R2.fastq.gz" 
    ## [237] "e03738-plate3-a_S213_R1.fastq.gz"  "e03738-plate3-a_S213_R2.fastq.gz" 
    ## [239] "e03738-plate3-b_S261_R1.fastq.gz"  "e03738-plate3-b_S261_R2.fastq.gz" 
    ## [241] "e03739-plate1-a_S22_R1.fastq.gz"   "e03739-plate1-a_S22_R2.fastq.gz"  
    ## [243] "e03739-plate1-b_S71_R1.fastq.gz"   "e03739-plate1-b_S71_R2.fastq.gz"  
    ## [245] "e03739-plate2-a_S118_R1.fastq.gz"  "e03739-plate2-a_S118_R2.fastq.gz" 
    ## [247] "e03739-plate2-b_S167_R1.fastq.gz"  "e03739-plate2-b_S167_R2.fastq.gz" 
    ## [249] "e03739-plate3-a_S214_R1.fastq.gz"  "e03739-plate3-a_S214_R2.fastq.gz" 
    ## [251] "e03739-plate3-b_S263_R1.fastq.gz"  "e03739-plate3-b_S263_R2.fastq.gz" 
    ## [253] "e03740-plate1-a_S23_R1.fastq.gz"   "e03740-plate1-a_S23_R2.fastq.gz"  
    ## [255] "e03740-plate1-b_S72_R1.fastq.gz"   "e03740-plate1-b_S72_R2.fastq.gz"  
    ## [257] "e03740-plate2-a_S119_R1.fastq.gz"  "e03740-plate2-a_S119_R2.fastq.gz" 
    ## [259] "e03740-plate2-b_S168_R1.fastq.gz"  "e03740-plate2-b_S168_R2.fastq.gz" 
    ## [261] "e03740-plate3-a_S215_R1.fastq.gz"  "e03740-plate3-a_S215_R2.fastq.gz" 
    ## [263] "e03740-plate3-b_S264_R1.fastq.gz"  "e03740-plate3-b_S264_R2.fastq.gz" 
    ## [265] "e03741-plate1-a_S25_R1.fastq.gz"   "e03741-plate1-a_S25_R2.fastq.gz"  
    ## [267] "e03741-plate1-b_S73_R1.fastq.gz"   "e03741-plate1-b_S73_R2.fastq.gz"  
    ## [269] "e03741-plate2-a_S121_R1.fastq.gz"  "e03741-plate2-a_S121_R2.fastq.gz" 
    ## [271] "e03741-plate2-b_S169_R1.fastq.gz"  "e03741-plate2-b_S169_R2.fastq.gz" 
    ## [273] "e03741-plate3-a_S217_R1.fastq.gz"  "e03741-plate3-a_S217_R2.fastq.gz" 
    ## [275] "e03741-plate3-b_S265_R1.fastq.gz"  "e03741-plate3-b_S265_R2.fastq.gz" 
    ## [277] "e03742-plate1-a_S26_R1.fastq.gz"   "e03742-plate1-a_S26_R2.fastq.gz"  
    ## [279] "e03742-plate1-b_S74_R1.fastq.gz"   "e03742-plate1-b_S74_R2.fastq.gz"  
    ## [281] "e03742-plate2-a_S122_R1.fastq.gz"  "e03742-plate2-a_S122_R2.fastq.gz" 
    ## [283] "e03742-plate2-b_S170_R1.fastq.gz"  "e03742-plate2-b_S170_R2.fastq.gz" 
    ## [285] "e03742-plate3-a_S218_R1.fastq.gz"  "e03742-plate3-a_S218_R2.fastq.gz" 
    ## [287] "e03742-plate3-b_S266_R1.fastq.gz"  "e03742-plate3-b_S266_R2.fastq.gz" 
    ## [289] "e03743-plate1-a_S27_R1.fastq.gz"   "e03743-plate1-a_S27_R2.fastq.gz"  
    ## [291] "e03743-plate1-b_S75_R1.fastq.gz"   "e03743-plate1-b_S75_R2.fastq.gz"  
    ## [293] "e03743-plate2-a_S123_R1.fastq.gz"  "e03743-plate2-a_S123_R2.fastq.gz" 
    ## [295] "e03743-plate2-b_S171_R1.fastq.gz"  "e03743-plate2-b_S171_R2.fastq.gz" 
    ## [297] "e03743-plate3-a_S219_R1.fastq.gz"  "e03743-plate3-a_S219_R2.fastq.gz" 
    ## [299] "e03743-plate3-b_S267_R1.fastq.gz"  "e03743-plate3-b_S267_R2.fastq.gz" 
    ## [301] "e03744-plate1-a_S28_R1.fastq.gz"   "e03744-plate1-a_S28_R2.fastq.gz"  
    ## [303] "e03744-plate1-b_S76_R1.fastq.gz"   "e03744-plate1-b_S76_R2.fastq.gz"  
    ## [305] "e03744-plate2-a_S124_R1.fastq.gz"  "e03744-plate2-a_S124_R2.fastq.gz" 
    ## [307] "e03744-plate2-b_S172_R1.fastq.gz"  "e03744-plate2-b_S172_R2.fastq.gz" 
    ## [309] "e03744-plate3-a_S220_R1.fastq.gz"  "e03744-plate3-a_S220_R2.fastq.gz" 
    ## [311] "e03744-plate3-b_S268_R1.fastq.gz"  "e03744-plate3-b_S268_R2.fastq.gz" 
    ## [313] "e03745-plate1-a_S29_R1.fastq.gz"   "e03745-plate1-a_S29_R2.fastq.gz"  
    ## [315] "e03745-plate1-b_S78_R1.fastq.gz"   "e03745-plate1-b_S78_R2.fastq.gz"  
    ## [317] "e03745-plate2-a_S125_R1.fastq.gz"  "e03745-plate2-a_S125_R2.fastq.gz" 
    ## [319] "e03745-plate2-b_S174_R1.fastq.gz"  "e03745-plate2-b_S174_R2.fastq.gz" 
    ## [321] "e03745-plate3-a_S221_R1.fastq.gz"  "e03745-plate3-a_S221_R2.fastq.gz" 
    ## [323] "e03745-plate3-b_S270_R1.fastq.gz"  "e03745-plate3-b_S270_R2.fastq.gz" 
    ## [325] "e03746-plate1-a_S30_R1.fastq.gz"   "e03746-plate1-a_S30_R2.fastq.gz"  
    ## [327] "e03746-plate1-b_S77_R1.fastq.gz"   "e03746-plate1-b_S77_R2.fastq.gz"  
    ## [329] "e03746-plate2-a_S126_R1.fastq.gz"  "e03746-plate2-a_S126_R2.fastq.gz" 
    ## [331] "e03746-plate2-b_S173_R1.fastq.gz"  "e03746-plate2-b_S173_R2.fastq.gz" 
    ## [333] "e03746-plate3-a_S222_R1.fastq.gz"  "e03746-plate3-a_S222_R2.fastq.gz" 
    ## [335] "e03746-plate3-b_S269_R1.fastq.gz"  "e03746-plate3-b_S269_R2.fastq.gz" 
    ## [337] "e03747-plate1-a_S31_R1.fastq.gz"   "e03747-plate1-a_S31_R2.fastq.gz"  
    ## [339] "e03747-plate1-b_S79_R1.fastq.gz"   "e03747-plate1-b_S79_R2.fastq.gz"  
    ## [341] "e03747-plate2-a_S127_R1.fastq.gz"  "e03747-plate2-a_S127_R2.fastq.gz" 
    ## [343] "e03747-plate2-b_S175_R1.fastq.gz"  "e03747-plate2-b_S175_R2.fastq.gz" 
    ## [345] "e03747-plate3-a_S223_R1.fastq.gz"  "e03747-plate3-a_S223_R2.fastq.gz" 
    ## [347] "e03747-plate3-b_S271_R1.fastq.gz"  "e03747-plate3-b_S271_R2.fastq.gz" 
    ## [349] "e03748-plate1-a_S32_R1.fastq.gz"   "e03748-plate1-a_S32_R2.fastq.gz"  
    ## [351] "e03748-plate1-b_S80_R1.fastq.gz"   "e03748-plate1-b_S80_R2.fastq.gz"  
    ## [353] "e03748-plate2-a_S128_R1.fastq.gz"  "e03748-plate2-a_S128_R2.fastq.gz" 
    ## [355] "e03748-plate2-b_S176_R1.fastq.gz"  "e03748-plate2-b_S176_R2.fastq.gz" 
    ## [357] "e03748-plate3-a_S224_R1.fastq.gz"  "e03748-plate3-a_S224_R2.fastq.gz" 
    ## [359] "e03748-plate3-b_S272_R1.fastq.gz"  "e03748-plate3-b_S272_R2.fastq.gz" 
    ## [361] "e03749-plate1-a_S34_R1.fastq.gz"   "e03749-plate1-a_S34_R2.fastq.gz"  
    ## [363] "e03749-plate1-b_S81_R1.fastq.gz"   "e03749-plate1-b_S81_R2.fastq.gz"  
    ## [365] "e03749-plate2-a_S130_R1.fastq.gz"  "e03749-plate2-a_S130_R2.fastq.gz" 
    ## [367] "e03749-plate2-b_S177_R1.fastq.gz"  "e03749-plate2-b_S177_R2.fastq.gz" 
    ## [369] "e03749-plate3-a_S226_R1.fastq.gz"  "e03749-plate3-a_S226_R2.fastq.gz" 
    ## [371] "e03749-plate3-b_S273_R1.fastq.gz"  "e03749-plate3-b_S273_R2.fastq.gz" 
    ## [373] "e03750-plate1-a_S35_R1.fastq.gz"   "e03750-plate1-a_S35_R2.fastq.gz"  
    ## [375] "e03750-plate1-b_S82_R1.fastq.gz"   "e03750-plate1-b_S82_R2.fastq.gz"  
    ## [377] "e03750-plate2-a_S131_R1.fastq.gz"  "e03750-plate2-a_S131_R2.fastq.gz" 
    ## [379] "e03750-plate2-b_S178_R1.fastq.gz"  "e03750-plate2-b_S178_R2.fastq.gz" 
    ## [381] "e03750-plate3-a_S227_R1.fastq.gz"  "e03750-plate3-a_S227_R2.fastq.gz" 
    ## [383] "e03750-plate3-b_S274_R1.fastq.gz"  "e03750-plate3-b_S274_R2.fastq.gz" 
    ## [385] "e03751-plate1-a_S36_R1.fastq.gz"   "e03751-plate1-a_S36_R2.fastq.gz"  
    ## [387] "e03751-plate1-b_S83_R1.fastq.gz"   "e03751-plate1-b_S83_R2.fastq.gz"  
    ## [389] "e03751-plate2-a_S132_R1.fastq.gz"  "e03751-plate2-a_S132_R2.fastq.gz" 
    ## [391] "e03751-plate2-b_S179_R1.fastq.gz"  "e03751-plate2-b_S179_R2.fastq.gz" 
    ## [393] "e03751-plate3-a_S228_R1.fastq.gz"  "e03751-plate3-a_S228_R2.fastq.gz" 
    ## [395] "e03751-plate3-b_S275_R1.fastq.gz"  "e03751-plate3-b_S275_R2.fastq.gz" 
    ## [397] "e03752-plate1-a_S37_R1.fastq.gz"   "e03752-plate1-a_S37_R2.fastq.gz"  
    ## [399] "e03752-plate1-b_S84_R1.fastq.gz"   "e03752-plate1-b_S84_R2.fastq.gz"  
    ## [401] "e03752-plate2-a_S133_R1.fastq.gz"  "e03752-plate2-a_S133_R2.fastq.gz" 
    ## [403] "e03752-plate2-b_S180_R1.fastq.gz"  "e03752-plate2-b_S180_R2.fastq.gz" 
    ## [405] "e03752-plate3-a_S229_R1.fastq.gz"  "e03752-plate3-a_S229_R2.fastq.gz" 
    ## [407] "e03752-plate3-b_S276_R1.fastq.gz"  "e03752-plate3-b_S276_R2.fastq.gz" 
    ## [409] "e03753-plate1-a_S38_R1.fastq.gz"   "e03753-plate1-a_S38_R2.fastq.gz"  
    ## [411] "e03753-plate1-b_S85_R1.fastq.gz"   "e03753-plate1-b_S85_R2.fastq.gz"  
    ## [413] "e03753-plate2-a_S134_R1.fastq.gz"  "e03753-plate2-a_S134_R2.fastq.gz" 
    ## [415] "e03753-plate2-b_S181_R1.fastq.gz"  "e03753-plate2-b_S181_R2.fastq.gz" 
    ## [417] "e03753-plate3-a_S230_R1.fastq.gz"  "e03753-plate3-a_S230_R2.fastq.gz" 
    ## [419] "e03753-plate3-b_S277_R1.fastq.gz"  "e03753-plate3-b_S277_R2.fastq.gz" 
    ## [421] "e03754-plate1-a_S39_R1.fastq.gz"   "e03754-plate1-a_S39_R2.fastq.gz"  
    ## [423] "e03754-plate1-b_S86_R1.fastq.gz"   "e03754-plate1-b_S86_R2.fastq.gz"  
    ## [425] "e03754-plate2-a_S135_R1.fastq.gz"  "e03754-plate2-a_S135_R2.fastq.gz" 
    ## [427] "e03754-plate2-b_S182_R1.fastq.gz"  "e03754-plate2-b_S182_R2.fastq.gz" 
    ## [429] "e03754-plate3-a_S231_R1.fastq.gz"  "e03754-plate3-a_S231_R2.fastq.gz" 
    ## [431] "e03754-plate3-b_S278_R1.fastq.gz"  "e03754-plate3-b_S278_R2.fastq.gz" 
    ## [433] "e03755-plate1-a_S40_R1.fastq.gz"   "e03755-plate1-a_S40_R2.fastq.gz"  
    ## [435] "e03755-plate1-b_S87_R1.fastq.gz"   "e03755-plate1-b_S87_R2.fastq.gz"  
    ## [437] "e03755-plate2-a_S136_R1.fastq.gz"  "e03755-plate2-a_S136_R2.fastq.gz" 
    ## [439] "e03755-plate2-b_S183_R1.fastq.gz"  "e03755-plate2-b_S183_R2.fastq.gz" 
    ## [441] "e03755-plate3-a_S232_R1.fastq.gz"  "e03755-plate3-a_S232_R2.fastq.gz" 
    ## [443] "e03755-plate3-b_S279_R1.fastq.gz"  "e03755-plate3-b_S279_R2.fastq.gz" 
    ## [445] "e03756-plate1-a_S41_R1.fastq.gz"   "e03756-plate1-a_S41_R2.fastq.gz"  
    ## [447] "e03756-plate1-b_S90_R1.fastq.gz"   "e03756-plate1-b_S90_R2.fastq.gz"  
    ## [449] "e03756-plate2-a_S137_R1.fastq.gz"  "e03756-plate2-a_S137_R2.fastq.gz" 
    ## [451] "e03756-plate2-b_S186_R1.fastq.gz"  "e03756-plate2-b_S186_R2.fastq.gz" 
    ## [453] "e03756-plate3-a_S233_R1.fastq.gz"  "e03756-plate3-a_S233_R2.fastq.gz" 
    ## [455] "e03756-plate3-b_S282_R1.fastq.gz"  "e03756-plate3-b_S282_R2.fastq.gz" 
    ## [457] "e03757-plate1-a_S42_R1.fastq.gz"   "e03757-plate1-a_S42_R2.fastq.gz"  
    ## [459] "e03757-plate1-b_S91_R1.fastq.gz"   "e03757-plate1-b_S91_R2.fastq.gz"  
    ## [461] "e03757-plate2-a_S138_R1.fastq.gz"  "e03757-plate2-a_S138_R2.fastq.gz" 
    ## [463] "e03757-plate2-b_S187_R1.fastq.gz"  "e03757-plate2-b_S187_R2.fastq.gz" 
    ## [465] "e03757-plate3-a_S234_R1.fastq.gz"  "e03757-plate3-a_S234_R2.fastq.gz" 
    ## [467] "e03757-plate3-b_S283_R1.fastq.gz"  "e03757-plate3-b_S283_R2.fastq.gz" 
    ## [469] "e03758-plate1-a_S43_R1.fastq.gz"   "e03758-plate1-a_S43_R2.fastq.gz"  
    ## [471] "e03758-plate1-b_S92_R1.fastq.gz"   "e03758-plate1-b_S92_R2.fastq.gz"  
    ## [473] "e03758-plate2-a_S139_R1.fastq.gz"  "e03758-plate2-a_S139_R2.fastq.gz" 
    ## [475] "e03758-plate2-b_S188_R1.fastq.gz"  "e03758-plate2-b_S188_R2.fastq.gz" 
    ## [477] "e03758-plate3-a_S235_R1.fastq.gz"  "e03758-plate3-a_S235_R2.fastq.gz" 
    ## [479] "e03758-plate3-b_S284_R1.fastq.gz"  "e03758-plate3-b_S284_R2.fastq.gz" 
    ## [481] "e03759-plate1-a_S44_R1.fastq.gz"   "e03759-plate1-a_S44_R2.fastq.gz"  
    ## [483] "e03759-plate1-b_S93_R1.fastq.gz"   "e03759-plate1-b_S93_R2.fastq.gz"  
    ## [485] "e03759-plate2-a_S140_R1.fastq.gz"  "e03759-plate2-a_S140_R2.fastq.gz" 
    ## [487] "e03759-plate2-b_S189_R1.fastq.gz"  "e03759-plate2-b_S189_R2.fastq.gz" 
    ## [489] "e03759-plate3-a_S236_R1.fastq.gz"  "e03759-plate3-a_S236_R2.fastq.gz" 
    ## [491] "e03759-plate3-b_S285_R1.fastq.gz"  "e03759-plate3-b_S285_R2.fastq.gz" 
    ## [493] "e03760-plate1-a_S45_R1.fastq.gz"   "e03760-plate1-a_S45_R2.fastq.gz"  
    ## [495] "e03760-plate1-b_S94_R1.fastq.gz"   "e03760-plate1-b_S94_R2.fastq.gz"  
    ## [497] "e03760-plate2-a_S141_R1.fastq.gz"  "e03760-plate2-a_S141_R2.fastq.gz" 
    ## [499] "e03760-plate2-b_S190_R1.fastq.gz"  "e03760-plate2-b_S190_R2.fastq.gz" 
    ## [501] "e03760-plate3-a_S237_R1.fastq.gz"  "e03760-plate3-a_S237_R2.fastq.gz" 
    ## [503] "e03760-plate3-b_S286_R1.fastq.gz"  "e03760-plate3-b_S286_R2.fastq.gz" 
    ## [505] "e03761-plate1-a_S46_R1.fastq.gz"   "e03761-plate1-a_S46_R2.fastq.gz"  
    ## [507] "e03761-plate1-b_S95_R1.fastq.gz"   "e03761-plate1-b_S95_R2.fastq.gz"  
    ## [509] "e03761-plate2-a_S142_R1.fastq.gz"  "e03761-plate2-a_S142_R2.fastq.gz" 
    ## [511] "e03761-plate2-b_S191_R1.fastq.gz"  "e03761-plate2-b_S191_R2.fastq.gz" 
    ## [513] "e03761-plate3-a_S238_R1.fastq.gz"  "e03761-plate3-a_S238_R2.fastq.gz" 
    ## [515] "e03761-plate3-b_S287_R1.fastq.gz"  "e03761-plate3-b_S287_R2.fastq.gz" 
    ## [517] "e03762-plate1-EB_S12_R1.fastq.gz"  "e03762-plate1-EB_S12_R2.fastq.gz" 
    ## [519] "e03762-plate2-EB_S108_R1.fastq.gz" "e03762-plate2-EB_S108_R2.fastq.gz"
    ## [521] "e03762-plate3-EB_S204_R1.fastq.gz" "e03762-plate3-EB_S204_R2.fastq.gz"
    ## [523] "e03763-plate1-EB_S24_R1.fastq.gz"  "e03763-plate1-EB_S24_R2.fastq.gz" 
    ## [525] "e03763-plate2-EB_S120_R1.fastq.gz" "e03763-plate2-EB_S120_R2.fastq.gz"
    ## [527] "e03763-plate3-EB_S216_R1.fastq.gz" "e03763-plate3-EB_S216_R2.fastq.gz"
    ## [529] "e03764-plate1-EB_S33_R1.fastq.gz"  "e03764-plate1-EB_S33_R2.fastq.gz" 
    ## [531] "e03764-plate2-EB_S129_R1.fastq.gz" "e03764-plate2-EB_S129_R2.fastq.gz"
    ## [533] "e03764-plate3-EB_S225_R1.fastq.gz" "e03764-plate3-EB_S225_R2.fastq.gz"
    ## [535] "e03765-plate1-EB_S47_R1.fastq.gz"  "e03765-plate1-EB_S47_R2.fastq.gz" 
    ## [537] "e03765-plate2-EB_S143_R1.fastq.gz" "e03765-plate2-EB_S143_R2.fastq.gz"
    ## [539] "e03765-plate3-EB_S239_R1.fastq.gz" "e03765-plate3-EB_S239_R2.fastq.gz"
    ## [541] "e03766-plate1-EB_S58_R1.fastq.gz"  "e03766-plate1-EB_S58_R2.fastq.gz" 
    ## [543] "e03766-plate2-EB_S154_R1.fastq.gz" "e03766-plate2-EB_S154_R2.fastq.gz"
    ## [545] "e03766-plate3-EB_S250_R1.fastq.gz" "e03766-plate3-EB_S250_R2.fastq.gz"
    ## [547] "e03767-plate1-EB_S70_R1.fastq.gz"  "e03767-plate1-EB_S70_R2.fastq.gz" 
    ## [549] "e03767-plate2-EB_S166_R1.fastq.gz" "e03767-plate2-EB_S166_R2.fastq.gz"
    ## [551] "e03767-plate3-EB_S262_R1.fastq.gz" "e03767-plate3-EB_S262_R2.fastq.gz"
    ## [553] "e03768-plate1-EB_S88_R1.fastq.gz"  "e03768-plate1-EB_S88_R2.fastq.gz" 
    ## [555] "e03768-plate2-EB_S184_R1.fastq.gz" "e03768-plate2-EB_S184_R2.fastq.gz"
    ## [557] "e03768-plate3-EB_S280_R1.fastq.gz" "e03768-plate3-EB_S280_R2.fastq.gz"
    ## [559] "e03769-plate1-FB_S89_R1.fastq.gz"  "e03769-plate1-FB_S89_R2.fastq.gz" 
    ## [561] "e03769-plate2-FB_S185_R1.fastq.gz" "e03769-plate2-FB_S185_R2.fastq.gz"
    ## [563] "e03769-plate3-FB_S281_R1.fastq.gz" "e03769-plate3-FB_S281_R2.fastq.gz"
    ## [565] "filtered"                          "NC-plate1_S48_R1.fastq.gz"        
    ## [567] "NC-plate1_S48_R2.fastq.gz"         "NC-plate2_S144_R1.fastq.gz"       
    ## [569] "NC-plate2_S144_R2.fastq.gz"        "NC-plate3_S240_R1.fastq.gz"       
    ## [571] "NC-plate3_S240_R2.fastq.gz"        "PC-plate1_S96_R1.fastq.gz"        
    ## [573] "PC-plate1_S96_R2.fastq.gz"         "PC-plate2_S192_R1.fastq.gz"       
    ## [575] "PC-plate2_S192_R2.fastq.gz"        "PC-plate3_S288_R1.fastq.gz"       
    ## [577] "PC-plate3_S288_R2.fastq.gz"

``` r
fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

``` r
plotQualityProfile(fnFs[1:2])
```

    ## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
    ## of ggplot2 3.3.4.
    ## ℹ The deprecated feature was likely used in the dada2 package.
    ##   Please report the issue at <https://github.com/benjjneb/dada2/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](fecal-mifish-dada2_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:2])
```

![](fecal-mifish-dada2_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

Looks like high quality sequencing data - which coincides with the MiSeq
specs as well.

``` r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

These parameters typically need to be tweaked substantially to deal with
data quality (and amplicon length). In this case, I’ve tested relaxing
the truncLen parameter to allow for more reads to merge and now am
increasing the maxEE and removing the truncLen to see if that helps
retain reads.

- removing the truncLen parameter decreased the number of reads
  retained.
- Kim suggested using 100 for both the fwd and rev reads.

``` r
# modify these parameters - sometimes substantially
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(100,100),
              maxN=0, maxEE=c(2,2), truncQ=c(2,2), rm.phix=TRUE,
              compress=TRUE, multithread=TRUE) 
head(out)
```

    ##                                  reads.in reads.out
    ## e03719-plate1-a_S1_R1.fastq.gz      39155     39017
    ## e03719-plate1-b_S49_R1.fastq.gz     37077     36937
    ## e03719-plate2-a_S97_R1.fastq.gz     39338     38973
    ## e03719-plate2-b_S145_R1.fastq.gz    47160     46766
    ## e03719-plate3-a_S193_R1.fastq.gz    26731     26614
    ## e03719-plate3-b_S241_R1.fastq.gz    32975     32782

Although it looks like these filters keep a fair number of reads, the
real test is with the merger step.

``` r
# ensure there are no samples without reads
# as these don't work in the next steps
out %>%
  as.data.frame() %>%
  filter(reads.out <1)
```

    ## [1] reads.in  reads.out
    ## <0 rows> (or 0-length row.names)

``` r
# filter the matrix to retain samples with >0 reads
```

Let’s move fwd with these for now… and come back if there are other
issues.

### Error rates

``` r
# fwd error rates
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 100780000 total bases in 1007800 reads from 25 samples will be used for learning the error rates.

``` r
# reverse error rates
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 100780000 total bases in 1007800 reads from 25 samples will be used for learning the error rates.

``` r
# plot the erors
p1 <- plotErrors(errF, nominalQ=TRUE)
p2 <- plotErrors(errR, nominalQ=TRUE)

p1
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](fecal-mifish-dada2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
p2
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    ## Transformation introduced infinite values in continuous y-axis

![](fecal-mifish-dada2_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

### Sample inference

``` r
# I tend to use the pseudo pool option 
# to retain rare ASVs across samples

# forwards
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 39017 reads in 7508 unique sequences.
    ## Sample 2 - 36937 reads in 6078 unique sequences.
    ## Sample 3 - 38973 reads in 10757 unique sequences.
    ## Sample 4 - 46766 reads in 11939 unique sequences.
    ## Sample 5 - 26614 reads in 6159 unique sequences.
    ## Sample 6 - 32782 reads in 6641 unique sequences.
    ## Sample 7 - 40836 reads in 6078 unique sequences.
    ## Sample 8 - 39476 reads in 5898 unique sequences.
    ## Sample 9 - 44111 reads in 12168 unique sequences.
    ## Sample 10 - 44482 reads in 10093 unique sequences.
    ## Sample 11 - 25363 reads in 5332 unique sequences.
    ## Sample 12 - 29479 reads in 4836 unique sequences.
    ## Sample 13 - 39151 reads in 7883 unique sequences.
    ## Sample 14 - 53206 reads in 10444 unique sequences.
    ## Sample 15 - 41422 reads in 13768 unique sequences.
    ## Sample 16 - 61733 reads in 17334 unique sequences.
    ## Sample 17 - 25791 reads in 5298 unique sequences.
    ## Sample 18 - 40901 reads in 9257 unique sequences.
    ## Sample 19 - 43354 reads in 6653 unique sequences.
    ## Sample 20 - 52145 reads in 7843 unique sequences.
    ## Sample 21 - 45145 reads in 7400 unique sequences.
    ## Sample 22 - 57210 reads in 8796 unique sequences.
    ## Sample 23 - 28576 reads in 4916 unique sequences.
    ## Sample 24 - 35955 reads in 5804 unique sequences.
    ## Sample 25 - 38375 reads in 10089 unique sequences.
    ## Sample 26 - 49693 reads in 12988 unique sequences.
    ## Sample 27 - 43738 reads in 11270 unique sequences.
    ## Sample 28 - 58206 reads in 15091 unique sequences.
    ## Sample 29 - 29354 reads in 7679 unique sequences.
    ## Sample 30 - 37192 reads in 9709 unique sequences.
    ## Sample 31 - 38644 reads in 7242 unique sequences.
    ## Sample 32 - 41439 reads in 8150 unique sequences.
    ## Sample 33 - 41280 reads in 8537 unique sequences.
    ## Sample 34 - 48845 reads in 11618 unique sequences.
    ## Sample 35 - 21675 reads in 4345 unique sequences.
    ## Sample 36 - 30671 reads in 5902 unique sequences.
    ## Sample 37 - 43037 reads in 8136 unique sequences.
    ## Sample 38 - 42905 reads in 7366 unique sequences.
    ## Sample 39 - 45206 reads in 9151 unique sequences.
    ## Sample 40 - 45682 reads in 7856 unique sequences.
    ## Sample 41 - 25655 reads in 4999 unique sequences.
    ## Sample 42 - 28199 reads in 4974 unique sequences.
    ## Sample 43 - 40184 reads in 9019 unique sequences.
    ## Sample 44 - 43380 reads in 9158 unique sequences.
    ## Sample 45 - 41199 reads in 15815 unique sequences.
    ## Sample 46 - 45616 reads in 16377 unique sequences.
    ## Sample 47 - 23246 reads in 6048 unique sequences.
    ## Sample 48 - 33058 reads in 7543 unique sequences.
    ## Sample 49 - 38693 reads in 5412 unique sequences.
    ## Sample 50 - 38741 reads in 5986 unique sequences.
    ## Sample 51 - 42978 reads in 6353 unique sequences.
    ## Sample 52 - 54868 reads in 8614 unique sequences.
    ## Sample 53 - 27696 reads in 4392 unique sequences.
    ## Sample 54 - 31577 reads in 5255 unique sequences.
    ## Sample 55 - 44388 reads in 8899 unique sequences.
    ## Sample 56 - 41192 reads in 7667 unique sequences.
    ## Sample 57 - 50886 reads in 13290 unique sequences.
    ## Sample 58 - 45257 reads in 11370 unique sequences.
    ## Sample 59 - 29400 reads in 6540 unique sequences.
    ## Sample 60 - 28052 reads in 5737 unique sequences.
    ## Sample 61 - 52070 reads in 7857 unique sequences.
    ## Sample 62 - 42583 reads in 8502 unique sequences.
    ## Sample 63 - 61516 reads in 10724 unique sequences.
    ## Sample 64 - 48410 reads in 11327 unique sequences.
    ## Sample 65 - 37551 reads in 6052 unique sequences.
    ## Sample 66 - 25714 reads in 5718 unique sequences.
    ## Sample 67 - 60579 reads in 6885 unique sequences.
    ## Sample 68 - 46753 reads in 5519 unique sequences.
    ## Sample 69 - 59452 reads in 11371 unique sequences.
    ## Sample 70 - 49504 reads in 8134 unique sequences.
    ## Sample 71 - 32504 reads in 4298 unique sequences.
    ## Sample 72 - 31783 reads in 3875 unique sequences.
    ## Sample 73 - 53713 reads in 8271 unique sequences.
    ## Sample 74 - 53867 reads in 6243 unique sequences.
    ## Sample 75 - 58977 reads in 9435 unique sequences.
    ## Sample 76 - 62290 reads in 7438 unique sequences.
    ## Sample 77 - 33817 reads in 5540 unique sequences.
    ## Sample 78 - 39820 reads in 4870 unique sequences.
    ## Sample 79 - 50850 reads in 7115 unique sequences.
    ## Sample 80 - 56367 reads in 9453 unique sequences.
    ## Sample 81 - 52403 reads in 12040 unique sequences.
    ## Sample 82 - 58087 reads in 16108 unique sequences.
    ## Sample 83 - 33499 reads in 5333 unique sequences.
    ## Sample 84 - 35017 reads in 6644 unique sequences.
    ## Sample 85 - 23523 reads in 4795 unique sequences.
    ## Sample 86 - 41036 reads in 5096 unique sequences.
    ## Sample 87 - 21703 reads in 5629 unique sequences.
    ## Sample 88 - 38492 reads in 6618 unique sequences.
    ## Sample 89 - 11867 reads in 2791 unique sequences.
    ## Sample 90 - 24073 reads in 3464 unique sequences.
    ## Sample 91 - 35678 reads in 4996 unique sequences.
    ## Sample 92 - 40378 reads in 5664 unique sequences.
    ## Sample 93 - 37525 reads in 7168 unique sequences.
    ## Sample 94 - 51034 reads in 8610 unique sequences.
    ## Sample 95 - 23132 reads in 3719 unique sequences.
    ## Sample 96 - 32169 reads in 4861 unique sequences.
    ## Sample 97 - 44286 reads in 5211 unique sequences.
    ## Sample 98 - 42665 reads in 5094 unique sequences.
    ## Sample 99 - 45018 reads in 7358 unique sequences.
    ## Sample 100 - 52002 reads in 9943 unique sequences.
    ## Sample 101 - 27881 reads in 3802 unique sequences.
    ## Sample 102 - 34061 reads in 4127 unique sequences.
    ## Sample 103 - 41451 reads in 5050 unique sequences.
    ## Sample 104 - 37127 reads in 6960 unique sequences.
    ## Sample 105 - 40735 reads in 8085 unique sequences.
    ## Sample 106 - 40887 reads in 10387 unique sequences.
    ## Sample 107 - 24513 reads in 3532 unique sequences.
    ## Sample 108 - 27344 reads in 5614 unique sequences.
    ## Sample 109 - 44606 reads in 10436 unique sequences.
    ## Sample 110 - 42728 reads in 9880 unique sequences.
    ## Sample 111 - 49920 reads in 14139 unique sequences.
    ## Sample 112 - 49529 reads in 14424 unique sequences.
    ## Sample 113 - 28161 reads in 7221 unique sequences.
    ## Sample 114 - 30352 reads in 7638 unique sequences.
    ## Sample 115 - 59561 reads in 10030 unique sequences.
    ## Sample 116 - 44247 reads in 7014 unique sequences.
    ## Sample 117 - 55172 reads in 15516 unique sequences.
    ## Sample 118 - 43140 reads in 12985 unique sequences.
    ## Sample 119 - 32860 reads in 6897 unique sequences.
    ## Sample 120 - 29931 reads in 6795 unique sequences.
    ## Sample 121 - 36762 reads in 12023 unique sequences.
    ## Sample 122 - 36650 reads in 8098 unique sequences.
    ## Sample 123 - 40851 reads in 12729 unique sequences.
    ## Sample 124 - 43349 reads in 9786 unique sequences.
    ## Sample 125 - 24801 reads in 7858 unique sequences.
    ## Sample 126 - 27020 reads in 5804 unique sequences.
    ## Sample 127 - 53712 reads in 13018 unique sequences.
    ## Sample 128 - 40677 reads in 8538 unique sequences.
    ## Sample 129 - 54761 reads in 14220 unique sequences.
    ## Sample 130 - 43149 reads in 12237 unique sequences.
    ## Sample 131 - 32504 reads in 7959 unique sequences.
    ## Sample 132 - 28099 reads in 6147 unique sequences.
    ## Sample 133 - 36437 reads in 7574 unique sequences.
    ## Sample 134 - 42262 reads in 9678 unique sequences.
    ## Sample 135 - 46617 reads in 9059 unique sequences.
    ## Sample 136 - 51694 reads in 13604 unique sequences.
    ## Sample 137 - 27563 reads in 5811 unique sequences.
    ## Sample 138 - 32329 reads in 8351 unique sequences.
    ## Sample 139 - 38252 reads in 10306 unique sequences.
    ## Sample 140 - 40772 reads in 11812 unique sequences.
    ## Sample 141 - 44615 reads in 11729 unique sequences.
    ## Sample 142 - 49632 reads in 15527 unique sequences.
    ## Sample 143 - 27776 reads in 7666 unique sequences.
    ## Sample 144 - 33470 reads in 9735 unique sequences.
    ## Sample 145 - 41351 reads in 15563 unique sequences.
    ## Sample 146 - 44775 reads in 16491 unique sequences.
    ## Sample 147 - 42641 reads in 16345 unique sequences.
    ## Sample 148 - 50089 reads in 18497 unique sequences.
    ## Sample 149 - 26550 reads in 9948 unique sequences.
    ## Sample 150 - 33633 reads in 12512 unique sequences.
    ## Sample 151 - 42068 reads in 6103 unique sequences.
    ## Sample 152 - 41221 reads in 5886 unique sequences.
    ## Sample 153 - 44877 reads in 7790 unique sequences.
    ## Sample 154 - 44685 reads in 8410 unique sequences.
    ## Sample 155 - 26474 reads in 4279 unique sequences.
    ## Sample 156 - 25647 reads in 4072 unique sequences.
    ## Sample 157 - 46097 reads in 10815 unique sequences.
    ## Sample 158 - 42779 reads in 9720 unique sequences.
    ## Sample 159 - 52226 reads in 13119 unique sequences.
    ## Sample 160 - 53694 reads in 13505 unique sequences.
    ## Sample 161 - 29822 reads in 7215 unique sequences.
    ## Sample 162 - 31911 reads in 7261 unique sequences.
    ## Sample 163 - 49539 reads in 8771 unique sequences.
    ## Sample 164 - 44674 reads in 7349 unique sequences.
    ## Sample 165 - 54966 reads in 15835 unique sequences.
    ## Sample 166 - 54402 reads in 13198 unique sequences.
    ## Sample 167 - 33900 reads in 7490 unique sequences.
    ## Sample 168 - 29446 reads in 5747 unique sequences.
    ## Sample 169 - 58603 reads in 6820 unique sequences.
    ## Sample 170 - 46100 reads in 5945 unique sequences.
    ## Sample 171 - 54824 reads in 10100 unique sequences.
    ## Sample 172 - 48925 reads in 7903 unique sequences.
    ## Sample 173 - 36922 reads in 4525 unique sequences.
    ## Sample 174 - 26665 reads in 3876 unique sequences.
    ## Sample 175 - 37468 reads in 11663 unique sequences.
    ## Sample 176 - 41745 reads in 12093 unique sequences.
    ## Sample 177 - 38882 reads in 13137 unique sequences.
    ## Sample 178 - 44012 reads in 15009 unique sequences.
    ## Sample 179 - 23112 reads in 7184 unique sequences.
    ## Sample 180 - 24411 reads in 7048 unique sequences.
    ## Sample 181 - 27011 reads in 4248 unique sequences.
    ## Sample 182 - 37356 reads in 7050 unique sequences.
    ## Sample 183 - 31307 reads in 8152 unique sequences.
    ## Sample 184 - 46768 reads in 13109 unique sequences.
    ## Sample 185 - 22896 reads in 3841 unique sequences.
    ## Sample 186 - 30751 reads in 5844 unique sequences.
    ## Sample 187 - 49898 reads in 6346 unique sequences.
    ## Sample 188 - 35695 reads in 5666 unique sequences.
    ## Sample 189 - 51323 reads in 16365 unique sequences.
    ## Sample 190 - 38137 reads in 10831 unique sequences.
    ## Sample 191 - 37650 reads in 8164 unique sequences.
    ## Sample 192 - 23120 reads in 4064 unique sequences.
    ## Sample 193 - 35111 reads in 6570 unique sequences.
    ## Sample 194 - 48600 reads in 8148 unique sequences.
    ## Sample 195 - 39973 reads in 8004 unique sequences.
    ## Sample 196 - 55696 reads in 11288 unique sequences.
    ## Sample 197 - 21682 reads in 4352 unique sequences.
    ## Sample 198 - 35741 reads in 6639 unique sequences.
    ## Sample 199 - 51821 reads in 10625 unique sequences.
    ## Sample 200 - 47264 reads in 8934 unique sequences.
    ## Sample 201 - 60464 reads in 21236 unique sequences.
    ## Sample 202 - 57728 reads in 20534 unique sequences.
    ## Sample 203 - 43145 reads in 10124 unique sequences.
    ## Sample 204 - 41317 reads in 10835 unique sequences.
    ## Sample 205 - 44413 reads in 13267 unique sequences.
    ## Sample 206 - 45887 reads in 11977 unique sequences.
    ## Sample 207 - 47360 reads in 14782 unique sequences.
    ## Sample 208 - 59482 reads in 17485 unique sequences.
    ## Sample 209 - 33395 reads in 9592 unique sequences.
    ## Sample 210 - 37951 reads in 10039 unique sequences.
    ## Sample 211 - 38604 reads in 7271 unique sequences.
    ## Sample 212 - 47020 reads in 9495 unique sequences.
    ## Sample 213 - 41506 reads in 9129 unique sequences.
    ## Sample 214 - 61432 reads in 20961 unique sequences.
    ## Sample 215 - 27833 reads in 4987 unique sequences.
    ## Sample 216 - 45004 reads in 11889 unique sequences.
    ## Sample 217 - 36429 reads in 8537 unique sequences.
    ## Sample 218 - 40332 reads in 8697 unique sequences.
    ## Sample 219 - 35524 reads in 7987 unique sequences.
    ## Sample 220 - 43113 reads in 9799 unique sequences.
    ## Sample 221 - 23098 reads in 5474 unique sequences.
    ## Sample 222 - 25886 reads in 6129 unique sequences.
    ## Sample 223 - 39074 reads in 5271 unique sequences.
    ## Sample 224 - 37210 reads in 6430 unique sequences.
    ## Sample 225 - 51425 reads in 9277 unique sequences.
    ## Sample 226 - 42973 reads in 11791 unique sequences.
    ## Sample 227 - 31676 reads in 5203 unique sequences.
    ## Sample 228 - 26898 reads in 3711 unique sequences.
    ## Sample 229 - 34192 reads in 4394 unique sequences.
    ## Sample 230 - 44782 reads in 6923 unique sequences.
    ## Sample 231 - 39721 reads in 10302 unique sequences.
    ## Sample 232 - 45746 reads in 11687 unique sequences.
    ## Sample 233 - 26014 reads in 4299 unique sequences.
    ## Sample 234 - 25964 reads in 4640 unique sequences.
    ## Sample 235 - 39413 reads in 8074 unique sequences.
    ## Sample 236 - 43770 reads in 8431 unique sequences.
    ## Sample 237 - 48502 reads in 10790 unique sequences.
    ## Sample 238 - 48092 reads in 10951 unique sequences.
    ## Sample 239 - 28727 reads in 6046 unique sequences.
    ## Sample 240 - 29293 reads in 6109 unique sequences.
    ## Sample 241 - 38577 reads in 5553 unique sequences.
    ## Sample 242 - 43137 reads in 6217 unique sequences.
    ## Sample 243 - 44529 reads in 7063 unique sequences.
    ## Sample 244 - 49960 reads in 12354 unique sequences.
    ## Sample 245 - 24860 reads in 3703 unique sequences.
    ## Sample 246 - 29488 reads in 5246 unique sequences.
    ## Sample 247 - 7973 reads in 2421 unique sequences.
    ## Sample 248 - 40416 reads in 7486 unique sequences.
    ## Sample 249 - 8135 reads in 2594 unique sequences.
    ## Sample 250 - 46338 reads in 10612 unique sequences.
    ## Sample 251 - 6257 reads in 1880 unique sequences.
    ## Sample 252 - 27915 reads in 5863 unique sequences.
    ## Sample 253 - 26851 reads in 7542 unique sequences.
    ## Sample 254 - 26987 reads in 8905 unique sequences.
    ## Sample 255 - 32943 reads in 9440 unique sequences.
    ## Sample 256 - 32363 reads in 11281 unique sequences.
    ## Sample 257 - 20715 reads in 5580 unique sequences.
    ## Sample 258 - 21466 reads in 6697 unique sequences.
    ## Sample 259 - 33791 reads in 18271 unique sequences.
    ## Sample 260 - 46535 reads in 25300 unique sequences.
    ## Sample 261 - 56133 reads in 28339 unique sequences.
    ## Sample 262 - 31747 reads in 16512 unique sequences.
    ## Sample 263 - 37944 reads in 19926 unique sequences.
    ## Sample 264 - 33812 reads in 17542 unique sequences.
    ## Sample 265 - 24122 reads in 10992 unique sequences.
    ## Sample 266 - 39312 reads in 21166 unique sequences.
    ## Sample 267 - 33689 reads in 18614 unique sequences.
    ## Sample 268 - 34797 reads in 18649 unique sequences.
    ## Sample 269 - 40345 reads in 21761 unique sequences.
    ## Sample 270 - 38556 reads in 17699 unique sequences.
    ## Sample 271 - 19919 reads in 8054 unique sequences.
    ## Sample 272 - 34790 reads in 18212 unique sequences.
    ## Sample 273 - 32849 reads in 16740 unique sequences.
    ## Sample 274 - 32228 reads in 15094 unique sequences.
    ## Sample 275 - 45131 reads in 23740 unique sequences.
    ## Sample 276 - 46620 reads in 24642 unique sequences.
    ## Sample 277 - 32535 reads in 15967 unique sequences.
    ## Sample 278 - 39235 reads in 20616 unique sequences.
    ## Sample 279 - 45445 reads in 22948 unique sequences.
    ## Sample 280 - 31907 reads in 17752 unique sequences.
    ## Sample 281 - 38041 reads in 21720 unique sequences.
    ## Sample 282 - 32718 reads in 16420 unique sequences.
    ## Sample 283 - 28596 reads in 14060 unique sequences.
    ## Sample 284 - 37839 reads in 19775 unique sequences.
    ## Sample 285 - 39764 reads in 21590 unique sequences.
    ## Sample 286 - 27422 reads in 12265 unique sequences.
    ## Sample 287 - 32051 reads in 13831 unique sequences.
    ## Sample 288 - 20355 reads in 9080 unique sequences.

``` r
# reverses
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 39017 reads in 6201 unique sequences.
    ## Sample 2 - 36937 reads in 7005 unique sequences.
    ## Sample 3 - 38973 reads in 11942 unique sequences.
    ## Sample 4 - 46766 reads in 14530 unique sequences.
    ## Sample 5 - 26614 reads in 5167 unique sequences.
    ## Sample 6 - 32782 reads in 7764 unique sequences.
    ## Sample 7 - 40836 reads in 6399 unique sequences.
    ## Sample 8 - 39476 reads in 6450 unique sequences.
    ## Sample 9 - 44111 reads in 15916 unique sequences.
    ## Sample 10 - 44482 reads in 11041 unique sequences.
    ## Sample 11 - 25363 reads in 6326 unique sequences.
    ## Sample 12 - 29479 reads in 5099 unique sequences.
    ## Sample 13 - 39151 reads in 8189 unique sequences.
    ## Sample 14 - 53206 reads in 9272 unique sequences.
    ## Sample 15 - 41422 reads in 18322 unique sequences.
    ## Sample 16 - 61733 reads in 15242 unique sequences.
    ## Sample 17 - 25791 reads in 5279 unique sequences.
    ## Sample 18 - 40901 reads in 7903 unique sequences.
    ## Sample 19 - 43354 reads in 6939 unique sequences.
    ## Sample 20 - 52145 reads in 8429 unique sequences.
    ## Sample 21 - 45145 reads in 7863 unique sequences.
    ## Sample 22 - 57210 reads in 9349 unique sequences.
    ## Sample 23 - 28576 reads in 4843 unique sequences.
    ## Sample 24 - 35955 reads in 5937 unique sequences.
    ## Sample 25 - 38375 reads in 8565 unique sequences.
    ## Sample 26 - 49693 reads in 11014 unique sequences.
    ## Sample 27 - 43738 reads in 10431 unique sequences.
    ## Sample 28 - 58206 reads in 14023 unique sequences.
    ## Sample 29 - 29354 reads in 6286 unique sequences.
    ## Sample 30 - 37192 reads in 8144 unique sequences.
    ## Sample 31 - 38644 reads in 6969 unique sequences.
    ## Sample 32 - 41439 reads in 8142 unique sequences.
    ## Sample 33 - 41280 reads in 9118 unique sequences.
    ## Sample 34 - 48845 reads in 13604 unique sequences.
    ## Sample 35 - 21675 reads in 4128 unique sequences.
    ## Sample 36 - 30671 reads in 5718 unique sequences.
    ## Sample 37 - 43037 reads in 8012 unique sequences.
    ## Sample 38 - 42905 reads in 7437 unique sequences.
    ## Sample 39 - 45206 reads in 9738 unique sequences.
    ## Sample 40 - 45682 reads in 8055 unique sequences.
    ## Sample 41 - 25655 reads in 4944 unique sequences.
    ## Sample 42 - 28199 reads in 4866 unique sequences.
    ## Sample 43 - 40184 reads in 8502 unique sequences.
    ## Sample 44 - 43380 reads in 8887 unique sequences.
    ## Sample 45 - 41199 reads in 17175 unique sequences.
    ## Sample 46 - 45616 reads in 13757 unique sequences.
    ## Sample 47 - 23246 reads in 5821 unique sequences.
    ## Sample 48 - 33058 reads in 7078 unique sequences.
    ## Sample 49 - 38693 reads in 5591 unique sequences.
    ## Sample 50 - 38741 reads in 6249 unique sequences.
    ## Sample 51 - 42978 reads in 6930 unique sequences.
    ## Sample 52 - 54868 reads in 9640 unique sequences.
    ## Sample 53 - 27696 reads in 4378 unique sequences.
    ## Sample 54 - 31577 reads in 5224 unique sequences.
    ## Sample 55 - 44388 reads in 8466 unique sequences.
    ## Sample 56 - 41192 reads in 7731 unique sequences.
    ## Sample 57 - 50886 reads in 12949 unique sequences.
    ## Sample 58 - 45257 reads in 12357 unique sequences.
    ## Sample 59 - 29400 reads in 6102 unique sequences.
    ## Sample 60 - 28052 reads in 5520 unique sequences.
    ## Sample 61 - 52070 reads in 8470 unique sequences.
    ## Sample 62 - 42583 reads in 8807 unique sequences.
    ## Sample 63 - 61516 reads in 13201 unique sequences.
    ## Sample 64 - 48410 reads in 12169 unique sequences.
    ## Sample 65 - 37551 reads in 6409 unique sequences.
    ## Sample 66 - 25714 reads in 5836 unique sequences.
    ## Sample 67 - 60579 reads in 7745 unique sequences.
    ## Sample 68 - 46753 reads in 6008 unique sequences.
    ## Sample 69 - 59452 reads in 15100 unique sequences.
    ## Sample 70 - 49504 reads in 9961 unique sequences.
    ## Sample 71 - 32504 reads in 4707 unique sequences.
    ## Sample 72 - 31783 reads in 4024 unique sequences.
    ## Sample 73 - 53713 reads in 8374 unique sequences.
    ## Sample 74 - 53867 reads in 7246 unique sequences.
    ## Sample 75 - 58977 reads in 10806 unique sequences.
    ## Sample 76 - 62290 reads in 8996 unique sequences.
    ## Sample 77 - 33817 reads in 5451 unique sequences.
    ## Sample 78 - 39820 reads in 5433 unique sequences.
    ## Sample 79 - 50850 reads in 7795 unique sequences.
    ## Sample 80 - 56367 reads in 9986 unique sequences.
    ## Sample 81 - 52403 reads in 12358 unique sequences.
    ## Sample 82 - 58087 reads in 18247 unique sequences.
    ## Sample 83 - 33499 reads in 5690 unique sequences.
    ## Sample 84 - 35017 reads in 7029 unique sequences.
    ## Sample 85 - 23523 reads in 3453 unique sequences.
    ## Sample 86 - 41036 reads in 5269 unique sequences.
    ## Sample 87 - 21703 reads in 4204 unique sequences.
    ## Sample 88 - 38492 reads in 7158 unique sequences.
    ## Sample 89 - 11867 reads in 2093 unique sequences.
    ## Sample 90 - 24073 reads in 3587 unique sequences.
    ## Sample 91 - 35678 reads in 5227 unique sequences.
    ## Sample 92 - 40378 reads in 6007 unique sequences.
    ## Sample 93 - 37525 reads in 7982 unique sequences.
    ## Sample 94 - 51034 reads in 10155 unique sequences.
    ## Sample 95 - 23132 reads in 3780 unique sequences.
    ## Sample 96 - 32169 reads in 4759 unique sequences.
    ## Sample 97 - 44286 reads in 5480 unique sequences.
    ## Sample 98 - 42665 reads in 5654 unique sequences.
    ## Sample 99 - 45018 reads in 7052 unique sequences.
    ## Sample 100 - 52002 reads in 10304 unique sequences.
    ## Sample 101 - 27881 reads in 3724 unique sequences.
    ## Sample 102 - 34061 reads in 4414 unique sequences.
    ## Sample 103 - 41451 reads in 5500 unique sequences.
    ## Sample 104 - 37127 reads in 6931 unique sequences.
    ## Sample 105 - 40735 reads in 9622 unique sequences.
    ## Sample 106 - 40887 reads in 11008 unique sequences.
    ## Sample 107 - 24513 reads in 3873 unique sequences.
    ## Sample 108 - 27344 reads in 5416 unique sequences.
    ## Sample 109 - 44606 reads in 9177 unique sequences.
    ## Sample 110 - 42728 reads in 9688 unique sequences.
    ## Sample 111 - 49920 reads in 15592 unique sequences.
    ## Sample 112 - 49529 reads in 16996 unique sequences.
    ## Sample 113 - 28161 reads in 6534 unique sequences.
    ## Sample 114 - 30352 reads in 7546 unique sequences.
    ## Sample 115 - 59561 reads in 11842 unique sequences.
    ## Sample 116 - 44247 reads in 14416 unique sequences.
    ## Sample 117 - 55172 reads in 21292 unique sequences.
    ## Sample 118 - 43140 reads in 23697 unique sequences.
    ## Sample 119 - 32860 reads in 8402 unique sequences.
    ## Sample 120 - 29931 reads in 12678 unique sequences.
    ## Sample 121 - 36762 reads in 8309 unique sequences.
    ## Sample 122 - 36650 reads in 6609 unique sequences.
    ## Sample 123 - 40851 reads in 9403 unique sequences.
    ## Sample 124 - 43349 reads in 9818 unique sequences.
    ## Sample 125 - 24801 reads in 5240 unique sequences.
    ## Sample 126 - 27020 reads in 4649 unique sequences.
    ## Sample 127 - 53712 reads in 12344 unique sequences.
    ## Sample 128 - 40677 reads in 8474 unique sequences.
    ## Sample 129 - 54761 reads in 16041 unique sequences.
    ## Sample 130 - 43149 reads in 12919 unique sequences.
    ## Sample 131 - 32504 reads in 7636 unique sequences.
    ## Sample 132 - 28099 reads in 6054 unique sequences.
    ## Sample 133 - 36437 reads in 6837 unique sequences.
    ## Sample 134 - 42262 reads in 8210 unique sequences.
    ## Sample 135 - 46617 reads in 8554 unique sequences.
    ## Sample 136 - 51694 reads in 14412 unique sequences.
    ## Sample 137 - 27563 reads in 5011 unique sequences.
    ## Sample 138 - 32329 reads in 6478 unique sequences.
    ## Sample 139 - 38252 reads in 7635 unique sequences.
    ## Sample 140 - 40772 reads in 9007 unique sequences.
    ## Sample 141 - 44615 reads in 9615 unique sequences.
    ## Sample 142 - 49632 reads in 14878 unique sequences.
    ## Sample 143 - 27776 reads in 5414 unique sequences.
    ## Sample 144 - 33470 reads in 7068 unique sequences.
    ## Sample 145 - 41351 reads in 11164 unique sequences.
    ## Sample 146 - 44775 reads in 10835 unique sequences.
    ## Sample 147 - 42641 reads in 14030 unique sequences.
    ## Sample 148 - 50089 reads in 15185 unique sequences.
    ## Sample 149 - 26550 reads in 6985 unique sequences.
    ## Sample 150 - 33633 reads in 8167 unique sequences.
    ## Sample 151 - 42068 reads in 6077 unique sequences.
    ## Sample 152 - 41221 reads in 6262 unique sequences.
    ## Sample 153 - 44877 reads in 7562 unique sequences.
    ## Sample 154 - 44685 reads in 9786 unique sequences.
    ## Sample 155 - 26474 reads in 4160 unique sequences.
    ## Sample 156 - 25647 reads in 4084 unique sequences.
    ## Sample 157 - 46097 reads in 9453 unique sequences.
    ## Sample 158 - 42779 reads in 8923 unique sequences.
    ## Sample 159 - 52226 reads in 14142 unique sequences.
    ## Sample 160 - 53694 reads in 14651 unique sequences.
    ## Sample 161 - 29822 reads in 6410 unique sequences.
    ## Sample 162 - 31911 reads in 6476 unique sequences.
    ## Sample 163 - 49539 reads in 7665 unique sequences.
    ## Sample 164 - 44674 reads in 7146 unique sequences.
    ## Sample 165 - 54966 reads in 20738 unique sequences.
    ## Sample 166 - 54402 reads in 16541 unique sequences.
    ## Sample 167 - 33900 reads in 7556 unique sequences.
    ## Sample 168 - 29446 reads in 5828 unique sequences.
    ## Sample 169 - 58603 reads in 7349 unique sequences.
    ## Sample 170 - 46100 reads in 6598 unique sequences.
    ## Sample 171 - 54824 reads in 11247 unique sequences.
    ## Sample 172 - 48925 reads in 9577 unique sequences.
    ## Sample 173 - 36922 reads in 4794 unique sequences.
    ## Sample 174 - 26665 reads in 4239 unique sequences.
    ## Sample 175 - 37468 reads in 10376 unique sequences.
    ## Sample 176 - 41745 reads in 11319 unique sequences.
    ## Sample 177 - 38882 reads in 13752 unique sequences.
    ## Sample 178 - 44012 reads in 14411 unique sequences.
    ## Sample 179 - 23112 reads in 6252 unique sequences.
    ## Sample 180 - 24411 reads in 6285 unique sequences.
    ## Sample 181 - 27011 reads in 4571 unique sequences.
    ## Sample 182 - 37356 reads in 7905 unique sequences.
    ## Sample 183 - 31307 reads in 8578 unique sequences.
    ## Sample 184 - 46768 reads in 18503 unique sequences.
    ## Sample 185 - 22896 reads in 4003 unique sequences.
    ## Sample 186 - 30751 reads in 6279 unique sequences.
    ## Sample 187 - 49898 reads in 6874 unique sequences.
    ## Sample 188 - 35695 reads in 6534 unique sequences.
    ## Sample 189 - 51323 reads in 21832 unique sequences.
    ## Sample 190 - 38137 reads in 14025 unique sequences.
    ## Sample 191 - 37650 reads in 9871 unique sequences.
    ## Sample 192 - 23120 reads in 4399 unique sequences.
    ## Sample 193 - 35111 reads in 5640 unique sequences.
    ## Sample 194 - 48600 reads in 8285 unique sequences.
    ## Sample 195 - 39973 reads in 6935 unique sequences.
    ## Sample 196 - 55696 reads in 12035 unique sequences.
    ## Sample 197 - 21682 reads in 3730 unique sequences.
    ## Sample 198 - 35741 reads in 6546 unique sequences.
    ## Sample 199 - 51821 reads in 11608 unique sequences.
    ## Sample 200 - 47264 reads in 9101 unique sequences.
    ## Sample 201 - 60464 reads in 25859 unique sequences.
    ## Sample 202 - 57728 reads in 27865 unique sequences.
    ## Sample 203 - 43145 reads in 11662 unique sequences.
    ## Sample 204 - 41317 reads in 14249 unique sequences.
    ## Sample 205 - 44413 reads in 10110 unique sequences.
    ## Sample 206 - 45887 reads in 9659 unique sequences.
    ## Sample 207 - 47360 reads in 14229 unique sequences.
    ## Sample 208 - 59482 reads in 17973 unique sequences.
    ## Sample 209 - 33395 reads in 7106 unique sequences.
    ## Sample 210 - 37951 reads in 7771 unique sequences.
    ## Sample 211 - 38604 reads in 8266 unique sequences.
    ## Sample 212 - 47020 reads in 11850 unique sequences.
    ## Sample 213 - 41506 reads in 9911 unique sequences.
    ## Sample 214 - 61432 reads in 31219 unique sequences.
    ## Sample 215 - 27833 reads in 5495 unique sequences.
    ## Sample 216 - 45004 reads in 16329 unique sequences.
    ## Sample 217 - 36429 reads in 7110 unique sequences.
    ## Sample 218 - 40332 reads in 7600 unique sequences.
    ## Sample 219 - 35524 reads in 6987 unique sequences.
    ## Sample 220 - 43113 reads in 8733 unique sequences.
    ## Sample 221 - 23098 reads in 4543 unique sequences.
    ## Sample 222 - 25886 reads in 5123 unique sequences.
    ## Sample 223 - 39074 reads in 5944 unique sequences.
    ## Sample 224 - 37210 reads in 8443 unique sequences.
    ## Sample 225 - 51425 reads in 11066 unique sequences.
    ## Sample 226 - 42973 reads in 16893 unique sequences.
    ## Sample 227 - 31676 reads in 5824 unique sequences.
    ## Sample 228 - 26898 reads in 3978 unique sequences.
    ## Sample 229 - 34192 reads in 4377 unique sequences.
    ## Sample 230 - 44782 reads in 7514 unique sequences.
    ## Sample 231 - 39721 reads in 9717 unique sequences.
    ## Sample 232 - 45746 reads in 14161 unique sequences.
    ## Sample 233 - 26014 reads in 3973 unique sequences.
    ## Sample 234 - 25964 reads in 5072 unique sequences.
    ## Sample 235 - 39413 reads in 7845 unique sequences.
    ## Sample 236 - 43770 reads in 8414 unique sequences.
    ## Sample 237 - 48502 reads in 12147 unique sequences.
    ## Sample 238 - 48092 reads in 11365 unique sequences.
    ## Sample 239 - 28727 reads in 5640 unique sequences.
    ## Sample 240 - 29293 reads in 5900 unique sequences.
    ## Sample 241 - 38577 reads in 5654 unique sequences.
    ## Sample 242 - 43137 reads in 6879 unique sequences.
    ## Sample 243 - 44529 reads in 7989 unique sequences.
    ## Sample 244 - 49960 reads in 17488 unique sequences.
    ## Sample 245 - 24860 reads in 3649 unique sequences.
    ## Sample 246 - 29488 reads in 5657 unique sequences.
    ## Sample 247 - 7973 reads in 1779 unique sequences.
    ## Sample 248 - 40416 reads in 6579 unique sequences.
    ## Sample 249 - 8135 reads in 2066 unique sequences.
    ## Sample 250 - 46338 reads in 11486 unique sequences.
    ## Sample 251 - 6257 reads in 1347 unique sequences.
    ## Sample 252 - 27915 reads in 4945 unique sequences.
    ## Sample 253 - 26851 reads in 5496 unique sequences.
    ## Sample 254 - 26987 reads in 8418 unique sequences.
    ## Sample 255 - 32943 reads in 7708 unique sequences.
    ## Sample 256 - 32363 reads in 11467 unique sequences.
    ## Sample 257 - 20715 reads in 3922 unique sequences.
    ## Sample 258 - 21466 reads in 6347 unique sequences.
    ## Sample 259 - 33791 reads in 27267 unique sequences.
    ## Sample 260 - 46535 reads in 33386 unique sequences.
    ## Sample 261 - 56133 reads in 35128 unique sequences.
    ## Sample 262 - 31747 reads in 13482 unique sequences.
    ## Sample 263 - 37944 reads in 17263 unique sequences.
    ## Sample 264 - 33812 reads in 15498 unique sequences.
    ## Sample 265 - 24122 reads in 23382 unique sequences.
    ## Sample 266 - 39312 reads in 32543 unique sequences.
    ## Sample 267 - 33689 reads in 32433 unique sequences.
    ## Sample 268 - 34797 reads in 14215 unique sequences.
    ## Sample 269 - 40345 reads in 16953 unique sequences.
    ## Sample 270 - 38556 reads in 13525 unique sequences.
    ## Sample 271 - 19919 reads in 6386 unique sequences.
    ## Sample 272 - 34790 reads in 15905 unique sequences.
    ## Sample 273 - 32849 reads in 14381 unique sequences.
    ## Sample 274 - 32228 reads in 20670 unique sequences.
    ## Sample 275 - 45131 reads in 36562 unique sequences.
    ## Sample 276 - 46620 reads in 41991 unique sequences.
    ## Sample 277 - 32535 reads in 12513 unique sequences.
    ## Sample 278 - 39235 reads in 17505 unique sequences.
    ## Sample 279 - 45445 reads in 17663 unique sequences.
    ## Sample 280 - 31907 reads in 28925 unique sequences.
    ## Sample 281 - 38041 reads in 29623 unique sequences.
    ## Sample 282 - 32718 reads in 22899 unique sequences.
    ## Sample 283 - 28596 reads in 11294 unique sequences.
    ## Sample 284 - 37839 reads in 17071 unique sequences.
    ## Sample 285 - 39764 reads in 17776 unique sequences.
    ## Sample 286 - 27422 reads in 8872 unique sequences.
    ## Sample 287 - 32051 reads in 9769 unique sequences.
    ## Sample 288 - 20355 reads in 6204 unique sequences.

Merge paired end reads

``` r
# too many reads are getting filtered out during the merge step!
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE,
                      minOverlap = 10,
                      maxMismatch = 0) # adding parameters to see if that helps...
```

    ## 37891 paired-reads (in 9 unique pairings) successfully merged out of 38388 (in 23 pairings) input.

    ## 34387 paired-reads (in 10 unique pairings) successfully merged out of 36087 (in 38 pairings) input.

    ## 30862 paired-reads (in 30 unique pairings) successfully merged out of 36589 (in 86 pairings) input.

    ## 36210 paired-reads (in 41 unique pairings) successfully merged out of 43844 (in 156 pairings) input.

    ## 24976 paired-reads (in 11 unique pairings) successfully merged out of 25924 (in 30 pairings) input.

    ## 28492 paired-reads (in 13 unique pairings) successfully merged out of 31492 (in 35 pairings) input.

    ## 39847 paired-reads (in 10 unique pairings) successfully merged out of 40437 (in 53 pairings) input.

    ## 37178 paired-reads (in 13 unique pairings) successfully merged out of 38824 (in 61 pairings) input.

    ## 33461 paired-reads (in 39 unique pairings) successfully merged out of 41249 (in 162 pairings) input.

    ## 37222 paired-reads (in 58 unique pairings) successfully merged out of 42485 (in 154 pairings) input.

    ## 22354 paired-reads (in 7 unique pairings) successfully merged out of 24456 (in 38 pairings) input.

    ## 27551 paired-reads (in 12 unique pairings) successfully merged out of 28749 (in 55 pairings) input.

    ## 36573 paired-reads (in 14 unique pairings) successfully merged out of 37143 (in 47 pairings) input.

    ## 49601 paired-reads (in 20 unique pairings) successfully merged out of 51947 (in 63 pairings) input.

    ## 28449 paired-reads (in 52 unique pairings) successfully merged out of 37577 (in 164 pairings) input.

    ## 45433 paired-reads (in 59 unique pairings) successfully merged out of 57879 (in 186 pairings) input.

    ## 24161 paired-reads (in 11 unique pairings) successfully merged out of 24852 (in 38 pairings) input.

    ## 36246 paired-reads (in 14 unique pairings) successfully merged out of 39544 (in 59 pairings) input.

    ## 42341 paired-reads (in 14 unique pairings) successfully merged out of 42810 (in 34 pairings) input.

    ## 50900 paired-reads (in 16 unique pairings) successfully merged out of 51400 (in 41 pairings) input.

    ## 43408 paired-reads (in 32 unique pairings) successfully merged out of 44455 (in 75 pairings) input.

    ## 55332 paired-reads (in 24 unique pairings) successfully merged out of 56458 (in 76 pairings) input.

    ## 27855 paired-reads (in 13 unique pairings) successfully merged out of 28188 (in 30 pairings) input.

    ## 35103 paired-reads (in 14 unique pairings) successfully merged out of 35512 (in 36 pairings) input.

    ## 36844 paired-reads (in 33 unique pairings) successfully merged out of 37616 (in 56 pairings) input.

    ## 47302 paired-reads (in 46 unique pairings) successfully merged out of 48419 (in 83 pairings) input.

    ## 39827 paired-reads (in 44 unique pairings) successfully merged out of 40889 (in 96 pairings) input.

    ## 51093 paired-reads (in 61 unique pairings) successfully merged out of 55800 (in 133 pairings) input.

    ## 28175 paired-reads (in 29 unique pairings) successfully merged out of 28752 (in 50 pairings) input.

    ## 35401 paired-reads (in 36 unique pairings) successfully merged out of 36115 (in 61 pairings) input.

    ## 37230 paired-reads (in 18 unique pairings) successfully merged out of 37938 (in 43 pairings) input.

    ## 39116 paired-reads (in 30 unique pairings) successfully merged out of 40384 (in 64 pairings) input.

    ## 37214 paired-reads (in 38 unique pairings) successfully merged out of 39971 (in 89 pairings) input.

    ## 40572 paired-reads (in 49 unique pairings) successfully merged out of 46468 (in 120 pairings) input.

    ## 20793 paired-reads (in 13 unique pairings) successfully merged out of 21208 (in 31 pairings) input.

    ## 29387 paired-reads (in 26 unique pairings) successfully merged out of 30016 (in 57 pairings) input.

    ## 40635 paired-reads (in 20 unique pairings) successfully merged out of 42240 (in 41 pairings) input.

    ## 40610 paired-reads (in 21 unique pairings) successfully merged out of 42205 (in 45 pairings) input.

    ## 41032 paired-reads (in 34 unique pairings) successfully merged out of 43891 (in 87 pairings) input.

    ## 42365 paired-reads (in 33 unique pairings) successfully merged out of 44815 (in 74 pairings) input.

    ## 24185 paired-reads (in 19 unique pairings) successfully merged out of 25169 (in 40 pairings) input.

    ## 26676 paired-reads (in 16 unique pairings) successfully merged out of 27748 (in 33 pairings) input.

    ## 38052 paired-reads (in 19 unique pairings) successfully merged out of 39148 (in 67 pairings) input.

    ## 39756 paired-reads (in 26 unique pairings) successfully merged out of 42163 (in 68 pairings) input.

    ## 24287 paired-reads (in 46 unique pairings) successfully merged out of 36955 (in 196 pairings) input.

    ## 27781 paired-reads (in 61 unique pairings) successfully merged out of 41301 (in 226 pairings) input.

    ## 21017 paired-reads (in 16 unique pairings) successfully merged out of 22331 (in 68 pairings) input.

    ## 29039 paired-reads (in 18 unique pairings) successfully merged out of 31880 (in 76 pairings) input.

    ## 38084 paired-reads (in 10 unique pairings) successfully merged out of 38309 (in 28 pairings) input.

    ## 37981 paired-reads (in 15 unique pairings) successfully merged out of 38288 (in 38 pairings) input.

    ## 41472 paired-reads (in 19 unique pairings) successfully merged out of 42440 (in 46 pairings) input.

    ## 52193 paired-reads (in 31 unique pairings) successfully merged out of 53791 (in 74 pairings) input.

    ## 27145 paired-reads (in 10 unique pairings) successfully merged out of 27412 (in 29 pairings) input.

    ## 30730 paired-reads (in 19 unique pairings) successfully merged out of 31012 (in 38 pairings) input.

    ## 42606 paired-reads (in 29 unique pairings) successfully merged out of 43523 (in 76 pairings) input.

    ## 39559 paired-reads (in 28 unique pairings) successfully merged out of 40461 (in 72 pairings) input.

    ## 41354 paired-reads (in 58 unique pairings) successfully merged out of 48192 (in 151 pairings) input.

    ## 37485 paired-reads (in 46 unique pairings) successfully merged out of 43018 (in 142 pairings) input.

    ## 27766 paired-reads (in 27 unique pairings) successfully merged out of 28663 (in 63 pairings) input.

    ## 26756 paired-reads (in 15 unique pairings) successfully merged out of 27451 (in 47 pairings) input.

    ## 50793 paired-reads (in 13 unique pairings) successfully merged out of 51252 (in 29 pairings) input.

    ## 40896 paired-reads (in 16 unique pairings) successfully merged out of 41636 (in 37 pairings) input.

    ## 56225 paired-reads (in 43 unique pairings) successfully merged out of 59916 (in 106 pairings) input.

    ## 42808 paired-reads (in 35 unique pairings) successfully merged out of 46546 (in 112 pairings) input.

    ## 36466 paired-reads (in 12 unique pairings) successfully merged out of 36786 (in 27 pairings) input.

    ## 24257 paired-reads (in 12 unique pairings) successfully merged out of 24944 (in 35 pairings) input.

    ## 59722 paired-reads (in 7 unique pairings) successfully merged out of 59890 (in 25 pairings) input.

    ## 45863 paired-reads (in 7 unique pairings) successfully merged out of 46014 (in 26 pairings) input.

    ## 50816 paired-reads (in 25 unique pairings) successfully merged out of 56481 (in 87 pairings) input.

    ## 44751 paired-reads (in 25 unique pairings) successfully merged out of 47794 (in 92 pairings) input.

    ## 31628 paired-reads (in 7 unique pairings) successfully merged out of 31850 (in 32 pairings) input.

    ## 31272 paired-reads (in 9 unique pairings) successfully merged out of 31396 (in 28 pairings) input.

    ## 52249 paired-reads (in 23 unique pairings) successfully merged out of 52832 (in 47 pairings) input.

    ## 53211 paired-reads (in 6 unique pairings) successfully merged out of 53373 (in 19 pairings) input.

    ## 55185 paired-reads (in 20 unique pairings) successfully merged out of 55627 (in 66 pairings) input.

    ## 60907 paired-reads (in 13 unique pairings) successfully merged out of 61727 (in 43 pairings) input.

    ## 32775 paired-reads (in 11 unique pairings) successfully merged out of 33064 (in 28 pairings) input.

    ## 39391 paired-reads (in 7 unique pairings) successfully merged out of 39501 (in 21 pairings) input.

    ## 49184 paired-reads (in 8 unique pairings) successfully merged out of 50193 (in 31 pairings) input.

    ## 53627 paired-reads (in 19 unique pairings) successfully merged out of 55388 (in 67 pairings) input.

    ## 42132 paired-reads (in 45 unique pairings) successfully merged out of 49942 (in 142 pairings) input.

    ## 43522 paired-reads (in 44 unique pairings) successfully merged out of 54602 (in 170 pairings) input.

    ## 31172 paired-reads (in 8 unique pairings) successfully merged out of 32820 (in 25 pairings) input.

    ## 32180 paired-reads (in 15 unique pairings) successfully merged out of 34092 (in 60 pairings) input.

    ## 23135 paired-reads (in 7 unique pairings) successfully merged out of 23296 (in 26 pairings) input.

    ## 40693 paired-reads (in 7 unique pairings) successfully merged out of 40795 (in 23 pairings) input.

    ## 19554 paired-reads (in 16 unique pairings) successfully merged out of 20916 (in 59 pairings) input.

    ## 35439 paired-reads (in 23 unique pairings) successfully merged out of 37567 (in 65 pairings) input.

    ## 11358 paired-reads (in 7 unique pairings) successfully merged out of 11631 (in 23 pairings) input.

    ## 23398 paired-reads (in 5 unique pairings) successfully merged out of 23757 (in 21 pairings) input.

    ## 35148 paired-reads (in 7 unique pairings) successfully merged out of 35331 (in 17 pairings) input.

    ## 39674 paired-reads (in 10 unique pairings) successfully merged out of 39879 (in 29 pairings) input.

    ## 33445 paired-reads (in 32 unique pairings) successfully merged out of 36454 (in 103 pairings) input.

    ## 47225 paired-reads (in 20 unique pairings) successfully merged out of 49735 (in 78 pairings) input.

    ## 22483 paired-reads (in 6 unique pairings) successfully merged out of 22820 (in 18 pairings) input.

    ## 31720 paired-reads (in 6 unique pairings) successfully merged out of 31900 (in 20 pairings) input.

    ## 43865 paired-reads (in 10 unique pairings) successfully merged out of 44047 (in 31 pairings) input.

    ## 42014 paired-reads (in 8 unique pairings) successfully merged out of 42386 (in 26 pairings) input.

    ## 41527 paired-reads (in 32 unique pairings) successfully merged out of 44107 (in 109 pairings) input.

    ## 44878 paired-reads (in 41 unique pairings) successfully merged out of 50272 (in 105 pairings) input.

    ## 26962 paired-reads (in 6 unique pairings) successfully merged out of 27486 (in 24 pairings) input.

    ## 33739 paired-reads (in 10 unique pairings) successfully merged out of 33851 (in 24 pairings) input.

    ## 41207 paired-reads (in 6 unique pairings) successfully merged out of 41282 (in 16 pairings) input.

    ## 35554 paired-reads (in 13 unique pairings) successfully merged out of 36326 (in 37 pairings) input.

    ## 35151 paired-reads (in 32 unique pairings) successfully merged out of 39230 (in 87 pairings) input.

    ## 33261 paired-reads (in 32 unique pairings) successfully merged out of 38689 (in 94 pairings) input.

    ## 23821 paired-reads (in 6 unique pairings) successfully merged out of 24242 (in 17 pairings) input.

    ## 25656 paired-reads (in 10 unique pairings) successfully merged out of 26634 (in 29 pairings) input.

    ## 42374 paired-reads (in 42 unique pairings) successfully merged out of 43294 (in 91 pairings) input.

    ## 39667 paired-reads (in 31 unique pairings) successfully merged out of 41412 (in 89 pairings) input.

    ## 40859 paired-reads (in 65 unique pairings) successfully merged out of 47019 (in 159 pairings) input.

    ## 39111 paired-reads (in 56 unique pairings) successfully merged out of 46309 (in 160 pairings) input.

    ## 26011 paired-reads (in 44 unique pairings) successfully merged out of 27181 (in 91 pairings) input.

    ## 27404 paired-reads (in 25 unique pairings) successfully merged out of 29261 (in 64 pairings) input.

    ## 54099 paired-reads (in 16 unique pairings) successfully merged out of 58191 (in 70 pairings) input.

    ## 38847 paired-reads (in 8 unique pairings) successfully merged out of 41885 (in 58 pairings) input.

    ## 39283 paired-reads (in 74 unique pairings) successfully merged out of 50933 (in 244 pairings) input.

    ## 27000 paired-reads (in 31 unique pairings) successfully merged out of 37496 (in 146 pairings) input.

    ## 28052 paired-reads (in 10 unique pairings) successfully merged out of 29212 (in 72 pairings) input.

    ## 24097 paired-reads (in 13 unique pairings) successfully merged out of 26449 (in 62 pairings) input.

    ## 35629 paired-reads (in 16 unique pairings) successfully merged out of 35980 (in 31 pairings) input.

    ## 34819 paired-reads (in 17 unique pairings) successfully merged out of 35632 (in 36 pairings) input.

    ## 38296 paired-reads (in 30 unique pairings) successfully merged out of 38674 (in 70 pairings) input.

    ## 38867 paired-reads (in 31 unique pairings) successfully merged out of 41496 (in 76 pairings) input.

    ## 23883 paired-reads (in 16 unique pairings) successfully merged out of 24090 (in 30 pairings) input.

    ## 25763 paired-reads (in 15 unique pairings) successfully merged out of 26300 (in 32 pairings) input.

    ## 51593 paired-reads (in 25 unique pairings) successfully merged out of 52240 (in 60 pairings) input.

    ## 38692 paired-reads (in 18 unique pairings) successfully merged out of 39678 (in 49 pairings) input.

    ## 46962 paired-reads (in 48 unique pairings) successfully merged out of 51776 (in 163 pairings) input.

    ## 33619 paired-reads (in 45 unique pairings) successfully merged out of 40573 (in 146 pairings) input.

    ## 29920 paired-reads (in 27 unique pairings) successfully merged out of 31165 (in 69 pairings) input.

    ## 26040 paired-reads (in 13 unique pairings) successfully merged out of 27324 (in 54 pairings) input.

    ## 35190 paired-reads (in 27 unique pairings) successfully merged out of 35793 (in 59 pairings) input.

    ## 40158 paired-reads (in 16 unique pairings) successfully merged out of 41373 (in 36 pairings) input.

    ## 44458 paired-reads (in 30 unique pairings) successfully merged out of 45059 (in 69 pairings) input.

    ## 43982 paired-reads (in 38 unique pairings) successfully merged out of 49269 (in 108 pairings) input.

    ## 26650 paired-reads (in 23 unique pairings) successfully merged out of 27055 (in 49 pairings) input.

    ## 30108 paired-reads (in 14 unique pairings) successfully merged out of 31426 (in 36 pairings) input.

    ## 36624 paired-reads (in 12 unique pairings) successfully merged out of 37095 (in 35 pairings) input.

    ## 38617 paired-reads (in 16 unique pairings) successfully merged out of 39643 (in 51 pairings) input.

    ## 40767 paired-reads (in 32 unique pairings) successfully merged out of 42985 (in 91 pairings) input.

    ## 40905 paired-reads (in 43 unique pairings) successfully merged out of 46537 (in 127 pairings) input.

    ## 25225 paired-reads (in 9 unique pairings) successfully merged out of 26678 (in 28 pairings) input.

    ## 31717 paired-reads (in 22 unique pairings) successfully merged out of 32499 (in 50 pairings) input.

    ## 38215 paired-reads (in 17 unique pairings) successfully merged out of 39167 (in 45 pairings) input.

    ## 42052 paired-reads (in 27 unique pairings) successfully merged out of 42936 (in 53 pairings) input.

    ## 33701 paired-reads (in 50 unique pairings) successfully merged out of 39287 (in 165 pairings) input.

    ## 39608 paired-reads (in 39 unique pairings) successfully merged out of 45945 (in 123 pairings) input.

    ## 22803 paired-reads (in 16 unique pairings) successfully merged out of 24568 (in 45 pairings) input.

    ## 31568 paired-reads (in 17 unique pairings) successfully merged out of 32321 (in 45 pairings) input.

    ## 41147 paired-reads (in 11 unique pairings) successfully merged out of 41677 (in 41 pairings) input.

    ## 40316 paired-reads (in 9 unique pairings) successfully merged out of 40834 (in 31 pairings) input.

    ## 41524 paired-reads (in 34 unique pairings) successfully merged out of 43896 (in 107 pairings) input.

    ## 39812 paired-reads (in 21 unique pairings) successfully merged out of 43330 (in 71 pairings) input.

    ## 25475 paired-reads (in 7 unique pairings) successfully merged out of 26125 (in 31 pairings) input.

    ## 24885 paired-reads (in 4 unique pairings) successfully merged out of 25290 (in 21 pairings) input.

    ## 44306 paired-reads (in 46 unique pairings) successfully merged out of 45003 (in 74 pairings) input.

    ## 40409 paired-reads (in 41 unique pairings) successfully merged out of 41312 (in 89 pairings) input.

    ## 45069 paired-reads (in 64 unique pairings) successfully merged out of 49757 (in 171 pairings) input.

    ## 46622 paired-reads (in 64 unique pairings) successfully merged out of 50778 (in 149 pairings) input.

    ## 27978 paired-reads (in 44 unique pairings) successfully merged out of 28620 (in 86 pairings) input.

    ## 30142 paired-reads (in 39 unique pairings) successfully merged out of 30732 (in 87 pairings) input.

    ## 47944 paired-reads (in 18 unique pairings) successfully merged out of 48680 (in 82 pairings) input.

    ## 43613 paired-reads (in 10 unique pairings) successfully merged out of 44036 (in 53 pairings) input.

    ## 40352 paired-reads (in 47 unique pairings) successfully merged out of 50908 (in 195 pairings) input.

    ## 44501 paired-reads (in 39 unique pairings) successfully merged out of 49558 (in 146 pairings) input.

    ## 30357 paired-reads (in 15 unique pairings) successfully merged out of 31008 (in 80 pairings) input.

    ## 27595 paired-reads (in 9 unique pairings) successfully merged out of 28024 (in 43 pairings) input.

    ## 57388 paired-reads (in 11 unique pairings) successfully merged out of 58181 (in 31 pairings) input.

    ## 45389 paired-reads (in 3 unique pairings) successfully merged out of 45561 (in 9 pairings) input.

    ## 47731 paired-reads (in 33 unique pairings) successfully merged out of 52983 (in 109 pairings) input.

    ## 45726 paired-reads (in 18 unique pairings) successfully merged out of 45959 (in 52 pairings) input.

    ## 36227 paired-reads (in 9 unique pairings) successfully merged out of 36615 (in 37 pairings) input.

    ## 26134 paired-reads (in 4 unique pairings) successfully merged out of 26205 (in 13 pairings) input.

    ## 35069 paired-reads (in 36 unique pairings) successfully merged out of 35893 (in 73 pairings) input.

    ## 38535 paired-reads (in 35 unique pairings) successfully merged out of 40049 (in 78 pairings) input.

    ## 31824 paired-reads (in 49 unique pairings) successfully merged out of 36043 (in 129 pairings) input.

    ## 33003 paired-reads (in 41 unique pairings) successfully merged out of 40179 (in 119 pairings) input.

    ## 21476 paired-reads (in 32 unique pairings) successfully merged out of 22218 (in 70 pairings) input.

    ## 22777 paired-reads (in 30 unique pairings) successfully merged out of 23487 (in 73 pairings) input.

    ## 26187 paired-reads (in 8 unique pairings) successfully merged out of 26550 (in 28 pairings) input.

    ## 34564 paired-reads (in 11 unique pairings) successfully merged out of 34936 (in 42 pairings) input.

    ## 24890 paired-reads (in 48 unique pairings) successfully merged out of 29633 (in 134 pairings) input.

    ## 33970 paired-reads (in 42 unique pairings) successfully merged out of 43504 (in 136 pairings) input.

    ## 21793 paired-reads (in 6 unique pairings) successfully merged out of 22458 (in 28 pairings) input.

    ## 28719 paired-reads (in 10 unique pairings) successfully merged out of 29081 (in 39 pairings) input.

    ## 49153 paired-reads (in 10 unique pairings) successfully merged out of 49537 (in 37 pairings) input.

    ## 34149 paired-reads (in 7 unique pairings) successfully merged out of 34610 (in 29 pairings) input.

    ## 32923 paired-reads (in 47 unique pairings) successfully merged out of 46631 (in 147 pairings) input.

    ## 29707 paired-reads (in 49 unique pairings) successfully merged out of 35323 (in 136 pairings) input.

    ## 31075 paired-reads (in 5 unique pairings) successfully merged out of 35908 (in 34 pairings) input.

    ## 21916 paired-reads (in 7 unique pairings) successfully merged out of 22633 (in 32 pairings) input.

    ## 34006 paired-reads (in 16 unique pairings) successfully merged out of 34625 (in 47 pairings) input.

    ## 47385 paired-reads (in 14 unique pairings) successfully merged out of 48041 (in 42 pairings) input.

    ## 37241 paired-reads (in 19 unique pairings) successfully merged out of 39108 (in 67 pairings) input.

    ## 49670 paired-reads (in 28 unique pairings) successfully merged out of 53467 (in 94 pairings) input.

    ## 20854 paired-reads (in 11 unique pairings) successfully merged out of 21308 (in 34 pairings) input.

    ## 34117 paired-reads (in 15 unique pairings) successfully merged out of 35051 (in 40 pairings) input.

    ## 45955 paired-reads (in 16 unique pairings) successfully merged out of 50216 (in 78 pairings) input.

    ## 44486 paired-reads (in 9 unique pairings) successfully merged out of 46251 (in 60 pairings) input.

    ## 40944 paired-reads (in 120 unique pairings) successfully merged out of 54389 (in 381 pairings) input.

    ## 35237 paired-reads (in 60 unique pairings) successfully merged out of 52191 (in 207 pairings) input.

    ## 35625 paired-reads (in 14 unique pairings) successfully merged out of 41114 (in 87 pairings) input.

    ## 30861 paired-reads (in 15 unique pairings) successfully merged out of 38964 (in 77 pairings) input.

    ## 41999 paired-reads (in 31 unique pairings) successfully merged out of 43068 (in 67 pairings) input.

    ## 43234 paired-reads (in 28 unique pairings) successfully merged out of 44529 (in 64 pairings) input.

    ## 38636 paired-reads (in 63 unique pairings) successfully merged out of 41814 (in 144 pairings) input.

    ## 49293 paired-reads (in 50 unique pairings) successfully merged out of 55934 (in 148 pairings) input.

    ## 31483 paired-reads (in 30 unique pairings) successfully merged out of 32228 (in 60 pairings) input.

    ## 35573 paired-reads (in 29 unique pairings) successfully merged out of 36578 (in 61 pairings) input.

    ## 35929 paired-reads (in 18 unique pairings) successfully merged out of 37745 (in 44 pairings) input.

    ## 45130 paired-reads (in 18 unique pairings) successfully merged out of 45933 (in 52 pairings) input.

    ## 35455 paired-reads (in 43 unique pairings) successfully merged out of 39999 (in 148 pairings) input.

    ## 40056 paired-reads (in 76 unique pairings) successfully merged out of 55105 (in 239 pairings) input.

    ## 26585 paired-reads (in 13 unique pairings) successfully merged out of 27380 (in 43 pairings) input.

    ## 35754 paired-reads (in 12 unique pairings) successfully merged out of 42068 (in 70 pairings) input.

    ## 34575 paired-reads (in 23 unique pairings) successfully merged out of 35610 (in 56 pairings) input.

    ## 38530 paired-reads (in 32 unique pairings) successfully merged out of 39429 (in 69 pairings) input.

    ## 32959 paired-reads (in 31 unique pairings) successfully merged out of 34515 (in 79 pairings) input.

    ## 40182 paired-reads (in 37 unique pairings) successfully merged out of 41935 (in 89 pairings) input.

    ## 22020 paired-reads (in 18 unique pairings) successfully merged out of 22509 (in 42 pairings) input.

    ## 24683 paired-reads (in 33 unique pairings) successfully merged out of 25267 (in 73 pairings) input.

    ## 37884 paired-reads (in 6 unique pairings) successfully merged out of 38579 (in 21 pairings) input.

    ## 33348 paired-reads (in 9 unique pairings) successfully merged out of 35099 (in 22 pairings) input.

    ## 45807 paired-reads (in 31 unique pairings) successfully merged out of 49888 (in 81 pairings) input.

    ## 31024 paired-reads (in 54 unique pairings) successfully merged out of 39971 (in 141 pairings) input.

    ## 29426 paired-reads (in 5 unique pairings) successfully merged out of 30911 (in 23 pairings) input.

    ## 26473 paired-reads (in 9 unique pairings) successfully merged out of 26552 (in 24 pairings) input.

    ## 32390 paired-reads (in 4 unique pairings) successfully merged out of 33702 (in 31 pairings) input.

    ## 43634 paired-reads (in 7 unique pairings) successfully merged out of 44272 (in 31 pairings) input.

    ## 28084 paired-reads (in 39 unique pairings) successfully merged out of 37102 (in 117 pairings) input.

    ## 36604 paired-reads (in 34 unique pairings) successfully merged out of 43280 (in 107 pairings) input.

    ## 23473 paired-reads (in 11 unique pairings) successfully merged out of 25400 (in 42 pairings) input.

    ## 24606 paired-reads (in 5 unique pairings) successfully merged out of 25449 (in 26 pairings) input.

    ## 37910 paired-reads (in 16 unique pairings) successfully merged out of 38637 (in 42 pairings) input.

    ## 42110 paired-reads (in 15 unique pairings) successfully merged out of 42854 (in 35 pairings) input.

    ## 43409 paired-reads (in 32 unique pairings) successfully merged out of 46761 (in 99 pairings) input.

    ## 42330 paired-reads (in 30 unique pairings) successfully merged out of 46211 (in 90 pairings) input.

    ## 27483 paired-reads (in 10 unique pairings) successfully merged out of 27882 (in 25 pairings) input.

    ## 27837 paired-reads (in 11 unique pairings) successfully merged out of 28570 (in 32 pairings) input.

    ## 37964 paired-reads (in 9 unique pairings) successfully merged out of 38217 (in 29 pairings) input.

    ## 41989 paired-reads (in 9 unique pairings) successfully merged out of 42313 (in 36 pairings) input.

    ## 41742 paired-reads (in 24 unique pairings) successfully merged out of 43682 (in 62 pairings) input.

    ## 38303 paired-reads (in 37 unique pairings) successfully merged out of 47097 (in 117 pairings) input.

    ## 24352 paired-reads (in 8 unique pairings) successfully merged out of 24577 (in 26 pairings) input.

    ## 27742 paired-reads (in 8 unique pairings) successfully merged out of 28106 (in 40 pairings) input.

    ## 7448 paired-reads (in 4 unique pairings) successfully merged out of 7703 (in 10 pairings) input.

    ## 38699 paired-reads (in 11 unique pairings) successfully merged out of 39677 (in 31 pairings) input.

    ## 7328 paired-reads (in 10 unique pairings) successfully merged out of 7560 (in 25 pairings) input.

    ## 40317 paired-reads (in 29 unique pairings) successfully merged out of 44434 (in 77 pairings) input.

    ## 5922 paired-reads (in 4 unique pairings) successfully merged out of 6065 (in 8 pairings) input.

    ## 26569 paired-reads (in 9 unique pairings) successfully merged out of 26968 (in 19 pairings) input.

    ## 26050 paired-reads (in 12 unique pairings) successfully merged out of 26334 (in 28 pairings) input.

    ## 25425 paired-reads (in 16 unique pairings) successfully merged out of 25962 (in 37 pairings) input.

    ## 30374 paired-reads (in 30 unique pairings) successfully merged out of 31843 (in 73 pairings) input.

    ## 26625 paired-reads (in 35 unique pairings) successfully merged out of 29968 (in 112 pairings) input.

    ## 20193 paired-reads (in 11 unique pairings) successfully merged out of 20336 (in 26 pairings) input.

    ## 20321 paired-reads (in 16 unique pairings) successfully merged out of 20653 (in 34 pairings) input.

    ## 1770 paired-reads (in 4 unique pairings) successfully merged out of 28521 (in 33 pairings) input.

    ## 11876 paired-reads (in 72 unique pairings) successfully merged out of 38146 (in 174 pairings) input.

    ## 12012 paired-reads (in 16 unique pairings) successfully merged out of 47787 (in 73 pairings) input.

    ## 4566 paired-reads (in 5 unique pairings) successfully merged out of 26994 (in 16 pairings) input.

    ## 9971 paired-reads (in 59 unique pairings) successfully merged out of 33400 (in 173 pairings) input.

    ## 604 paired-reads (in 2 unique pairings) successfully merged out of 30324 (in 18 pairings) input.

    ## 1429 paired-reads (in 4 unique pairings) successfully merged out of 21381 (in 25 pairings) input.

    ## 8745 paired-reads (in 51 unique pairings) successfully merged out of 31742 (in 138 pairings) input.

    ## 1890 paired-reads (in 4 unique pairings) successfully merged out of 27163 (in 10 pairings) input.

    ## 914 paired-reads (in 2 unique pairings) successfully merged out of 30546 (in 17 pairings) input.

    ## 9611 paired-reads (in 54 unique pairings) successfully merged out of 34945 (in 151 pairings) input.

    ## 4435 paired-reads (in 6 unique pairings) successfully merged out of 34983 (in 43 pairings) input.

    ## 16891 paired-reads (in 17 unique pairings) successfully merged out of 18355 (in 42 pairings) input.

    ## 10117 paired-reads (in 61 unique pairings) successfully merged out of 30902 (in 176 pairings) input.

    ## 5665 paired-reads (in 7 unique pairings) successfully merged out of 29542 (in 29 pairings) input.

    ## 14062 paired-reads (in 8 unique pairings) successfully merged out of 27196 (in 24 pairings) input.

    ## 10753 paired-reads (in 71 unique pairings) successfully merged out of 38365 (in 192 pairings) input.

    ## 5510 paired-reads (in 5 unique pairings) successfully merged out of 39566 (in 15 pairings) input.

    ## 5227 paired-reads (in 4 unique pairings) successfully merged out of 28440 (in 22 pairings) input.

    ## 8953 paired-reads (in 58 unique pairings) successfully merged out of 34463 (in 175 pairings) input.

    ## 6066 paired-reads (in 7 unique pairings) successfully merged out of 41038 (in 32 pairings) input.

    ## 4279 paired-reads (in 9 unique pairings) successfully merged out of 26098 (in 32 pairings) input.

    ## 11402 paired-reads (in 52 unique pairings) successfully merged out of 31308 (in 160 pairings) input.

    ## 12421 paired-reads (in 14 unique pairings) successfully merged out of 27998 (in 49 pairings) input.

    ## 6243 paired-reads (in 4 unique pairings) successfully merged out of 24944 (in 16 pairings) input.

    ## 9774 paired-reads (in 50 unique pairings) successfully merged out of 33691 (in 135 pairings) input.

    ## 3657 paired-reads (in 6 unique pairings) successfully merged out of 35577 (in 24 pairings) input.

    ## 25923 paired-reads (in 5 unique pairings) successfully merged out of 26382 (in 18 pairings) input.

    ## 29980 paired-reads (in 11 unique pairings) successfully merged out of 30814 (in 40 pairings) input.

    ## 19337 paired-reads (in 4 unique pairings) successfully merged out of 19623 (in 13 pairings) input.

``` r
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                      sequence
    ## 1    CACCGCGGTCATACGATTAACCCAAACTAACGGGCCTACGGCGTAAAGCGTGTAAAAGATCTATTAATACTAAAGTTAAAATTTAACCAAGCCGTAAAAAGCTACCGTTAATATAAAATAAGCTACGAAAGTGACTTTATTAACTCTGATTACACGATAGCTAAGACC
    ## 2 CACCGCGGTTATACGAGAGGCCCAAATTGATGAAAAACGGCGTAAAGCGTGGTTAAGAAAAAAAAGAGAAAATATGGCCGAACAGCTTCAAAGCAGTTATACGCATCCGAAGTCACGAAGAACAATCACGAAAGTTGCCCTAAAACCTCCGATTCCACGAAAGCCATAAAA
    ## 3   CACCGCGGTTATACGAGAGGCCCAAGTTGACAACCACCGGCGTAAAGAGTGGTTATGGAAACACTATACTAAAGCCGAACACCCTCTAGGCTGTCATACGCACCTGAGGCTACGAAGCTCCCCCACGAAAGTGGCTTTATTCCTCCCTGAACCCACGACAGCTACGATA
    ## 5   CACCGCGGTTATACGAGAGGCCCAAGTTGACAACCACCGGCGTAAAGAGTGGTTATGGAAACACTATACTAAAGCCGAACACCCTCTAGGCTGTCATACGCACCTGAGACTACGAAGCTCCCCCACGAAAGTGGCTTTATTCCTCCCTGAACCCACGACAGCTACGATA
    ## 6    CACCGCGGTCATACGATTAACCCAAACTAACGGGCCTACGGCGTAAAGCGTGTAAAAGATCTATTAATACTAAAGTTAAAATTTAACCAAGCCGTAAAAAGCTACCGTTAATATAAAATAAGCTACGAAAGTGACTTTATTAACTCTGATTACACGATAGCTAAGACT
    ## 7   CACCGCAGTCATTTGATTAACCCAAATTAATAGGCCTGTGGTGTAAAACGTGTTAAAGATATTTCCATACAAAAGTTAAGACTTAACTAAGCCGTAAGAAGCTACAATGACCACAAAAAGAAGCTACGGGAGTGACTTCAATCCTTCTGACTACACAATAACTAAAACC
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1     31055       1       1     32         0      0      1   TRUE
    ## 2      5452       2       2     29         0      0      2   TRUE
    ## 3       654       3       3     31         0      0      1   TRUE
    ## 5       318       3       4     31         0      0      1   TRUE
    ## 6       216       1       9     32         0      0      1   TRUE
    ## 7       144       4       5     31         0      0      2   TRUE

Make a sequence table

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  288 1638

``` r
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 113 116 121 124 128 132 136 138 140 147 148 149 152 156 158 159 161 163 164 165 
    ##   1   4   1   1   2   1   4   1   1   3   1   3   2   2   1   1   1   2   1   2 
    ## 166 167 168 169 170 171 172 173 174 175 177 178 179 180 181 182 183 186 187 188 
    ##   5  20 320 808 271  30  59   8   6   1   3   2   6  29  12   1   1   1   1  13 
    ## 189 190 
    ##   1   5

Let’s remove the singletons and off-target sequences

``` r
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 162:182] # expanded the range based on Kim's data analysis
```

Remove chimeras

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 713 bimeras out of 1586 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1] 288 873

Calculate frequency of chimeras

``` r
sum(seqtab.nochim)/sum(seqtab2)
```

    ## [1] 0.9926118

Track reads through the pipeline

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
```

    ##                  input filtered denoisedF denoisedR merged nonchim
    ## e03719-plate1-a  39155    39017     38584     38712  37891   37839
    ## e03719-plate1-b  37077    36937     36409     36363  34387   34385
    ## e03719-plate2-a  39338    38973     37354     37485  30862   30731
    ## e03719-plate2-b  47160    46766     44713     45053  36210   36016
    ## e03719-plate3-a  26731    26614     26132     26255  24976   24944
    ## e03719-plate3-b  32975    32782     31961     31919  28492   28486
    ## e03720-plate1-a  40942    40836     40609     40564  39847   39828
    ## e03720-plate1-b  39602    39476     39034     39099  37178   37135
    ## e03720-plate2-a  44603    44111     42040     42510  33461   33425
    ## e03720-plate2-b  44825    44482     43072     43307  37222   37066
    ## e03720-plate3-a  25528    25363     24725     24823  22354   22340
    ## e03720-plate3-b  29589    29479     29059     28989  27551   27545
    ## e03721-plate1-a  39337    39151     38585     37385  36573   36437
    ## e03721-plate1-b  53415    53206     52397     52496  49601   49364
    ## e03721-plate2-a  41972    41422     38942     38879  28449   28347
    ## e03721-plate2-b  62343    61733     58767     60075  45433   45077
    ## e03721-plate3-a  25918    25791     25377     25049  24161   24043
    ## e03721-plate3-b  41127    40901     39907     40280  36246   36106
    ## e03722-plate1-a  43498    43354     43153     42928  42341   42228
    ## e03722-plate1-b  52330    52145     51951     51508  50900   50731
    ## e03722-plate2-a  45339    45145     44772     44663  43408   43251
    ## e03722-plate2-b  57448    57210     56811     56722  55332   55187
    ## e03722-plate3-a  28697    28576     28376     28309  27855   27735
    ## e03722-plate3-b  36079    35955     35771     35633  35103   34998
    ## e03723-plate1-a  38558    38375     38112     37783  36844   35488
    ## e03723-plate1-b  49943    49693     49127     48862  47302   45392
    ## e03723-plate2-a  43963    43738     42964     41218  39827   38524
    ## e03723-plate2-b  58530    58206     57205     56497  51093   49181
    ## e03723-plate3-a  29478    29354     29139     28856  28175   27208
    ## e03723-plate3-b  37369    37192     36895     36241  35401   34140
    ## e03724-plate1-a  38833    38644     38355     38119  37230   37073
    ## e03724-plate1-b  41691    41439     41015     40655  39116   38896
    ## e03724-plate2-a  41544    41280     40458     40505  37214   37109
    ## e03724-plate2-b  49263    48845     47365     47363  40572   40260
    ## e03724-plate3-a  21798    21675     21480     21356  20793   20766
    ## e03724-plate3-b  30831    30671     30356     30230  29387   29267
    ## e03725-plate1-a  43271    43037     42745     42428  40635   40421
    ## e03725-plate1-b  43127    42905     42634     42364  40610   40439
    ## e03725-plate2-a  45477    45206     44452     44347  41032   40656
    ## e03725-plate2-b  45898    45682     45249     45074  42365   42013
    ## e03725-plate3-a  25780    25655     25421     25325  24185   24040
    ## e03725-plate3-b  28341    28199     28024     27837  26676   26586
    ## e03726-plate1-a  40444    40184     39683     39446  38052   37945
    ## e03726-plate1-b  43693    43380     42751     42563  39756   39604
    ## e03726-plate2-a  41935    41199     38035     39063  24287   24220
    ## e03726-plate2-b  46282    45616     42398     43601  27781   27316
    ## e03726-plate3-a  23454    23246     22684     22671  21017   20787
    ## e03726-plate3-b  33302    33058     32330     32396  29039   29003
    ## e03727-plate1-a  38789    38693     38562     38400  38084   38046
    ## e03727-plate1-b  38874    38741     38574     38384  37981   37897
    ## e03727-plate2-a  43095    42978     42639     42645  41472   41414
    ## e03727-plate2-b  55098    54868     54391     54052  52193   52089
    ## e03727-plate3-a  27775    27696     27550     27506  27145   27088
    ## e03727-plate3-b  31682    31577     31388     31117  30730   30610
    ## e03728-plate1-a  44570    44388     44017     43765  42606   42256
    ## e03728-plate1-b  41385    41192     40840     40686  39559   39344
    ## e03728-plate2-a  51364    50886     49004     49419  41354   41086
    ## e03728-plate2-b  45657    45257     43712     44044  37485   37306
    ## e03728-plate3-a  29563    29400     29028     28885  27766   27539
    ## e03728-plate3-b  28160    28052     27765     27649  26756   26661
    ## e03729-plate1-a  52240    52070     51814     51421  50793   50468
    ## e03729-plate1-b  42818    42583     42147     41825  40896   40817
    ## e03729-plate2-a  61821    61516     60488     60537  56225   55703
    ## e03729-plate2-b  48869    48410     47169     47327  42808   42506
    ## e03729-plate3-a  37693    37551     37287     36951  36466   36346
    ## e03729-plate3-b  25877    25714     25321     25169  24257   24188
    ## e03730-plate1-a  60690    60579     60406     59957  59722   59652
    ## e03730-plate1-b  46845    46753     46582     46064  45863   45820
    ## e03730-plate2-a  59834    59452     57840     57237  50816   50768
    ## e03730-plate2-b  49747    49504     48621     48170  44751   44696
    ## e03730-plate3-a  32588    32504     32329     31921  31628   31596
    ## e03730-plate3-b  31838    31783     31685     31449  31272   31237
    ## e03731-plate1-a  53900    53713     53297     53078  52249   52072
    ## e03731-plate1-b  53966    53867     53748     53454  53211   53211
    ## e03731-plate2-a  59266    58977     58170     55889  55185   55109
    ## e03731-plate2-b  62424    62290     62002     61875  60907   60907
    ## e03731-plate3-a  33930    33817     33522     33250  32775   32705
    ## e03731-plate3-b  39892    39820     39709     39560  39391   39391
    ## e03732-plate1-a  51003    50850     50420     50456  49184   48666
    ## e03732-plate1-b  56624    56367     55825     55720  53627   53597
    ## e03732-plate2-a  52835    52403     50542     51220  42132   41766
    ## e03732-plate2-b  58807    58087     55419     56462  43522   43078
    ## e03732-plate3-a  33619    33499     33049     33121  31172   31172
    ## e03732-plate3-b  35204    35017     34410     34445  32180   32161
    ## e03733-plate1-a  23604    23523     23394     23374  23135   23135
    ## e03733-plate1-b  41111    41036     40922     40849  40693   40660
    ## e03733-plate2-a  21862    21703     21161     21290  19554   19548
    ## e03733-plate2-b  38683    38492     37855     37950  35439   35422
    ## e03733-plate3-a  11918    11867     11715     11722  11358   11345
    ## e03733-plate3-b  24118    24073     23914     23848  23398   23398
    ## e03734-plate1-a  35767    35678     35491     35448  35148   35061
    ## e03734-plate1-b  40469    40378     40214     39976  39674   39668
    ## e03734-plate2-a  37730    37525     36752     36951  33445   33374
    ## e03734-plate2-b  51250    51034     50248     50121  47225   47221
    ## e03734-plate3-a  23214    23132     22955     22918  22483   22483
    ## e03734-plate3-b  32235    32169     32033     31989  31720   31720
    ## e03735-plate1-a  44355    44286     44167     44119  43865   43865
    ## e03735-plate1-b  42743    42665     42524     42471  42014   42014
    ## e03735-plate2-a  45201    45018     44351     44557  41527   41487
    ## e03735-plate2-b  52261    52002     50702     51140  44878   44696
    ## e03735-plate3-a  27934    27881     27695     27598  26962   26962
    ## e03735-plate3-b  34118    34061     33938     33919  33739   33681
    ## e03736-plate1-a  41517    41451     41368     41325  41207   41207
    ## e03736-plate1-b  37295    37127     36807     36519  35554   35484
    ## e03736-plate2-a  41015    40735     39697     39829  35151   35132
    ## e03736-plate2-b  41322    40887     39454     39656  33261   33145
    ## e03736-plate3-a  24563    24513     24340     24335  23821   23821
    ## e03736-plate3-b  27486    27344     26926     26909  25656   25578
    ## e03737-plate1-a  44874    44606     43975     43712  42374   41196
    ## e03737-plate1-b  43002    42728     41980     41907  39667   38872
    ## e03737-plate2-a  50460    49920     48135     47949  40859   39799
    ## e03737-plate2-b  50131    49529     47363     47656  39111   38384
    ## e03737-plate3-a  28341    28161     27540     27608  26011   25391
    ## e03737-plate3-b  30564    30352     29684     29678  27404   27010
    ## e03738-plate1-a  59787    59561     58812     58608  54099   53659
    ## e03738-plate1-b  44578    44247     43824     42055  38847   38594
    ## e03738-plate2-a  55813    55172     52532     52564  39283   38670
    ## e03738-plate2-b  44207    43140     40875     38460  27000   26916
    ## e03738-plate3-a  33076    32860     32034     29455  28052   27747
    ## e03738-plate3-b  30410    29931     29096     26749  24097   23971
    ## e03739-plate1-a  36966    36762     36394     36182  35629   34255
    ## e03739-plate1-b  36781    36650     36143     35985  34819   34416
    ## e03739-plate2-a  41137    40851     40081     38928  38296   36688
    ## e03739-plate2-b  43605    43349     42447     41963  38867   38407
    ## e03739-plate3-a  24930    24801     24390     24344  23883   23038
    ## e03739-plate3-b  27125    27020     26625     26538  25763   25463
    ## e03740-plate1-a  54175    53712     53301     52416  51593   50927
    ## e03740-plate1-b  40957    40677     40211     39960  38692   38487
    ## e03740-plate2-a  55352    54761     53295     52455  46962   46378
    ## e03740-plate2-b  43629    43149     41324     41737  33619   33412
    ## e03740-plate3-a  32764    32504     32083     31453  29920   29545
    ## e03740-plate3-b  28296    28099     27630     27631  26040   25944
    ## e03741-plate1-a  36559    36437     36177     35957  35190   34633
    ## e03741-plate1-b  42426    42262     41871     41586  40158   39862
    ## e03741-plate2-a  46786    46617     46151     45258  44458   43826
    ## e03741-plate2-b  52077    51694     50250     50072  43982   43437
    ## e03741-plate3-a  27657    27563     27339     27184  26650   26342
    ## e03741-plate3-b  32427    32329     31886     31694  30108   29869
    ## e03742-plate1-a  38414    38252     37734     37351  36624   36276
    ## e03742-plate1-b  40975    40772     40215     39922  38617   38224
    ## e03742-plate2-a  44891    44615     43763     43451  40767   40281
    ## e03742-plate2-b  50054    49632     47638     47456  40905   40156
    ## e03742-plate3-a  27902    27776     27274     27102  25225   24996
    ## e03742-plate3-b  33657    33470     32942     32761  31717   31302
    ## e03743-plate1-a  41842    41351     40041     39874  38215   38077
    ## e03743-plate1-b  45226    44775     43796     43499  42052   41872
    ## e03743-plate2-a  43240    42641     40155     41029  33701   33371
    ## e03743-plate2-b  50715    50089     47487     48045  39608   39306
    ## e03743-plate3-a  26846    26550     25370     25616  22803   22733
    ## e03743-plate3-b  33931    33633     32788     32827  31568   31451
    ## e03744-plate1-a  42163    42068     41868     41804  41147   40943
    ## e03744-plate1-b  41331    41221     40995     40969  40316   40104
    ## e03744-plate2-a  45065    44877     44192     44375  41524   41195
    ## e03744-plate2-b  44913    44685     43711     43938  39812   39452
    ## e03744-plate3-a  26538    26474     26244     26289  25475   25448
    ## e03744-plate3-b  25709    25647     25441     25418  24885   24885
    ## e03745-plate1-a  46339    46097     45765     45177  44306   43277
    ## e03745-plate1-b  43012    42779     42225     41673  40409   39746
    ## e03745-plate2-a  52610    52226     50811     50628  45069   44206
    ## e03745-plate2-b  54128    53694     52231     51663  46622   45939
    ## e03745-plate3-a  29980    29822     29407     28867  27978   27230
    ## e03745-plate3-b  32076    31911     31524     30964  30142   29588
    ## e03746-plate1-a  49697    49539     49169     48906  47944   47714
    ## e03746-plate1-b  44816    44674     44412     44177  43613   43455
    ## e03746-plate2-a  55617    54966     52304     52497  40352   40041
    ## e03746-plate2-b  54954    54402     52360     50645  44501   44175
    ## e03746-plate3-a  34109    33900     33132     31252  30357   30207
    ## e03746-plate3-b  29586    29446     28960     28242  27595   27477
    ## e03747-plate1-a  58714    58603     58361     58306  57388   57361
    ## e03747-plate1-b  46220    46100     45967     45607  45389   45389
    ## e03747-plate2-a  55149    54824     53517     53775  47731   47695
    ## e03747-plate2-b  49168    48925     48250     46096  45726   45726
    ## e03747-plate3-a  37008    36922     36743     36719  36227   36145
    ## e03747-plate3-b  26732    26665     26512     26251  26134   26130
    ## e03748-plate1-a  37876    37468     37048     36039  35069   34227
    ## e03748-plate1-b  42247    41745     41044     40341  38535   37755
    ## e03748-plate2-a  39382    38882     37354     36676  31824   31084
    ## e03748-plate2-b  44709    44012     41572     41335  33003   32307
    ## e03748-plate3-a  23341    23112     22765     22370  21476   20937
    ## e03748-plate3-b  24632    24411     24044     23639  22777   22327
    ## e03749-plate1-a  27126    27011     26734     26728  26187   26187
    ## e03749-plate1-b  37561    37356     36791     35100  34564   34534
    ## e03749-plate2-a  31584    31307     30099     30393  24890   24755
    ## e03749-plate2-b  47331    46768     44436     44985  33970   33882
    ## e03749-plate3-a  22984    22896     22636     22606  21793   21793
    ## e03749-plate3-b  30900    30751     30334     29201  28719   28667
    ## e03750-plate1-a  50006    49898     49714     49636  49153   49122
    ## e03750-plate1-b  35813    35695     35328     34790  34149   34149
    ## e03750-plate2-a  52123    51323     47999     48581  32923   32662
    ## e03750-plate2-b  38556    38137     36347     36314  29707   29543
    ## e03750-plate3-a  37946    37650     36360     36690  31075   31075
    ## e03750-plate3-b  23206    23120     22827     22797  21916   21916
    ## e03751-plate1-a  35236    35111     34844     34802  34006   33878
    ## e03751-plate1-b  48756    48600     48338     48212  47385   47257
    ## e03751-plate2-a  40163    39973     39392     39515  37241   36990
    ## e03751-plate2-b  56052    55696     54654     54150  49670   49233
    ## e03751-plate3-a  21781    21682     21449     21445  20854   20843
    ## e03751-plate3-b  35863    35741     35402     35279  34117   34042
    ## e03752-plate1-a  52075    51821     50766     50962  45955   45431
    ## e03752-plate1-b  47459    47264     46644     46660  44486   43934
    ## e03752-plate2-a  61243    60464     57218     56312  40944   38883
    ## e03752-plate2-b  58655    57728     53745     54616  35237   34803
    ## e03752-plate3-a  43422    43145     41815     41907  35625   35129
    ## e03752-plate3-b  41715    41317     39680     39960  30861   30447
    ## e03753-plate1-a  44607    44413     44041     43240  41999   41002
    ## e03753-plate1-b  46092    45887     45239     44938  43234   42339
    ## e03753-plate2-a  47701    47360     45895     42625  38636   37815
    ## e03753-plate2-b  60046    59482     57313     57211  49293   48411
    ## e03753-plate3-a  33511    33395     33071     32364  31483   30846
    ## e03753-plate3-b  38136    37951     37348     36957  35573   34834
    ## e03754-plate1-a  38737    38604     38059     38105  35929   35588
    ## e03754-plate1-b  47256    47020     46434     46306  45130   44202
    ## e03754-plate2-a  41741    41506     40396     40769  35455   34928
    ## e03754-plate2-b  62318    61432     57666     56896  40056   39203
    ## e03754-plate3-a  27923    27833     27534     27567  26585   26327
    ## e03754-plate3-b  45330    45004     43158     43053  35754   35561
    ## e03755-plate1-a  36637    36429     36009     35870  34575   33878
    ## e03755-plate1-b  40511    40332     39854     39733  38530   37766
    ## e03755-plate2-a  35746    35524     34895     34902  32959   32392
    ## e03755-plate2-b  43352    43113     42417     42345  40182   39483
    ## e03755-plate3-a  23227    23098     22776     22696  22020   21596
    ## e03755-plate3-b  26013    25886     25579     25445  24683   24043
    ## e03756-plate1-a  39171    39074     38827     38714  37884   37874
    ## e03756-plate1-b  37370    37210     36495     35399  33348   33341
    ## e03756-plate2-a  51719    51425     50301     50561  45807   45741
    ## e03756-plate2-b  43449    42973     40724     41374  31024   30938
    ## e03756-plate3-a  31813    31676     31124     31231  29426   29419
    ## e03756-plate3-b  26955    26898     26766     26594  26473   26457
    ## e03757-plate1-a  34281    34192     33825     33950  32390   32390
    ## e03757-plate1-b  44947    44782     44487     44430  43634   43634
    ## e03757-plate2-a  40081    39721     37693     38438  28084   27930
    ## e03757-plate2-b  46223    45746     43957     44398  36604   36493
    ## e03757-plate3-a  26133    26014     25551     25724  23473   23473
    ## e03757-plate3-b  26099    25964     25621     25631  24606   24606
    ## e03758-plate1-a  39644    39413     39090     38856  37910   37732
    ## e03758-plate1-b  44002    43770     43362     43093  42110   41924
    ## e03758-plate2-a  48862    48502     47442     47382  43409   43186
    ## e03758-plate2-b  48464    48092     46868     46989  42330   42225
    ## e03758-plate3-a  28881    28727     28372     28093  27483   27400
    ## e03758-plate3-b  29482    29293     28904     28814  27837   27751
    ## e03759-plate1-a  38662    38577     38381     38355  37964   37958
    ## e03759-plate1-b  43264    43137     42868     42446  41989   41930
    ## e03759-plate2-a  44677    44529     43956     44052  41742   41712
    ## e03759-plate2-b  50423    49960     47936     48384  38303   38270
    ## e03759-plate3-a  24931    24860     24726     24665  24352   24344
    ## e03759-plate3-b  29590    29488     29092     28259  27742   27727
    ## e03760-plate1-a   8006     7973      7844      7784   7448    7423
    ## e03760-plate1-b  40522    40416     40044     39930  38699   38407
    ## e03760-plate2-a   8180     8135      7908      7660   7328    7289
    ## e03760-plate2-b  46616    46338     45116     45120  40317   40014
    ## e03760-plate3-a   6278     6257      6134      6154   5922    5904
    ## e03760-plate3-b  27997    27915     27531     27199  26569   26400
    ## e03761-plate1-a  27012    26851     26649     26445  26050   25906
    ## e03761-plate1-b  27227    26987     26574     26157  25425   25202
    ## e03761-plate2-a  33169    32943     32291     32143  30374   30109
    ## e03761-plate2-b  32793    32363     30927     30580  26625   26058
    ## e03761-plate3-a  20812    20715     20552     20405  20193   20041
    ## e03761-plate3-b  21673    21466     21100     20820  20321   20159
    ## e03762-plate1-EB 34995    33791     30579     30598   1770    1770
    ## e03762-plate2-EB 48031    46535     41934     40112  11876   11679
    ## e03762-plate3-EB 57523    56133     50681     50872  12012   11839
    ## e03763-plate1-EB 32769    31747     27769     30038   4566    4566
    ## e03763-plate2-EB 38917    37944     34311     35930   9971    9257
    ## e03763-plate3-EB 34827    33812     31005     32304    604     604
    ## e03764-plate1-EB 24679    24122     22278     22695   1429    1429
    ## e03764-plate2-EB 40365    39312     35571     33726   8745    8514
    ## e03764-plate3-EB 34647    33689     30962     28452   1890    1890
    ## e03765-plate1-EB 35975    34797     31371     32978    914     914
    ## e03765-plate2-EB 41435    40345     35999     37911   9611    8364
    ## e03765-plate3-EB 39424    38556     35733     36917   4435    4435
    ## e03766-plate1-EB 20164    19919     19117     18670  16891   16803
    ## e03766-plate2-EB 35633    34790     31601     33198  10117    9807
    ## e03766-plate3-EB 33499    32849     30003     31613   5665    5665
    ## e03767-plate1-EB 33047    32228     28662     29480  14062   14055
    ## e03767-plate2-EB 46311    45131     40799     40884  10753   10670
    ## e03767-plate3-EB 47737    46620     42860     41537   5510    5510
    ## e03768-plate1-EB 33438    32535     29251     30524   5227    5227
    ## e03768-plate2-EB 40192    39235     35450     37062   8953    8794
    ## e03768-plate3-EB 46461    45445     41767     43774   6066    5153
    ## e03769-plate1-FB 33442    31907     27974     28317   4279    4279
    ## e03769-plate2-FB 39853    38041     33609     33744  11402   11308
    ## e03769-plate3-FB 33825    32718     29357     29934  12421   12252
    ## NC-plate1        29458    28596     25550     27172   6243    6243
    ## NC-plate2        38753    37839     34459     36042   9774    9357
    ## NC-plate3        40945    39764     36114     38401   3657    1235
    ## PC-plate1        27680    27422     27190     26455  25923   25923
    ## PC-plate2        32362    32051     31643     30968  29980   29980
    ## PC-plate3        20539    20355     20147     19691  19337   19337

## Export files for taxonomy and samples/ASVs

``` r
 #make fasta file with ASVs
    asv_seqs=colnames(seqtab.nochim)
    for(i in 1:length(asv_seqs))
    {
        write.table(paste(">ASV",i, sep=""),file="csv_outputs/NFSfecal_MiFish_ASV_seqtab_nochim.csv", append=TRUE, col.names = F, row.names = F, quote=F)
        write.table(paste(asv_seqs[i], sep=""),file="csv_outputs/NFSfecal_MiFish_ASV_seqtab_nochim.csv", append=TRUE, col.names = F, row.names = F, quote=F)
    }
```

That’s the input for the FASTA blastn search.

Goal: change ASV headers to numbered ASVs that correspond to those
output in the FASTA file.

``` r
# Make map between brief names and full sequences
briefToSeq <- colnames(seqtab.nochim)
names(briefToSeq) <- paste0("ASV", seq(ncol(seqtab.nochim))) # Seq1, Seq2, ...
# Make new sequence table with brief names
st.brief <- seqtab.nochim
colnames(st.brief) <- names(briefToSeq)

# export the pool seq table with brief names:
write.csv(st.brief, file="csv_outputs/NFSfecal_MiFish_ASVtable.csv")
```

Now move on to the blastn taxonomy search.
