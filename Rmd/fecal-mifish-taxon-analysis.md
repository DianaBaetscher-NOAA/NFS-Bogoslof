fecal-mifish-taxonomy-analysis
================
dsb
2024-09-25

Using input data from the Rmd `fecal-mifish-dada2.Rmd` and
blastn/taxonkit output from sedna and the eDNA VM.

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
library(tidyr)
library(stringr)
library(ggplot2)
```

``` r
taxonomy <-read.delim("blastn/NFSfecal_mifish_shorttaxlineage.txt", header = FALSE, na.strings=c(""," ","NA"))

head(taxonomy)
```

    ##     V1                            V2      V3  V4 V5 V6 V7  V8    V9   V10
    ## 1 ASV1 gi|294514819|ref|NC_008415.3| 100.000 168  0  0  1 168   325   492
    ## 2 ASV1  gi|1475685114|gb|MG916809.1|  99.405 168  1  0  1 168 11645 11812
    ## 3 ASV2 gi|294514819|ref|NC_008415.3|  99.405 168  1  0  1 168   325   492
    ## 4 ASV2  gi|1475685114|gb|MG916809.1|  98.810 168  2  0  1 168 11645 11812
    ## 5 ASV3  gi|906357379|dbj|AB969999.1| 100.000 169  0  0  1 169     1   169
    ## 6 ASV3  gi|985566935|dbj|LC091602.1| 100.000 166  0  0  1 166     1   166
    ##        V11 V12 V13    V14    V15
    ## 1 4.51e-80 311 N/A  34884  34884
    ## 2 2.10e-78 305 N/A  34884  34884
    ## 3 2.10e-78 305 N/A  34884  34884
    ## 4 9.76e-77 300 N/A  34884  34884
    ## 5 1.26e-80 313 N/A 557346 557346
    ## 6 5.88e-79 307 N/A 557346 557346
    ##                                                                                                                                                                                                                                                                                                                                               V16
    ## 1     cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Teleostomi;Euteleostomi;Sarcopterygii;Dipnotetrapodomorpha;Tetrapoda;Amniota;Mammalia;Theria;Eutheria;Boreoeutheria;Laurasiatheria;Carnivora;Caniformia;Pinnipedia;Otariidae;Callorhinus;Callorhinus ursinus
    ## 2     cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Teleostomi;Euteleostomi;Sarcopterygii;Dipnotetrapodomorpha;Tetrapoda;Amniota;Mammalia;Theria;Eutheria;Boreoeutheria;Laurasiatheria;Carnivora;Caniformia;Pinnipedia;Otariidae;Callorhinus;Callorhinus ursinus
    ## 3     cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Teleostomi;Euteleostomi;Sarcopterygii;Dipnotetrapodomorpha;Tetrapoda;Amniota;Mammalia;Theria;Eutheria;Boreoeutheria;Laurasiatheria;Carnivora;Caniformia;Pinnipedia;Otariidae;Callorhinus;Callorhinus ursinus
    ## 4     cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Teleostomi;Euteleostomi;Sarcopterygii;Dipnotetrapodomorpha;Tetrapoda;Amniota;Mammalia;Theria;Eutheria;Boreoeutheria;Laurasiatheria;Carnivora;Caniformia;Pinnipedia;Otariidae;Callorhinus;Callorhinus ursinus
    ## 5 cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Teleostomi;Euteleostomi;Actinopterygii;Actinopteri;Neopterygii;Teleostei;Osteoglossocephalai;Clupeocephala;Euteleosteomorpha;Protacanthopterygii;Argentiniformes;Bathylagidae;Leuroglossus;Leuroglossus schmidti
    ## 6 cellular organisms;Eukaryota;Opisthokonta;Metazoa;Eumetazoa;Bilateria;Deuterostomia;Chordata;Craniata;Vertebrata;Gnathostomata;Teleostomi;Euteleostomi;Actinopterygii;Actinopteri;Neopterygii;Teleostei;Osteoglossocephalai;Clupeocephala;Euteleosteomorpha;Protacanthopterygii;Argentiniformes;Bathylagidae;Leuroglossus;Leuroglossus schmidti
    ##                                                                                              V17
    ## 1                Eukaryota;Chordata;Mammalia;Carnivora;Otariidae;Callorhinus;Callorhinus ursinus
    ## 2                Eukaryota;Chordata;Mammalia;Carnivora;Otariidae;Callorhinus;Callorhinus ursinus
    ## 3                Eukaryota;Chordata;Mammalia;Carnivora;Otariidae;Callorhinus;Callorhinus ursinus
    ## 4                Eukaryota;Chordata;Mammalia;Carnivora;Otariidae;Callorhinus;Callorhinus ursinus
    ## 5 Eukaryota;Chordata;Actinopteri;Argentiniformes;Bathylagidae;Leuroglossus;Leuroglossus schmidti
    ## 6 Eukaryota;Chordata;Actinopteri;Argentiniformes;Bathylagidae;Leuroglossus;Leuroglossus schmidti

# clean up the header a bit

``` r
# use the full taxonomy rather than the seq id to collapse identical entries
tax_df <- taxonomy %>%
  filter(V4 > 100) %>% # make sure all retained matches are >100 bp
  select(-V2, -(V5:V16)) %>%  #remove unnecessary columns
  group_by(V1, V17) %>% # group by the sequence key and the full taxonomy to reduce duplicate entries
  unique() %>% # doing that reduced the number of entries from 146k to 17k
  rename(qseqid=V1, perc_id=V3, length=V4, taxonomy=V17) %>% #rename headers
  filter(!str_detect(taxonomy, "environmental")) %>% # filter out any environmental samples
  filter(!str_detect(taxonomy, "synthetic")) %>% # filter out any synthetic "samples"
  filter(!str_detect(taxonomy, "bacterium"))
  #filter(perc_id >= 98) # seems like some of the matches below 98% are dubious (jellyfish and herring <1% different??)
```

``` r
tax_df %>% 
  ungroup() %>%
  select(qseqid) %>%
  unique()
```

    ## # A tibble: 695 × 1
    ##    qseqid
    ##    <chr> 
    ##  1 ASV1  
    ##  2 ASV2  
    ##  3 ASV3  
    ##  4 ASV4  
    ##  5 ASV5  
    ##  6 ASV6  
    ##  7 ASV7  
    ##  8 ASV8  
    ##  9 ASV9  
    ## 10 ASV10 
    ## # ℹ 685 more rows

### taxonomy clean-up

``` r
# formatting the taxonomy variables
taxon_df <- tax_df %>%
  filter(str_detect(taxonomy, ";")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")
```

Manually dealing with messy hybrid data:

``` r
fixed_tax_df <- taxon_df %>%
  filter(!str_detect(species, " x ")) %>% # remove hybrids
  filter(!species %in% c("Homo sapiens", "Canis lupus", "Oncorhynchus kawamurae", "Ammodytes japonicus", "Sebastes baramenuke", "Sebastes cheni", "Sebastes inermis")) %>% # remove species outside of their range
  filter(!str_detect(species, " sp. ")) %>%
  filter(!genus %in% c("Sus", "Bos", "Gorilla")) %>% # remove non-marine mammals
  filter(!order %in% c("Primates"))
```

There are four categories: 1. sequences that match a single species
unambiguously (the minority)

Sequences that match multiple species are divided in three categories:
2. top matches \> 2% identity than second-ranked matches 3. top matches
\< 2% identity than second-ranked matches 4. Multiple top matches with
the same % identity

``` r
# 1. sequences that are unambiguously a single species
single_spp_seqs <- fixed_tax_df %>% 
  group_by(qseqid) %>%
  add_tally(name = "n_taxa") %>%
  filter(n_taxa == 1)
  
single_spp_seqs
```

    ## # A tibble: 125 × 11
    ## # Groups:   qseqid [125]
    ##    qseqid perc_id length kingdom  phylum class order family genus species n_taxa
    ##    <chr>    <dbl>  <int> <chr>    <chr>  <chr> <chr> <chr>  <chr> <chr>    <int>
    ##  1 ASV9     100      169 Eukaryo… Chord… Acti… "Per… Hexag… Pleu… Pleuro…      1
    ##  2 ASV26    100      168 Eukaryo… Chord… Acti… ""    Embio… Cyma… Cymato…      1
    ##  3 ASV42    100      169 Eukaryo… Chord… Acti… "Per… Hexag… Ophi… Ophiod…      1
    ##  4 ASV48     98.2    170 Eukaryo… Chord… Acti… "Sal… Salmo… Onco… Oncorh…      1
    ##  5 ASV92    100      173 Eukaryo… Chord… Acti… "Aul… Notos… Scop… Scopel…      1
    ##  6 ASV116    98.3    180 Eukaryo… Chord… Acti… "Aci… Acipe… Acip… Acipen…      1
    ##  7 ASV122    99.4    169 Eukaryo… Chord… Acti… "Per… Hexag… Pleu… Pleuro…      1
    ##  8 ASV125   100      173 Eukaryo… Chord… Acti… "Clu… Engra… Engr… Engrau…      1
    ##  9 ASV126   100      168 Eukaryo… Chord… Mamm… "Car… Otari… Eume… Eumeto…      1
    ## 10 ASV128    99.4    169 Eukaryo… Chord… Acti… "Per… Hexag… Pleu… Pleuro…      1
    ## # ℹ 115 more rows

``` r
# remove the single-species seqs from the dataframe and then rank the hits by % identity for the remaining seqs
seq_id_diff <- fixed_tax_df %>%
  anti_join(., single_spp_seqs) %>%
  select(-length) %>%
  group_by(qseqid, species, genus, family, order, class, phylum, kingdom) %>%
    mutate(seq_percID = max(perc_id)) %>%
    group_by(qseqid, species, genus, family, order, class, phylum, kingdom, seq_percID) %>%
  summarise(max(seq_percID)) %>% # take just the top hit for each taxon (for each sequence)
  select(-`max(seq_percID)`) %>%
  ungroup() %>%
  group_by(qseqid) %>%
      mutate(id_rank = rank(desc(seq_percID), ties.method = "min")) %>% # rank the taxonomic hits per sequence by % id
       mutate(top_perc = max(seq_percID)) %>% # designate the highest % id for the best taxonomic hit in each sequence (in some, but not all cases, this is 100%)   
      mutate(diff = top_perc - seq_percID) %>% # calculate the difference between the % identity of the top hit and each subsequent taxonomic hit
        arrange(diff)
```

    ## Joining with `by = join_by(qseqid, perc_id, length, kingdom, phylum, class,
    ## order, family, genus, species)`
    ## `summarise()` has grouped output by 'qseqid', 'species', 'genus', 'family',
    ## 'order', 'class', 'phylum', 'kingdom'. You can override using the `.groups`
    ## argument.

``` r
seq_id_diff %>%
  filter(diff > 0)
```

    ## # A tibble: 3,496 × 12
    ## # Groups:   qseqid [442]
    ##    qseqid species     genus family order class phylum kingdom seq_percID id_rank
    ##    <chr>  <chr>       <chr> <chr>  <chr> <chr> <chr>  <chr>        <dbl>   <int>
    ##  1 ASV629 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… Eukary…       99.4       2
    ##  2 ASV29  Oncorhynch… Onco… Salmo… Salm… Acti… Chord… Eukary…       99.4       2
    ##  3 ASV72  Oncorhynch… Onco… Salmo… Salm… Acti… Chord… Eukary…       99.4       2
    ##  4 ASV245 Microgadus… Micr… Gadid… Gadi… Acti… Chord… Eukary…       98.8       2
    ##  5 ASV408 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… Eukary…       98.8       2
    ##  6 ASV43  Liopsetta … Liop… Pleur… Pleu… Acti… Chord… Eukary…       99.4       3
    ##  7 ASV43  Myzopsetta… Myzo… Pleur… Pleu… Acti… Chord… Eukary…       99.4       3
    ##  8 ASV43  Myzopsetta… Myzo… Pleur… Pleu… Acti… Chord… Eukary…       99.4       3
    ##  9 ASV43  Platichthy… Plat… Pleur… Pleu… Acti… Chord… Eukary…       99.4       3
    ## 10 ASV43  Platichthy… Plat… Pleur… Pleu… Acti… Chord… Eukary…       99.4       3
    ## # ℹ 3,486 more rows
    ## # ℹ 2 more variables: top_perc <dbl>, diff <dbl>

Now I have the single best entry for each species for each sequence
ranked and with the difference between the first and second ranked
entries calculated.

For sequences with multiple top hits, where the difference between
ranked taxa = 0, I will end up defaulting to genus- or family-level ID
(or carrying the individual species info around in some capacity). I
will do the same for any sequences where the difference betweeen the
first and second ranked taxa is \< 2%.

Figure out which differences are \> 2% and eliminate those first?

``` r
# filter out any taxa that are >2% less matching identity than the top taxonomic hit for a given sequence
to_remove_low_perc_hits <- seq_id_diff %>%
  ungroup() %>%
  group_by(qseqid) %>%
  filter(diff > 2)

keepers <- seq_id_diff %>%
  anti_join(to_remove_low_perc_hits)
```

    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`

``` r
# this data frame includes only those taxonomic hits that should be considered.
# so now I need to determine whether they should be assigned to genus, family, order, etc. 
singletons <- keepers %>%
  select(qseqid) %>%
  tally() %>%
  filter(n == 1)

# these are the seqs that now have only a single match
singleton_df <- singletons %>%
  left_join(keepers) %>%
  select(-n) %>%
  bind_rows(single_spp_seqs) %>% # combine the single spp data
  mutate(taxonomic_level = "species") %>%
  mutate(taxon = species)
```

    ## Joining with `by = join_by(qseqid)`

``` r
## Genus-level matches
# remove the singletons from the bigger df 
single_genus <- keepers %>%
  anti_join(singleton_df)%>% # at best, these should be genus-level matches
  group_by(qseqid, genus) %>%
  tally() %>%
  ungroup() %>%
  group_by(qseqid) %>%
  tally() %>%
  filter(n == 1) %>% # seqs that match a single genus
  select(-n) %>%
  left_join(., keepers) %>%
  mutate(taxonomic_level = "genus") %>%
  mutate(taxon = genus)
```

    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid)`

``` r
## Family-level matches
single_family <- keepers %>%
  anti_join(singleton_df)%>%
  anti_join(single_genus) %>%
  group_by(qseqid, family) %>%
  tally() %>%
  ungroup() %>%
  group_by(qseqid) %>%
  tally() %>%
  filter(n == 1) %>% # seqs that match a single family
  select(-n) %>%
  left_join(., keepers) %>%
  mutate(taxonomic_level = "family") %>%
  mutate(taxon = family)
```

    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid)`

``` r
## Order-level matches
single_order <- keepers %>%
  anti_join(singleton_df)%>%
  anti_join(single_genus) %>%
  anti_join(single_family) %>%
  group_by(qseqid, order) %>%
  tally() %>%
  ungroup() %>%
  group_by(qseqid) %>%
  tally() %>%
  filter(n == 1) %>% # seqs that match a single order
  select(-n) %>%
  left_join(., keepers) %>%
  mutate(taxonomic_level = "order") %>%
  mutate(taxon = order)
```

    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid)`

``` r
## Class-level matches
single_class <- keepers %>%
  anti_join(singleton_df)%>%
  anti_join(single_genus) %>%
  anti_join(single_family) %>%
  anti_join(single_order) %>%
  group_by(qseqid, class) %>%
  tally() %>%
  ungroup() %>%
  group_by(qseqid) %>%
  tally() %>% 
  filter(n == 1) %>% # seqs that match a single class
  select(-n) %>%
  left_join(., keepers) %>%
  mutate(taxonomic_level = "class") %>%
  mutate(taxon = class)
```

    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid, species, genus, family, order, class,
    ## phylum, kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(qseqid)`

``` r
# no higher level taxonomy was relevant in this dataset.
```

only 6 ASVs stuck at class

Modify the singleton_df to include the right variable headers

``` r
single_spp <- singleton_df %>%
  select(-perc_id, -length, -n_taxa) %>%
  mutate(taxonomic_level = "species") %>%
  mutate(taxon = species)
```

``` r
# recombine the full data set now that the appropriate level of taxonomy has been determined
sorted_tax_df <- bind_rows(single_class, single_order, single_family, single_genus, single_spp)
```

Kim has code to verify that species are within reasonable
range/distributions. – I’m adding that below:

# now let’s take a closer look at the assignments we are getting from inital blastn

## what non-fish are here?

``` r
not_Actinopteri <- sorted_tax_df %>%
  filter(class != "Actinopteri") %>%
  select(species, genus, family, order, class, phylum) %>%
  unique()

not_Actinopteri
```

    ## # A tibble: 3 × 6
    ##   species             genus       family      order          class        phylum
    ##   <chr>               <chr>       <chr>       <chr>          <chr>        <chr> 
    ## 1 Callorhinus ursinus Callorhinus Otariidae   Carnivora      Mammalia     Chord…
    ## 2 Hydrolagus colliei  Hydrolagus  Chimaeridae Chimaeriformes Chondrichth… Chord…
    ## 3 Eumetopias jubatus  Eumetopias  Otariidae   Carnivora      Mammalia     Chord…

all of those make sense.

## remove terrestrial/freshwater/out-of-range non-fish from data set

``` r
# not_Actinopteri_keepers <- not_Actinopteri %>% 
#   #class = Asteroidea - species = Asterias amurensis (North Pacific seastar) - seem reasonable to keep
#   #class = Aves - keeping murre and murrelet
#   filter(order != "Galliformes") %>%  ### remove chicken, junglefowl 
#   filter(species != "Harpagornis moorei") %>% ## definitely no extinct eagle from New Zealand... 
#   #family = Petromyzontidae - lamprey
#   filter(!(family == "Petromyzontidae" & species != "Lethenteron camtschaticum")) %>% ### remove lamprey except Arctic lamprey 
#   #picoplankton
#   filter(species != "Bathycoccus prasinos") %>% ### remove a picoplankton
#   #class = mammals
#   #order = Artiodactyla
#   filter(!family == "Bovidae") %>% ### remove bovids
#   filter(!(family == "Delphinidae" & species != "Lagenorhynchus obliquidens")) %>% ### remove dolphins expect Pacific white-sided 
#   filter(!family == "Suidae") %>% ### remove pigs 
#   filter(!family == "Cervidae") %>% ### remove moose 
#   #order = Carnivora
#   filter(!species == "Pusa sibirica") %>% ### remove Baikal seal 
#   filter(!species == "Pusa caspica") %>% ### remove Caspian seal 
#   filter(!family == "Hominidae") %>% ### remove humans
#   filter(!family == "Canidae") %>% ### remove dog
#   select(species) %>%
#   rename(Species = species) %>%
#   mutate(in_range = "yes")
```

## now look at the fish and figure out what taxa are in/out of our range

``` r
to_check_range <- sorted_tax_df %>%
  filter(class == "Actinopteri") %>%
  select(species, genus, family, order, class, phylum) %>%
  unique()
```

## check ranges for species using rfishbase

## also at this step, check

``` r
# library(remotes)
# remotes::install_github("ropensci/rfishbase")
# library(rfishbase)
# 
# #one-time download of all fishbase tables... this takes a bit 
# rfishbase::fb_import()
# 
# fb_tables("fishbase")
# 
# #first, validate species names using rfishbase synonyms
# spp_df <- synonyms(to_check_range$species)
# 
# syn <- spp_df %>% 
#   filter(Status == "synonym")
#   
# to_check_range <- to_check_range %>% 
#   mutate(validated_name = ifelse(species %in% syn$synonym, syn$Species, species))
#     
# to_check_range %>%
#   filter(species != validated_name)

# THE DISTRIBUTION FUNCTION ISN'T WORKING FOR ME?
#get distribution info 
# spp_distribution <- distribution(to_check_range$validated_name) %>%
#   select(Species, FAO) %>%
#   unique()

#add column to designate if we will consider a species as "in range"- for this study, this will be NE Pacific and Arctic Ocean 
# spp_distribution_range <- spp_distribution %>%
#   mutate(in_range = ifelse(is.na(FAO), NA, "no"),
#          in_range = ifelse(FAO == "Pacific, Northeast", "yes", in_range),
#          in_range = ifelse(FAO == "Arctic Ocean", "yes", in_range))
# 
# #keep just a list of spp names and yes/no/NA for "in range"  - this way we can keep track of what spp didn't have any reference information in fishbase to determine range 
# spp_range <- spp_distribution_range %>%
#   select(Species, in_range) %>%
#   unique()

#how many entries do not have range info? 
# range_na <- spp_range %>%
#   filter(is.na(in_range))
```

When a valid name was not found, the presence of a species in the study
area was checked using the GBIF database (<https://www.gbif.org/>).

# some species do not have range info - manually determine if these species should be considered in range

``` r
# range_na <- range_na %>%
#   mutate(in_range = ifelse(Species == "Ammodytes japonicus", "no", in_range),
#          #in_range = ifelse(Species == "Cleisthenes herzensteini", "no", in_range),
#          #in_range = ifelse(Species == "Gadus ogac", "no", in_range),
#          #in_range = ifelse(Species == "Myoxocephalus aenaeas", "no", in_range),
#          #in_range = ifelse(Species == "Kareius bicoloratus", "no", in_range),
#          #in_range = ifelse(Species == "Cottocomephorus grewingki", "no", in_range),
#          in_range = ifelse(Species == "Sebastes cheni", "no", in_range)) #,
#          #in_range = ifelse(Species == "Embassichthys bathybius", "yes", in_range),    ##new name is Microstomus bathybius
#          #in_range = ifelse(Species == "Pungitius kaibarae", "no", in_range),
#          #in_range = ifelse(Species == "Ulcina olrikii", "yes", in_range),  ### arctic alligatorfish - Aspidophoroides olrikii
#          #in_range = ifelse(Species == "Polypera greeni", "yes", in_range))  ## Liparis greeni
```

``` r
to_check_range
```

    ## # A tibble: 271 × 6
    ##    species                      genus           family        order class phylum
    ##    <chr>                        <chr>           <chr>         <chr> <chr> <chr> 
    ##  1 Liopsetta glacialis          Liopsetta       Pleuronectid… Pleu… Acti… Chord…
    ##  2 Myzopsetta ferruginea        Myzopsetta      Pleuronectid… Pleu… Acti… Chord…
    ##  3 Platichthys flesus           Platichthys     Pleuronectid… Pleu… Acti… Chord…
    ##  4 Platichthys stellatus        Platichthys     Pleuronectid… Pleu… Acti… Chord…
    ##  5 Cleisthenes herzensteini     Cleisthenes     Pleuronectid… Pleu… Acti… Chord…
    ##  6 Cleisthenes pinetorum        Cleisthenes     Pleuronectid… Pleu… Acti… Chord…
    ##  7 Hippoglossoides dubius       Hippoglossoides Pleuronectid… Pleu… Acti… Chord…
    ##  8 Hippoglossoides elassodon    Hippoglossoides Pleuronectid… Pleu… Acti… Chord…
    ##  9 Hippoglossoides platessoides Hippoglossoides Pleuronectid… Pleu… Acti… Chord…
    ## 10 Hippoglossoides robustus     Hippoglossoides Pleuronectid… Pleu… Acti… Chord…
    ## # ℹ 261 more rows

Create output taxonomy data frames

``` r
uncollapsed_taxonomy <- sorted_tax_df %>%
  select(-top_perc, -id_rank) %>%
  unique()


# quick look at the taxonomy before collapsing it:
uncollapsed_taxonomy
```

    ## # A tibble: 5,093 × 12
    ##    qseqid species       genus family order class phylum kingdom seq_percID  diff
    ##    <chr>  <chr>         <chr> <chr>  <chr> <chr> <chr>  <chr>        <dbl> <dbl>
    ##  1 ASV386 Liopsetta gl… Liop… Pleur… Pleu… Acti… Chord… Eukary…      100   0    
    ##  2 ASV386 Myzopsetta f… Myzo… Pleur… Pleu… Acti… Chord… Eukary…      100   0    
    ##  3 ASV386 Platichthys … Plat… Pleur… Pleu… Acti… Chord… Eukary…      100   0    
    ##  4 ASV386 Platichthys … Plat… Pleur… Pleu… Acti… Chord… Eukary…      100   0    
    ##  5 ASV386 Cleisthenes … Clei… Pleur… Pleu… Acti… Chord… Eukary…       99.4 0.592
    ##  6 ASV386 Cleisthenes … Clei… Pleur… Pleu… Acti… Chord… Eukary…       99.4 0.592
    ##  7 ASV386 Hippoglossoi… Hipp… Pleur… Pleu… Acti… Chord… Eukary…       99.4 0.592
    ##  8 ASV386 Hippoglossoi… Hipp… Pleur… Pleu… Acti… Chord… Eukary…       99.4 0.592
    ##  9 ASV386 Hippoglossoi… Hipp… Pleur… Pleu… Acti… Chord… Eukary…       99.4 0.592
    ## 10 ASV386 Hippoglossoi… Hipp… Pleur… Pleu… Acti… Chord… Eukary…       99.4 0.592
    ## # ℹ 5,083 more rows
    ## # ℹ 2 more variables: taxonomic_level <chr>, taxon <chr>

``` r
# and then collapse that down to just a single taxon per ASV
collapsed_taxonomy <- uncollapsed_taxonomy %>%
  select(qseqid, taxon, taxonomic_level) %>%
  unique()

# then remove those high level taxonomic assignments that aren't informative
collapsed_taxonomy %>%
  #filter(!taxonomic_level %in% c("class", "order")) %>%
  select(-qseqid) %>%
  unique()
```

    ## # A tibble: 70 × 2
    ##    taxon          taxonomic_level
    ##    <chr>          <chr>          
    ##  1 Actinopteri    class          
    ##  2 Perciformes    order          
    ##  3 Hexagrammidae  family         
    ##  4 Gadidae        family         
    ##  5 Sebastidae     family         
    ##  6 Pleuronectidae family         
    ##  7 Ammodytidae    family         
    ##  8 Stichaeidae    family         
    ##  9 Clupeidae      family         
    ## 10 Cottidae       family         
    ## # ℹ 60 more rows

``` r
collapsed_taxonomy %>%
  select(-qseqid) %>%
  unique() %>%
  group_by(taxonomic_level) %>%
  tally()
```

    ## # A tibble: 5 × 2
    ##   taxonomic_level     n
    ##   <chr>           <int>
    ## 1 class               1
    ## 2 family             10
    ## 3 genus              23
    ## 4 order               1
    ## 5 species            35

Conservatively, we have 35 species, 23 genera, 10 family-level taxonomic
assignments.

#### Read in ASV table

``` r
asv_table <- read.csv("csv_outputs/NFSfecal_MiFish_ASVtable.csv") %>%
  rename(sample = X)

# reformat
sample_df <- asv_table %>%
  pivot_longer(2:length(asv_table), names_to = "ASV", values_to = "reads") %>%
  filter(reads > 0) %>%
  separate(sample, into = c("ablg", "plate", "replicate"), remove = F)
```

    ## Warning: Expected 3 pieces. Missing pieces filled with `NA` in 66 rows [4159, 4160,
    ## 4161, 4162, 4163, 4164, 4165, 4166, 4167, 4168, 4169, 4170, 4171, 4172, 4173,
    ## 4174, 4175, 4176, 4177, 4178, ...].

``` r
head(sample_df)
```

    ## # A tibble: 6 × 6
    ##   sample          ablg   plate  replicate ASV   reads
    ##   <chr>           <chr>  <chr>  <chr>     <chr> <int>
    ## 1 e03719-plate1-a e03719 plate1 a         ASV1  31055
    ## 2 e03719-plate1-a e03719 plate1 a         ASV3    654
    ## 3 e03719-plate1-a e03719 plate1 a         ASV6    318
    ## 4 e03719-plate1-a e03719 plate1 a         ASV23   216
    ## 5 e03719-plate1-a e03719 plate1 a         ASV27  5452
    ## 6 e03719-plate1-a e03719 plate1 a         ASV40   144

``` r
sample_df %>%
  filter(plate == "plate1") %>%
  filter(reads > 10) %>%
  ggplot(aes(x = ablg, y = reads, fill = ASV)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90)
  ) +
  labs( 
    title = "NFS fecal samples - plate 1")
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/plate-1-1.png)<!-- -->

``` r
sample_df %>%
  filter(plate == "plate2") %>%
  filter(reads > 10) %>%
  ggplot(aes(x = ablg, y = reads, fill = ASV)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90)
  ) +
  labs( 
    title = "NFS fecal samples - plate 2")
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/plate-2-1.png)<!-- -->

Just quickly, what are the top species?

``` r
combined_df <- sample_df %>%
  left_join(., collapsed_taxonomy, by = c("ASV" = "qseqid"))


combined_df %>%
  arrange(desc(reads)) %>%
  filter(!taxon %in% c("Callorhinus ursinus", "Acipenser", "Acipenser fulvescens", "Actinopteri", "uncultured bacterium", "Eumetopias jubatus")) %>%
  group_by(taxon) %>%
  summarise(sum(reads)) %>%
  arrange(desc(`sum(reads)`))
```

    ## # A tibble: 66 × 2
    ##    taxon                       `sum(reads)`
    ##    <chr>                              <int>
    ##  1 Leuroglossus                     1038077
    ##  2 Oncorhynchus                      767699
    ##  3 Pleurogrammus monopterygius       256001
    ##  4 Gadidae                           201526
    ##  5 Sebastidae                         52174
    ##  6 <NA>                               25722
    ##  7 Stenobrachius                      19463
    ##  8 Cymatogaster aggregata             17452
    ##  9 Pleuronectidae                     13628
    ## 10 Stichaeidae                         8244
    ## # ℹ 56 more rows

``` r
combined_df %>%
  filter(!taxon %in% c("Callorhinus ursinus", "Acipenser", "Acipenser fulvescens", "Actinopteri", "uncultured bacterium", "Eumetopias jubatus")) %>%
  filter(reads > 10) %>%
  group_by(taxon) %>%
  tally() %>%
  arrange(desc(n))
```

    ## # A tibble: 61 × 2
    ##    taxon                           n
    ##    <chr>                       <int>
    ##  1 Leuroglossus                  721
    ##  2 Oncorhynchus                  481
    ##  3 <NA>                          295
    ##  4 Sebastidae                    278
    ##  5 Gadidae                       238
    ##  6 Pleurogrammus monopterygius   206
    ##  7 Pleuronectidae                120
    ##  8 Stenobrachius                 112
    ##  9 Cymatogaster aggregata        105
    ## 10 Stichaeidae                    78
    ## # ℹ 51 more rows

How many sturgeon reads were in the samples? (tag-jumping)

``` r
combined_df %>%
  filter(taxon == "Acipenser" | taxon == "Acipenser fulvescens") %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  arrange(desc(reads)) %>%
  group_by(plate) %>%
  summarise(max(reads))
```

    ## # A tibble: 3 × 2
    ##   plate  `max(reads)`
    ##   <chr>         <int>
    ## 1 plate1          184
    ## 2 plate2          501
    ## 3 plate3          181

Ok, so as few as 3 reads and as many as 501 in the samples. Not super
clean… but we’ll proceed for now. Notably, plates 1 and 3 are very
similar, and plate 2 has more apparent cross-contamination/tag-jumping.

``` r
# what about mean values
combined_df %>%
  filter(taxon == "Acipenser" | taxon == "Acipenser fulvescens") %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  arrange(desc(reads)) %>%
  group_by(plate) %>%
  summarise(mean(reads))
```

    ## # A tibble: 3 × 2
    ##   plate  `mean(reads)`
    ##   <chr>          <dbl>
    ## 1 plate1          94.5
    ## 2 plate2          70.5
    ## 3 plate3          73.7

Interesting! So the mean is actually higher in plate 1.

What’s the distribution of read depth across plates though?

``` r
combined_df %>%
  group_by(plate) %>%
  summarise(mean(reads))
```

    ## # A tibble: 3 × 2
    ##   plate  `mean(reads)`
    ##   <chr>          <dbl>
    ## 1 plate1         4278.
    ## 2 plate2         1328.
    ## 3 plate3         3013.

``` r
combined_df %>%
  group_by(plate) %>%
  summarise(mean(reads))
```

    ## # A tibble: 3 × 2
    ##   plate  `mean(reads)`
    ##   <chr>          <dbl>
    ## 1 plate1         4278.
    ## 2 plate2         1328.
    ## 3 plate3         3013.

``` r
#range?
combined_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  group_by(plate) %>%
  summarise(max(reads))
```

    ## # A tibble: 3 × 2
    ##   plate  `max(reads)`
    ##   <chr>         <int>
    ## 1 plate1        58759
    ## 2 plate2        59408
    ## 3 plate3        38293

plate 1 has a much higher mean.

``` r
clean_df <- combined_df %>%
  filter(!is.na(taxon)) %>% # remove entries for taxa that we filtered out (human, dog, bacteria, etc.)
  select(-ASV) %>%
  group_by(sample, taxon) %>%
  mutate(taxon_reads = sum(reads)) %>%
  select(-reads) %>% 
  unique() %>% # remove duplicate entries for taxa with multiple ASVs
  ungroup() %>%
  group_by(sample) %>%
  mutate(total_reads = sum(taxon_reads))
  

clean_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(total_reads > 5000) %>% # no samples had fewer than 5000 reads. Great!
  ggplot(aes(x = sample, y = taxon_reads, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_wrap(~plate) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90)
  ) 
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/need-to-combine-multiple-ASVs-for-the-same-species-1.png)<!-- -->

``` r
ggsave("pdf_outputs/fecal_samples_by_plate.pdf", width = 8, height = 5)
```

``` r
# how many samples fell below the 1000 read threshold?
clean_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(plate == "plate1") %>%
  ggplot(aes(x = ablg, y = taxon_reads, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(row = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90)
  ) +
  labs(y = "reads",
       title = "plate1")
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/consistency-across-extractions-from-the-same-sample?-1.png)<!-- -->
To actually evaluate this question, we’ll want to do a Bray-Curtis or
similar… probably would be useful to consult with Kim.

``` r
# how many samples fell below the 1000 read threshold?
clean_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(plate == "plate2") %>%
  ggplot(aes(x = ablg, y = taxon_reads, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(row = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90)
  ) +
  labs(y = "reads",
       title = "plate2")
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
# how many samples fell below the 1000 read threshold?
clean_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(plate == "plate3") %>%
  ggplot(aes(x = ablg, y = taxon_reads, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(row = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90)
  ) +
  labs(y = "reads",
       title = "plate3")
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Could also finagle these as proportions:

``` r
# first estimate the % of tag-jumping on a per-sample basis
prop_df <- clean_df %>%
  mutate(taxon_prop = taxon_reads/total_reads)
  
 
# and then establish the prop of sturgeon
prop_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(taxon == "Acipenser" | taxon == "Acipenser fulvescens") %>%
  group_by(plate) %>%
  summarise(max(taxon_prop))
```

    ## # A tibble: 3 × 2
    ##   plate  `max(taxon_prop)`
    ##   <chr>              <dbl>
    ## 1 plate1           0.00344
    ## 2 plate2           0.0168 
    ## 3 plate3           0.00889

mean proportion of tag-jumping/sturgeon is 0.2% across all three plate.
max proportion = 1.7% in plate 2, but otherwise, under 1% How many
samples are \>1%?

``` r
prop_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(taxon == "Acipenser" | taxon == "Acipenser fulvescens") %>%
  filter(taxon_prop > 0.01)
```

    ## # A tibble: 2 × 9
    ## # Groups:   sample [2]
    ##   sample     ablg  plate replicate taxon taxonomic_level taxon_reads total_reads
    ##   <chr>      <chr> <chr> <chr>     <chr> <chr>                 <int>       <int>
    ## 1 e03752-pl… e037… plat… b         Acip… genus                   494       34712
    ## 2 e03754-pl… e037… plat… b         Acip… genus                   655       38947
    ## # ℹ 1 more variable: taxon_prop <dbl>

Just two samples: 1.4% and 1.6% of reads.

So for filtering: If we used 0.2% cut-off vs. a 1% cut-off, what do we
lose in terms of species/data?

``` r
# as an experiment, it might be worth analyzing this both ways...
prop_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(taxon != "Acipenser" | taxon != "Acipenser fulvescens") %>%
  filter(taxon_prop < 0.002)
```

    ## # A tibble: 1,112 × 9
    ## # Groups:   sample [231]
    ##    sample    ablg  plate replicate taxon taxonomic_level taxon_reads total_reads
    ##    <chr>     <chr> <chr> <chr>     <chr> <chr>                 <int>       <int>
    ##  1 e03719-p… e037… plat… b         Onco… genus                    37       34216
    ##  2 e03719-p… e037… plat… b         Sten… genus                    38       34216
    ##  3 e03719-p… e037… plat… a         Acip… genus                    13       30615
    ##  4 e03719-p… e037… plat… a         Cyma… species                  20       30615
    ##  5 e03719-p… e037… plat… a         Cott… family                   18       30615
    ##  6 e03719-p… e037… plat… a         Stic… family                   56       30615
    ##  7 e03719-p… e037… plat… a         Gast… genus                    19       30615
    ##  8 e03719-p… e037… plat… a         Blep… genus                    19       30615
    ##  9 e03719-p… e037… plat… a         Enop… species                  24       30615
    ## 10 e03719-p… e037… plat… a         Hexa… family                    9       30615
    ## # ℹ 1,102 more rows
    ## # ℹ 1 more variable: taxon_prop <dbl>

``` r
prop_df %>%
  filter(ablg %in% c("PC") & !replicate %in% c("EB", "FB")) %>%
  filter(taxon == "Acipenser" | taxon == "Acipenser fulvescens")
```

    ## # A tibble: 4 × 9
    ## # Groups:   sample [3]
    ##   sample    ablg  plate  replicate taxon taxonomic_level taxon_reads total_reads
    ##   <chr>     <chr> <chr>  <chr>     <chr> <chr>                 <int>       <int>
    ## 1 PC-plate1 PC    plate1 <NA>      Acip… genus                 24991       25664
    ## 2 PC-plate2 PC    plate2 <NA>      Acip… genus                 28940       29980
    ## 3 PC-plate2 PC    plate2 <NA>      Acip… species                 329       29980
    ## 4 PC-plate3 PC    plate3 <NA>      Acip… genus                 18906       19337
    ## # ℹ 1 more variable: taxon_prop <dbl>

Filtering out species that occur at \< 0.2% of the reads:

``` r
prop_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(!taxon %in% c("Acipenser", "Acipenser fulvescens", "uncultured bacterium")) %>%
  filter(taxon_prop > 0.002) %>%
  filter(plate == "plate1") %>%
  ggplot(aes(x = ablg, y = taxon_prop, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(row = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5))
  ) +
  labs(y = "read proportion",
       title = "plate1, filtering threshold 0.2%",
       x = "sample")
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/plate1-w-filter-1.png)<!-- -->

``` r
ggsave("pdf_outputs/plate1_prop_02filter.pdf", width = 8, height = 5)
```

``` r
prop_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(!taxon %in% c("Acipenser", "Acipenser fulvescens", "uncultured bacterium")) %>%
  filter(taxon_reads > 10) %>%
  filter(plate == "plate1") %>%
  ggplot(aes(x = ablg, y = taxon_prop, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(row = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5))
  ) +
  labs(y = "read proportion",
       title = "plate1, 10 reads filtering threshold",
       x = "sample")
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/plate1-10-read-filter-1.png)<!-- -->

``` r
ggsave("pdf_outputs/plate1_prop_10reads_filter.pdf", width = 8, height = 5)
```

``` r
prop_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(!taxon %in% c("Acipenser", "Acipenser fulvescens", "uncultured bacterium")) %>%
  #filter(taxon_reads > 10) %>%
  filter(plate == "plate1") %>%
  ggplot(aes(x = ablg, y = taxon_prop, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(row = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5))
  ) +
  labs(y = "read proportion",
       title = "plate1, no filtering threshold",
       x = "sample")
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/plate1-wo-read-filter-1.png)<!-- -->

``` r
ggsave("pdf_outputs/plate1_prop_no_filter.pdf", width = 8, height = 5)
```

``` r
prop_df %>%
  filter(!ablg %in% c("PC", "NC") & !replicate %in% c("EB", "FB")) %>%
  filter(!taxon %in% c("Acipenser", "Acipenser fulvescens", "uncultured bacterium")) %>%
  filter(taxon_prop > 0.01) %>%
  filter(plate == "plate1") %>%
  ggplot(aes(x = ablg, y = taxon_prop, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(row = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5))
  ) +
  labs(y = "read proportion",
       title = "plate1, 1% filtering threshold",
       x = "sample")
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/plate1-1-percent-filter-1.png)<!-- -->

``` r
combined_df %>%
  filter(reads > 20) %>%
  #filter(!taxon %in% c("Callorhinus ursinus", "Acipenser", "Acipenser fulvescens", "Actinopteri", "uncultured bacterium", "Eumetopias jubatus")) %>%
  ggplot(aes(x = ablg, y = reads, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(plate), cols = vars(replicate), space = "free", scales = "free") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90)
  ) 
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/do-technical-plate-reps-look-identical-1.png)<!-- -->

``` r
combined_df %>%
  filter(replicate %in% c(
    "EB"
  )) %>%
  filter(reads > 20) %>%
   ggplot(aes(x = ablg, y = reads, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(plate), cols = vars(replicate), space = "free", scales = "free") +
  theme(
   # legend.position = "none",
    axis.text.x = element_text(angle = 90)
  )
```

![](fecal-mifish-taxon-analysis_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->
