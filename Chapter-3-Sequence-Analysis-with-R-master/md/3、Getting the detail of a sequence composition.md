
# Getting the detail of a sequence composition
After we have retrieved a sequence, we need to know more about it, such as the nucleotide or amino acid frequency, and Guanine and Cytosine nucleotide bases (GC content). The contents of a sequence also help in determining certain properties of the entire molecule, for example, the acid basicity hydrophobicity in proteins based on the amino acids present in the sequence. An interesting aspect is the GC content in a nucleotide. It refers to the fraction of Guanine (G) and Cytosine (C) in the sequence, and certain genomes, especially among the bacteria, show a significant difference on this scale and variations in terms of genomic regions. To illustrate, some Actinobacteria can have more than 70 percent of GC, whereas some Proteobacteria can have less that 20 percent of GC. Furthermore, the GC content is also used to predict the annealing temperature of the sequence during PCR experiments. These aspects make the
analysis of the contents in a sequence important. 


```R
library('seqinr')
```


```R
choosebank("genbank")
```


```R
# Gets the M.tuberculi sequences
q1 <- query(listname = "actino", query="SP=Mycobacterium tuberculosis AND K=rpoB") 
q1
```


    943 SQ for SP=Mycobacterium tuberculosis AND K=rpoB



```R
# Gets the E.coli sequences
q2 <-query(listname = "proteo", query="SP=Escherichia coli AND K=rpoB")
q2
```


    492 SQ for SP=Escherichia coli AND K=rpoB



```R
q1$req   # 866 Sequences as on Jan 29th 2018
```


    [[1]]
          name     length      frame     ncbicg 
    "AB711167"      "128"        "2"       "11" 
    
    [[2]]
          name     length      frame     ncbicg 
    "AB711168"      "128"        "2"       "11" 
    
    [[3]]
          name     length      frame     ncbicg 
    "AB711169"      "128"        "2"       "11" 
    
    [[4]]
          name     length      frame     ncbicg 
    "AB711170"      "128"        "2"       "11" 
    
    [[5]]
          name     length      frame     ncbicg 
    "AB711171"      "128"        "2"       "11" 
    
    [[6]]
          name     length      frame     ncbicg 
    "AB711172"      "128"        "2"       "11" 
    
    [[7]]
          name     length      frame     ncbicg 
    "AB711173"      "128"        "2"       "11" 
    
    [[8]]
          name     length      frame     ncbicg 
    "AB711174"      "128"        "2"       "11" 
    
    [[9]]
          name     length      frame     ncbicg 
    "AB711175"      "128"        "2"       "11" 
    
    [[10]]
          name     length      frame     ncbicg 
    "AB711176"      "128"        "2"       "11" 
    
    [[11]]
          name     length      frame     ncbicg 
    "AB711177"      "128"        "2"       "11" 
    
    [[12]]
          name     length      frame     ncbicg 
    "AB711178"      "128"        "2"       "11" 
    
    [[13]]
               name          length           frame          ncbicg 
    "AE000516.RPOB"          "3537"             "0"            "11" 
    
    [[14]]
          name     length      frame     ncbicg 
    "AF055891"      "278"        "0"       "11" 
    
    [[15]]
          name     length      frame     ncbicg 
    "AF055892"      "278"        "0"       "11" 
    
    [[16]]
          name     length      frame     ncbicg 
    "AF055893"      "278"        "0"       "11" 
    
    [[17]]
          name     length      frame     ncbicg 
    "AF057454"      "306"        "2"       "11" 
    
    [[18]]
          name     length      frame     ncbicg 
    "AF060353"      "705"        "2"       "11" 
    
    [[19]]
          name     length      frame     ncbicg 
    "AF112973"      "112"        "2"       "11" 
    
    [[20]]
          name     length      frame     ncbicg 
    "AF143771"      "156"        "0"       "11" 
    
    [[21]]
          name     length      frame     ncbicg 
    "AF146567"      "147"        "0"       "11" 
    
    [[22]]
          name     length      frame     ncbicg 
    "AF147030"      "144"        "0"       "11" 
    
    [[23]]
          name     length      frame     ncbicg 
    "AF147031"      "156"        "0"       "11" 
    
    [[24]]
          name     length      frame     ncbicg 
    "AF147033"      "156"        "0"       "11" 
    
    [[25]]
          name     length      frame     ncbicg 
    "AF147034"      "153"        "0"       "11" 
    
    [[26]]
          name     length      frame     ncbicg 
    "AF177294"      "334"        "0"       "11" 
    
    [[27]]
          name     length      frame     ncbicg 
    "AF312232"      "157"        "2"       "11" 
    
    [[28]]
          name     length      frame     ncbicg 
    "AF312233"      "157"        "2"       "11" 
    
    [[29]]
          name     length      frame     ncbicg 
    "AF312234"      "157"        "2"       "11" 
    
    [[30]]
          name     length      frame     ncbicg 
    "AF312235"      "157"        "2"       "11" 
    
    [[31]]
          name     length      frame     ncbicg 
    "AF312236"      "157"        "2"       "11" 
    
    [[32]]
          name     length      frame     ncbicg 
    "AF360399"      "278"        "0"       "11" 
    
    [[33]]
          name     length      frame     ncbicg 
    "AF360400"      "278"        "0"       "11" 
    
    [[34]]
          name     length      frame     ncbicg 
    "AF360401"      "278"        "0"       "11" 
    
    [[35]]
          name     length      frame     ncbicg 
    "AF503174"       "81"        "0"       "11" 
    
    [[36]]
          name     length      frame     ncbicg 
    "AF503175"       "81"        "0"       "11" 
    
    [[37]]
          name     length      frame     ncbicg 
    "AF515787"      "162"        "0"       "11" 
    
    [[38]]
          name     length      frame     ncbicg 
    "AF515788"      "120"        "0"       "11" 
    
    [[39]]
          name     length      frame     ncbicg 
    "AF515789"      "120"        "0"       "11" 
    
    [[40]]
          name     length      frame     ncbicg 
    "AF515790"      "120"        "0"       "11" 
    
    [[41]]
          name     length      frame     ncbicg 
    "AF515791"      "120"        "0"       "11" 
    
    [[42]]
          name     length      frame     ncbicg 
    "AF515792"      "120"        "0"       "11" 
    
    [[43]]
          name     length      frame     ncbicg 
    "AF515793"      "120"        "0"       "11" 
    
    [[44]]
          name     length      frame     ncbicg 
    "AF532616"      "252"        "0"       "11" 
    
    [[45]]
          name     length      frame     ncbicg 
    "AF532617"      "267"        "0"       "11" 
    
    [[46]]
          name     length      frame     ncbicg 
    "AJ297922"       "84"        "0"       "11" 
    
    [[47]]
          name     length      frame     ncbicg 
    "AJ297923"       "84"        "0"       "11" 
    
    [[48]]
          name     length      frame     ncbicg 
    "AJ297924"       "81"        "0"       "11" 
    
    [[49]]
          name     length      frame     ncbicg 
    "AJ297925"       "84"        "0"       "11" 
    
    [[50]]
          name     length      frame     ncbicg 
    "AJ297926"       "84"        "0"       "11" 
    
    [[51]]
          name     length      frame     ncbicg 
    "AJ297927"       "84"        "0"       "11" 
    
    [[52]]
          name     length      frame     ncbicg 
    "AJ297928"       "84"        "0"       "11" 
    
    [[53]]
          name     length      frame     ncbicg 
    "AJ297929"       "87"        "0"       "11" 
    
    [[54]]
               name          length           frame          ncbicg 
    "AJ318813.RPOB"           "615"             "0"            "11" 
    
    [[55]]
          name     length      frame     ncbicg 
    "AJ318814"      "633"        "0"       "11" 
    
    [[56]]
          name     length      frame     ncbicg 
    "AJ318815"      "618"        "0"       "11" 
    
    [[57]]
          name     length      frame     ncbicg 
    "AJ318816"      "637"        "0"       "11" 
    
    [[58]]
          name     length      frame     ncbicg 
    "AJ318817"      "618"        "0"       "11" 
    
    [[59]]
          name     length      frame     ncbicg 
    "AJ318818"      "610"        "0"       "11" 
    
    [[60]]
          name     length      frame     ncbicg 
    "AJ318819"      "610"        "0"       "11" 
    
    [[61]]
          name     length      frame     ncbicg 
    "AJ318820"      "105"        "0"       "11" 
    
    [[62]]
          name     length      frame     ncbicg 
    "AJ318821"      "639"        "0"       "11" 
    
    [[63]]
          name     length      frame     ncbicg 
    "AJ418714"       "81"        "0"       "11" 
    
    [[64]]
          name     length      frame     ncbicg 
    "AJ457091"       "81"        "0"       "11" 
    
    [[65]]
          name     length      frame     ncbicg 
    "AJ865088"       "84"        "0"       "11" 
    
    [[66]]
          name     length      frame     ncbicg 
    "AJ868293"      "225"        "0"       "11" 
    
    [[67]]
          name     length      frame     ncbicg 
    "AJ870393"       "84"        "0"       "11" 
    
    [[68]]
          name     length      frame     ncbicg 
    "AJ870394"       "87"        "0"       "11" 
    
    [[69]]
          name     length      frame     ncbicg 
    "AJ870395"       "87"        "0"       "11" 
    
    [[70]]
          name     length      frame     ncbicg 
    "AJ870396"       "84"        "0"       "11" 
    
    [[71]]
          name     length      frame     ncbicg 
    "AJ870397"       "84"        "0"       "11" 
    
    [[72]]
          name     length      frame     ncbicg 
    "AJ874716"      "225"        "0"       "11" 
    
    [[73]]
          name     length      frame     ncbicg 
    "AJ874717"      "225"        "0"       "11" 
    
    [[74]]
          name     length      frame     ncbicg 
    "AJ890105"       "87"        "0"       "11" 
    
    [[75]]
               name          length           frame          ncbicg 
    "AL123456.RPOB"          "3519"             "0"            "11" 
    
    [[76]]
               name          length           frame          ncbicg 
    "AP012340.RPOB"          "3291"             "0"            "11" 
    
    [[77]]
               name          length           frame          ncbicg 
    "AP014573.RPOB"          "3291"             "0"            "11" 
    
    [[78]]
               name          length           frame          ncbicg 
    "AP017901.RPOB"          "3537"             "0"            "11" 
    
    [[79]]
               name          length           frame          ncbicg 
    "AP018033.RPOB"          "3519"             "0"            "11" 
    
    [[80]]
               name          length           frame          ncbicg 
    "AP018034.RPOB"          "3519"             "0"            "11" 
    
    [[81]]
               name          length           frame          ncbicg 
    "AP018035.RPOB"          "3519"             "0"            "11" 
    
    [[82]]
               name          length           frame          ncbicg 
    "AP018036.RPOB"          "3519"             "0"            "11" 
    
    [[83]]
          name     length      frame     ncbicg 
    "AY095448"       "81"        "0"       "11" 
    
    [[84]]
          name     length      frame     ncbicg 
    "AY147208"      "411"        "1"       "11" 
    
    [[85]]
          name     length      frame     ncbicg 
    "AY147209"      "411"        "1"       "11" 
    
    [[86]]
          name     length      frame     ncbicg 
    "AY147210"      "411"        "1"       "11" 
    
    [[87]]
          name     length      frame     ncbicg 
    "AY147211"      "411"        "1"       "11" 
    
    [[88]]
          name     length      frame     ncbicg 
    "AY147212"      "411"        "1"       "11" 
    
    [[89]]
          name     length      frame     ncbicg 
    "AY147213"      "411"        "1"       "11" 
    
    [[90]]
          name     length      frame     ncbicg 
    "AY147214"      "411"        "1"       "11" 
    
    [[91]]
          name     length      frame     ncbicg 
    "AY147215"      "411"        "1"       "11" 
    
    [[92]]
          name     length      frame     ncbicg 
    "AY147216"      "411"        "1"       "11" 
    
    [[93]]
          name     length      frame     ncbicg 
    "AY147217"      "411"        "1"       "11" 
    
    [[94]]
          name     length      frame     ncbicg 
    "AY147218"      "411"        "1"       "11" 
    
    [[95]]
          name     length      frame     ncbicg 
    "AY155355"      "483"        "0"       "11" 
    
    [[96]]
          name     length      frame     ncbicg 
    "AY155356"      "483"        "0"       "11" 
    
    [[97]]
          name     length      frame     ncbicg 
    "AY155357"      "483"        "0"       "11" 
    
    [[98]]
          name     length      frame     ncbicg 
    "AY155358"      "483"        "0"       "11" 
    
    [[99]]
          name     length      frame     ncbicg 
    "AY155359"      "483"        "0"       "11" 
    
    [[100]]
          name     length      frame     ncbicg 
    "AY155360"      "483"        "0"       "11" 
    
    [[101]]
          name     length      frame     ncbicg 
    "AY155361"      "483"        "0"       "11" 
    
    [[102]]
          name     length      frame     ncbicg 
    "AY271363"      "159"        "0"       "11" 
    
    [[103]]
          name     length      frame     ncbicg 
    "AY271364"       "87"        "0"       "11" 
    
    [[104]]
          name     length      frame     ncbicg 
    "AY271365"      "244"        "0"       "11" 
    
    [[105]]
          name     length      frame     ncbicg 
    "AY271366"      "396"        "0"       "11" 
    
    [[106]]
          name     length      frame     ncbicg 
    "AY271367"      "291"        "0"       "11" 
    
    [[107]]
          name     length      frame     ncbicg 
    "AY271368"      "174"        "0"       "11" 
    
    [[108]]
          name     length      frame     ncbicg 
    "AY271369"      "357"        "0"       "11" 
    
    [[109]]
          name     length      frame     ncbicg 
    "AY271370"      "282"        "0"       "11" 
    
    [[110]]
          name     length      frame     ncbicg 
    "AY271371"      "318"        "0"       "11" 
    
    [[111]]
          name     length      frame     ncbicg 
    "AY271372"      "330"        "0"       "11" 
    
    [[112]]
          name     length      frame     ncbicg 
    "AY271373"      "306"        "0"       "11" 
    
    [[113]]
          name     length      frame     ncbicg 
    "AY271374"      "318"        "0"       "11" 
    
    [[114]]
          name     length      frame     ncbicg 
    "AY271375"      "318"        "0"       "11" 
    
    [[115]]
          name     length      frame     ncbicg 
    "AY271376"      "329"        "0"       "11" 
    
    [[116]]
          name     length      frame     ncbicg 
    "AY271377"      "282"        "0"       "11" 
    
    [[117]]
          name     length      frame     ncbicg 
    "AY271900"      "189"        "1"       "11" 
    
    [[118]]
          name     length      frame     ncbicg 
    "AY280807"      "342"        "0"       "11" 
    
    [[119]]
          name     length      frame     ncbicg 
    "AY280808"       "87"        "0"       "11" 
    
    [[120]]
          name     length      frame     ncbicg 
    "AY280809"       "87"        "0"       "11" 
    
    [[121]]
          name     length      frame     ncbicg 
    "AY280810"      "369"        "0"       "11" 
    
    [[122]]
          name     length      frame     ncbicg 
    "AY280811"      "228"        "0"       "11" 
    
    [[123]]
          name     length      frame     ncbicg 
    "AY280812"       "87"        "0"       "11" 
    
    [[124]]
          name     length      frame     ncbicg 
    "AY280813"      "318"        "0"       "11" 
    
    [[125]]
          name     length      frame     ncbicg 
    "AY280814"      "375"        "0"       "11" 
    
    [[126]]
          name     length      frame     ncbicg 
    "AY280815"       "87"        "0"       "11" 
    
    [[127]]
          name     length      frame     ncbicg 
    "AY280816"      "339"        "0"       "11" 
    
    [[128]]
          name     length      frame     ncbicg 
    "AY280817"      "303"        "0"       "11" 
    
    [[129]]
          name     length      frame     ncbicg 
    "AY280818"      "336"        "0"       "11" 
    
    [[130]]
          name     length      frame     ncbicg 
    "AY280819"      "168"        "0"       "11" 
    
    [[131]]
          name     length      frame     ncbicg 
    "AY280820"      "387"        "0"       "11" 
    
    [[132]]
          name     length      frame     ncbicg 
    "AY280821"      "309"        "0"       "11" 
    
    [[133]]
          name     length      frame     ncbicg 
    "AY280822"      "210"        "0"       "11" 
    
    [[134]]
          name     length      frame     ncbicg 
    "AY280823"      "201"        "0"       "11" 
    
    [[135]]
          name     length      frame     ncbicg 
    "AY280824"      "297"        "0"       "11" 
    
    [[136]]
          name     length      frame     ncbicg 
    "AY280825"      "291"        "0"       "11" 
    
    [[137]]
          name     length      frame     ncbicg 
    "AY280826"      "294"        "0"       "11" 
    
    [[138]]
          name     length      frame     ncbicg 
    "AY280827"      "378"        "0"       "11" 
    
    [[139]]
          name     length      frame     ncbicg 
    "AY280828"      "219"        "0"       "11" 
    
    [[140]]
          name     length      frame     ncbicg 
    "AY280829"      "366"        "0"       "11" 
    
    [[141]]
          name     length      frame     ncbicg 
    "AY280830"       "87"        "0"       "11" 
    
    [[142]]
          name     length      frame     ncbicg 
    "AY280831"      "327"        "0"       "11" 
    
    [[143]]
          name     length      frame     ncbicg 
    "AY280832"       "87"        "0"       "11" 
    
    [[144]]
          name     length      frame     ncbicg 
    "AY280833"      "360"        "0"       "11" 
    
    [[145]]
          name     length      frame     ncbicg 
    "AY280834"      "345"        "0"       "11" 
    
    [[146]]
          name     length      frame     ncbicg 
    "AY280835"      "324"        "0"       "11" 
    
    [[147]]
          name     length      frame     ncbicg 
    "AY280836"      "258"        "0"       "11" 
    
    [[148]]
          name     length      frame     ncbicg 
    "AY280837"       "87"        "0"       "11" 
    
    [[149]]
          name     length      frame     ncbicg 
    "AY280838"      "192"        "0"       "11" 
    
    [[150]]
          name     length      frame     ncbicg 
    "AY280839"      "297"        "0"       "11" 
    
    [[151]]
          name     length      frame     ncbicg 
    "AY280841"      "351"        "0"       "11" 
    
    [[152]]
          name     length      frame     ncbicg 
    "AY280842"      "339"        "0"       "11" 
    
    [[153]]
          name     length      frame     ncbicg 
    "AY280843"      "288"        "0"       "11" 
    
    [[154]]
          name     length      frame     ncbicg 
    "AY280844"      "360"        "0"       "11" 
    
    [[155]]
          name     length      frame     ncbicg 
    "AY280845"      "351"        "0"       "11" 
    
    [[156]]
          name     length      frame     ncbicg 
    "AY280846"      "306"        "0"       "11" 
    
    [[157]]
          name     length      frame     ncbicg 
    "AY308001"      "330"        "0"       "11" 
    
    [[158]]
          name     length      frame     ncbicg 
    "AY308002"      "378"        "0"       "11" 
    
    [[159]]
          name     length      frame     ncbicg 
    "AY308003"      "195"        "0"       "11" 
    
    [[160]]
          name     length      frame     ncbicg 
    "AY308004"      "294"        "0"       "11" 
    
    [[161]]
          name     length      frame     ncbicg 
    "AY308005"      "375"        "0"       "11" 
    
    [[162]]
          name     length      frame     ncbicg 
    "AY308006"      "390"        "0"       "11" 
    
    [[163]]
          name     length      frame     ncbicg 
    "AY308007"      "150"        "0"       "11" 
    
    [[164]]
          name     length      frame     ncbicg 
    "AY308008"      "198"        "0"       "11" 
    
    [[165]]
          name     length      frame     ncbicg 
    "AY308009"      "153"        "0"       "11" 
    
    [[166]]
          name     length      frame     ncbicg 
    "AY308010"      "294"        "0"       "11" 
    
    [[167]]
          name     length      frame     ncbicg 
    "AY308011"      "180"        "0"       "11" 
    
    [[168]]
          name     length      frame     ncbicg 
    "AY308012"      "135"        "0"       "11" 
    
    [[169]]
          name     length      frame     ncbicg 
    "AY308013"      "285"        "0"       "11" 
    
    [[170]]
          name     length      frame     ncbicg 
    "AY308014"      "177"        "0"       "11" 
    
    [[171]]
          name     length      frame     ncbicg 
    "AY308015"      "369"        "0"       "11" 
    
    [[172]]
          name     length      frame     ncbicg 
    "AY325125"      "159"        "0"       "11" 
    
    [[173]]
          name     length      frame     ncbicg 
    "AY325126"      "156"        "0"       "11" 
    
    [[174]]
          name     length      frame     ncbicg 
    "AY325127"      "159"        "0"       "11" 
    
    [[175]]
          name     length      frame     ncbicg 
    "AY325128"      "147"        "0"       "11" 
    
    [[176]]
          name     length      frame     ncbicg 
    "AY544973"      "398"        "0"       "11" 
    
    [[177]]
          name     length      frame     ncbicg 
    "AY544974"      "398"        "0"       "11" 
    
    [[178]]
          name     length      frame     ncbicg 
    "AY587520"      "248"        "0"       "11" 
    
    [[179]]
          name     length      frame     ncbicg 
    "AY663788"      "198"        "0"       "11" 
    
    [[180]]
          name     length      frame     ncbicg 
    "AY787173"      "367"        "2"       "11" 
    
    [[181]]
          name     length      frame     ncbicg 
    "AY793000"       "81"        "0"       "11" 
    
    [[182]]
          name     length      frame     ncbicg 
    "AY793001"       "81"        "0"       "11" 
    
    [[183]]
          name     length      frame     ncbicg 
    "AY793002"       "81"        "0"       "11" 
    
    [[184]]
          name     length      frame     ncbicg 
    "AY793003"       "81"        "0"       "11" 
    
    [[185]]
          name     length      frame     ncbicg 
    "AY793004"       "81"        "0"       "11" 
    
    [[186]]
          name     length      frame     ncbicg 
    "AY793005"       "84"        "0"       "11" 
    
    [[187]]
          name     length      frame     ncbicg 
    "AY793006"       "81"        "0"       "11" 
    
    [[188]]
          name     length      frame     ncbicg 
    "AY793007"       "81"        "0"       "11" 
    
    [[189]]
          name     length      frame     ncbicg 
    "AY793008"       "81"        "0"       "11" 
    
    [[190]]
          name     length      frame     ncbicg 
    "AY793009"       "81"        "0"       "11" 
    
    [[191]]
          name     length      frame     ncbicg 
    "AY793010"       "81"        "0"       "11" 
    
    [[192]]
          name     length      frame     ncbicg 
    "AY793011"       "81"        "0"       "11" 
    
    [[193]]
          name     length      frame     ncbicg 
    "AY793012"       "81"        "0"       "11" 
    
    [[194]]
          name     length      frame     ncbicg 
    "AY793013"       "81"        "0"       "11" 
    
    [[195]]
          name     length      frame     ncbicg 
    "AY793014"       "81"        "0"       "11" 
    
    [[196]]
          name     length      frame     ncbicg 
    "AY793015"       "81"        "0"       "11" 
    
    [[197]]
          name     length      frame     ncbicg 
    "AY793016"       "81"        "0"       "11" 
    
    [[198]]
          name     length      frame     ncbicg 
    "AY793017"       "81"        "0"       "11" 
    
    [[199]]
          name     length      frame     ncbicg 
    "AY793018"       "81"        "0"       "11" 
    
    [[200]]
          name     length      frame     ncbicg 
    "AY793019"       "81"        "0"       "11" 
    
    [[201]]
          name     length      frame     ncbicg 
    "AY800477"       "81"        "0"       "11" 
    
    [[202]]
          name     length      frame     ncbicg 
    "AY800478"       "81"        "0"       "11" 
    
    [[203]]
          name     length      frame     ncbicg 
    "AY800479"       "81"        "0"       "11" 
    
    [[204]]
          name     length      frame     ncbicg 
    "AY800480"       "81"        "0"       "11" 
    
    [[205]]
          name     length      frame     ncbicg 
    "AY800481"       "81"        "0"       "11" 
    
    [[206]]
          name     length      frame     ncbicg 
    "AY800482"       "81"        "0"       "11" 
    
    [[207]]
          name     length      frame     ncbicg 
    "AY800483"       "81"        "0"       "11" 
    
    [[208]]
          name     length      frame     ncbicg 
    "AY800484"       "81"        "0"       "11" 
    
    [[209]]
          name     length      frame     ncbicg 
    "AY800485"       "81"        "0"       "11" 
    
    [[210]]
          name     length      frame     ncbicg 
    "AY800486"       "81"        "0"       "11" 
    
    [[211]]
          name     length      frame     ncbicg 
    "AY800487"       "81"        "0"       "11" 
    
    [[212]]
          name     length      frame     ncbicg 
    "AY800488"       "81"        "0"       "11" 
    
    [[213]]
          name     length      frame     ncbicg 
    "AY800489"       "84"        "0"       "11" 
    
    [[214]]
          name     length      frame     ncbicg 
    "AY800490"       "81"        "0"       "11" 
    
    [[215]]
          name     length      frame     ncbicg 
    "AY800491"       "81"        "0"       "11" 
    
    [[216]]
          name     length      frame     ncbicg 
    "AY800492"       "81"        "0"       "11" 
    
    [[217]]
          name     length      frame     ncbicg 
    "AY800493"       "81"        "0"       "11" 
    
    [[218]]
          name     length      frame     ncbicg 
    "AY800494"       "81"        "0"       "11" 
    
    [[219]]
          name     length      frame     ncbicg 
    "AY800495"       "81"        "0"       "11" 
    
    [[220]]
          name     length      frame     ncbicg 
    "AY800496"       "81"        "0"       "11" 
    
    [[221]]
          name     length      frame     ncbicg 
    "AY800497"       "81"        "0"       "11" 
    
    [[222]]
          name     length      frame     ncbicg 
    "AY800498"       "81"        "0"       "11" 
    
    [[223]]
          name     length      frame     ncbicg 
    "AY800499"       "81"        "0"       "11" 
    
    [[224]]
          name     length      frame     ncbicg 
    "AY800500"       "81"        "0"       "11" 
    
    [[225]]
          name     length      frame     ncbicg 
    "AY800501"       "81"        "0"       "11" 
    
    [[226]]
          name     length      frame     ncbicg 
    "AY800502"       "81"        "0"       "11" 
    
    [[227]]
          name     length      frame     ncbicg 
    "AY800503"       "81"        "0"       "11" 
    
    [[228]]
          name     length      frame     ncbicg 
    "AY800504"       "81"        "0"       "11" 
    
    [[229]]
          name     length      frame     ncbicg 
    "AY800505"       "81"        "0"       "11" 
    
    [[230]]
          name     length      frame     ncbicg 
    "AY800506"       "81"        "0"       "11" 
    
    [[231]]
          name     length      frame     ncbicg 
    "AY800507"       "81"        "0"       "11" 
    
    [[232]]
          name     length      frame     ncbicg 
    "AY800508"       "81"        "0"       "11" 
    
    [[233]]
          name     length      frame     ncbicg 
    "AY800509"       "81"        "0"       "11" 
    
    [[234]]
          name     length      frame     ncbicg 
    "AY800510"       "81"        "0"       "11" 
    
    [[235]]
          name     length      frame     ncbicg 
    "AY800511"       "81"        "0"       "11" 
    
    [[236]]
          name     length      frame     ncbicg 
    "AY800512"       "81"        "0"       "11" 
    
    [[237]]
          name     length      frame     ncbicg 
    "AY800513"       "81"        "0"       "11" 
    
    [[238]]
          name     length      frame     ncbicg 
    "AY800514"       "81"        "0"       "11" 
    
    [[239]]
          name     length      frame     ncbicg 
    "AY800515"       "81"        "0"       "11" 
    
    [[240]]
          name     length      frame     ncbicg 
    "AY800516"       "81"        "0"       "11" 
    
    [[241]]
          name     length      frame     ncbicg 
    "AY800517"       "81"        "0"       "11" 
    
    [[242]]
          name     length      frame     ncbicg 
    "AY800518"       "81"        "0"       "11" 
    
    [[243]]
          name     length      frame     ncbicg 
    "AY800519"       "81"        "0"       "11" 
    
    [[244]]
          name     length      frame     ncbicg 
    "AY800520"       "81"        "0"       "11" 
    
    [[245]]
          name     length      frame     ncbicg 
    "AY800521"       "81"        "0"       "11" 
    
    [[246]]
          name     length      frame     ncbicg 
    "AY800522"       "81"        "0"       "11" 
    
    [[247]]
          name     length      frame     ncbicg 
    "AY800523"       "81"        "0"       "11" 
    
    [[248]]
          name     length      frame     ncbicg 
    "AY800524"       "81"        "0"       "11" 
    
    [[249]]
          name     length      frame     ncbicg 
    "AY800525"       "81"        "0"       "11" 
    
    [[250]]
          name     length      frame     ncbicg 
    "AY800526"       "81"        "0"       "11" 
    
    [[251]]
          name     length      frame     ncbicg 
    "AY800527"       "81"        "0"       "11" 
    
    [[252]]
          name     length      frame     ncbicg 
    "AY800528"       "81"        "0"       "11" 
    
    [[253]]
          name     length      frame     ncbicg 
    "AY800529"       "81"        "0"       "11" 
    
    [[254]]
          name     length      frame     ncbicg 
    "AY800530"       "81"        "0"       "11" 
    
    [[255]]
          name     length      frame     ncbicg 
    "AY800531"       "81"        "0"       "11" 
    
    [[256]]
          name     length      frame     ncbicg 
    "AY800532"       "81"        "0"       "11" 
    
    [[257]]
          name     length      frame     ncbicg 
    "AY800533"       "81"        "0"       "11" 
    
    [[258]]
          name     length      frame     ncbicg 
    "AY800534"       "81"        "0"       "11" 
    
    [[259]]
          name     length      frame     ncbicg 
    "AY800535"       "81"        "0"       "11" 
    
    [[260]]
          name     length      frame     ncbicg 
    "AY800536"       "81"        "0"       "11" 
    
    [[261]]
          name     length      frame     ncbicg 
    "AY800537"       "81"        "0"       "11" 
    
    [[262]]
          name     length      frame     ncbicg 
    "AY800538"       "81"        "0"       "11" 
    
    [[263]]
          name     length      frame     ncbicg 
    "AY800539"       "81"        "0"       "11" 
    
    [[264]]
          name     length      frame     ncbicg 
    "AY800540"       "81"        "0"       "11" 
    
    [[265]]
          name     length      frame     ncbicg 
    "AY800541"       "81"        "0"       "11" 
    
    [[266]]
          name     length      frame     ncbicg 
    "AY800542"       "81"        "0"       "11" 
    
    [[267]]
          name     length      frame     ncbicg 
    "AY800543"       "81"        "0"       "11" 
    
    [[268]]
          name     length      frame     ncbicg 
    "AY800544"       "81"        "0"       "11" 
    
    [[269]]
          name     length      frame     ncbicg 
    "AY800545"       "81"        "0"       "11" 
    
    [[270]]
          name     length      frame     ncbicg 
    "AY800546"       "81"        "0"       "11" 
    
    [[271]]
          name     length      frame     ncbicg 
    "AY800547"       "81"        "0"       "11" 
    
    [[272]]
          name     length      frame     ncbicg 
    "AY800548"       "81"        "0"       "11" 
    
    [[273]]
          name     length      frame     ncbicg 
    "AY800549"       "81"        "0"       "11" 
    
    [[274]]
          name     length      frame     ncbicg 
    "AY800550"       "81"        "0"       "11" 
    
    [[275]]
          name     length      frame     ncbicg 
    "AY800551"       "81"        "0"       "11" 
    
    [[276]]
          name     length      frame     ncbicg 
    "AY800552"       "81"        "0"       "11" 
    
    [[277]]
          name     length      frame     ncbicg 
    "AY800553"       "81"        "0"       "11" 
    
    [[278]]
          name     length      frame     ncbicg 
    "AY800554"       "81"        "0"       "11" 
    
    [[279]]
          name     length      frame     ncbicg 
    "AY800555"       "81"        "0"       "11" 
    
    [[280]]
          name     length      frame     ncbicg 
    "AY800556"       "84"        "0"       "11" 
    
    [[281]]
          name     length      frame     ncbicg 
    "AY819713"      "313"        "0"       "11" 
    
    [[282]]
          name     length      frame     ncbicg 
    "AY819714"      "303"        "0"       "11" 
    
    [[283]]
          name     length      frame     ncbicg 
    "AY823310"      "385"        "2"       "11" 
    
    [[284]]
          name     length      frame     ncbicg 
    "AY823311"      "385"        "2"       "11" 
    
    [[285]]
          name     length      frame     ncbicg 
    "AY823312"      "385"        "2"       "11" 
    
    [[286]]
          name     length      frame     ncbicg 
    "AY823313"      "376"        "2"       "11" 
    
    [[287]]
          name     length      frame     ncbicg 
    "AY823314"      "385"        "2"       "11" 
    
    [[288]]
          name     length      frame     ncbicg 
    "AY823315"      "385"        "2"       "11" 
    
    [[289]]
          name     length      frame     ncbicg 
    "AY823316"      "385"        "2"       "11" 
    
    [[290]]
          name     length      frame     ncbicg 
    "AY823317"      "385"        "2"       "11" 
    
    [[291]]
          name     length      frame     ncbicg 
    "AY823318"      "385"        "2"       "11" 
    
    [[292]]
          name     length      frame     ncbicg 
    "AY898740"      "295"        "0"       "11" 
    
    [[293]]
          name     length      frame     ncbicg 
    "AY898744"      "295"        "0"       "11" 
    
    [[294]]
               name          length           frame          ncbicg 
    "CP000611.RPOB"          "3519"             "0"            "11" 
    
    [[295]]
               name          length           frame          ncbicg 
    "CP002871.RPOB"          "3519"             "0"            "11" 
    
    [[296]]
               name          length           frame          ncbicg 
    "CP002882.RPOB"          "3519"             "0"            "11" 
    
    [[297]]
               name          length           frame          ncbicg 
    "CP002883.RPOB"          "3519"             "0"            "11" 
    
    [[298]]
               name          length           frame          ncbicg 
    "CP002884.RPOB"          "3519"             "0"            "11" 
    
    [[299]]
               name          length           frame          ncbicg 
    "CP002885.RPOB"          "3519"             "0"            "11" 
    
    [[300]]
               name          length           frame          ncbicg 
    "CP002992.RPOB"          "3519"             "0"            "11" 
    
    [[301]]
               name          length           frame          ncbicg 
    "CP003234.RPOB"           "129"             "0"            "11" 
    
    [[302]]
               name          length           frame          ncbicg 
    "CP005082.RPOB"          "3519"             "0"            "11" 
    
    [[303]]
               name          length           frame          ncbicg 
    "CP005386.RPOB"          "3516"             "0"            "11" 
    
    [[304]]
               name          length           frame          ncbicg 
    "CP005387.RPOB"          "3519"             "0"            "11" 
    
    [[305]]
               name          length           frame          ncbicg 
    "CP007027.RPOB"          "3519"             "0"            "11" 
    
    [[306]]
               name          length           frame          ncbicg 
    "CP007299.RPOB"          "3519"             "0"            "11" 
    
    [[307]]
                name           length            frame           ncbicg 
    "CP012506.PE697"           "3519"              "0"             "11" 
    
    [[308]]
               name          length           frame          ncbicg 
    "CP013475.RPOB"          "3291"             "0"            "11" 
    
    [[309]]
                name           length            frame           ncbicg 
    "CP016794.PE693"           "3291"              "0"             "11" 
    
    [[310]]
                name           length            frame           ncbicg 
    "CP016888.PE692"           "3291"              "0"             "11" 
    
    [[311]]
               name          length           frame          ncbicg 
    "CP024614.RPOB"          "1557"             "0"            "11" 
    
    [[312]]
          name     length      frame     ncbicg 
    "DQ118534"       "87"        "0"       "11" 
    
    [[313]]
          name     length      frame     ncbicg 
    "DQ205438"      "261"        "0"       "11" 
    
    [[314]]
          name     length      frame     ncbicg 
    "DQ205439"      "261"        "0"       "11" 
    
    [[315]]
          name     length      frame     ncbicg 
    "DQ205440"      "261"        "0"       "11" 
    
    [[316]]
          name     length      frame     ncbicg 
    "DQ205441"      "261"        "0"       "11" 
    
    [[317]]
          name     length      frame     ncbicg 
    "DQ985197"      "780"        "0"       "11" 
    
    [[318]]
          name     length      frame     ncbicg 
    "DQ985198"      "780"        "0"       "11" 
    
    [[319]]
          name     length      frame     ncbicg 
    "DQ985199"      "780"        "0"       "11" 
    
    [[320]]
          name     length      frame     ncbicg 
    "DQ985200"      "780"        "0"       "11" 
    
    [[321]]
          name     length      frame     ncbicg 
    "DQ985201"      "780"        "0"       "11" 
    
    [[322]]
          name     length      frame     ncbicg 
    "DQ985202"      "780"        "0"       "11" 
    
    [[323]]
          name     length      frame     ncbicg 
    "DQ985203"      "780"        "0"       "11" 
    
    [[324]]
          name     length      frame     ncbicg 
    "DQ985204"      "780"        "0"       "11" 
    
    [[325]]
          name     length      frame     ncbicg 
    "DQ985205"      "780"        "0"       "11" 
    
    [[326]]
          name     length      frame     ncbicg 
    "DQ985206"      "780"        "0"       "11" 
    
    [[327]]
          name     length      frame     ncbicg 
    "DQ985207"      "780"        "0"       "11" 
    
    [[328]]
          name     length      frame     ncbicg 
    "DQ985208"      "780"        "0"       "11" 
    
    [[329]]
          name     length      frame     ncbicg 
    "DQ985209"      "780"        "0"       "11" 
    
    [[330]]
          name     length      frame     ncbicg 
    "DQ985210"      "780"        "0"       "11" 
    
    [[331]]
          name     length      frame     ncbicg 
    "DQ985211"      "780"        "0"       "11" 
    
    [[332]]
          name     length      frame     ncbicg 
    "DQ985212"      "780"        "0"       "11" 
    
    [[333]]
          name     length      frame     ncbicg 
    "DQ985213"      "780"        "0"       "11" 
    
    [[334]]
          name     length      frame     ncbicg 
    "DQ985214"      "780"        "0"       "11" 
    
    [[335]]
          name     length      frame     ncbicg 
    "DQ985215"      "780"        "0"       "11" 
    
    [[336]]
          name     length      frame     ncbicg 
    "DQ985216"      "780"        "0"       "11" 
    
    [[337]]
          name     length      frame     ncbicg 
    "DQ985217"      "780"        "0"       "11" 
    
    [[338]]
          name     length      frame     ncbicg 
    "DQ985218"      "780"        "0"       "11" 
    
    [[339]]
          name     length      frame     ncbicg 
    "DQ985219"      "780"        "0"       "11" 
    
    [[340]]
          name     length      frame     ncbicg 
    "DQ985220"      "780"        "0"       "11" 
    
    [[341]]
          name     length      frame     ncbicg 
    "DQ985221"      "780"        "0"       "11" 
    
    [[342]]
          name     length      frame     ncbicg 
    "DQ985222"      "780"        "0"       "11" 
    
    [[343]]
          name     length      frame     ncbicg 
    "DQ985223"      "780"        "0"       "11" 
    
    [[344]]
          name     length      frame     ncbicg 
    "DQ985224"      "780"        "0"       "11" 
    
    [[345]]
          name     length      frame     ncbicg 
    "DQ985225"      "780"        "0"       "11" 
    
    [[346]]
          name     length      frame     ncbicg 
    "DQ985226"      "780"        "0"       "11" 
    
    [[347]]
          name     length      frame     ncbicg 
    "DQ985227"      "780"        "0"       "11" 
    
    [[348]]
          name     length      frame     ncbicg 
    "DQ985228"      "780"        "0"       "11" 
    
    [[349]]
          name     length      frame     ncbicg 
    "EF064790"       "80"        "1"       "11" 
    
    [[350]]
          name     length      frame     ncbicg 
    "EF628293"      "304"        "0"       "11" 
    
    [[351]]
          name     length      frame     ncbicg 
    "EF628300"      "386"        "2"       "11" 
    
    [[352]]
          name     length      frame     ncbicg 
    "EF628301"      "379"        "1"       "11" 
    
    [[353]]
          name     length      frame     ncbicg 
    "EF628302"      "385"        "2"       "11" 
    
    [[354]]
          name     length      frame     ncbicg 
    "EF628303"      "383"        "2"       "11" 
    
    [[355]]
          name     length      frame     ncbicg 
    "EF628304"      "372"        "2"       "11" 
    
    [[356]]
          name     length      frame     ncbicg 
    "EF628305"      "400"        "2"       "11" 
    
    [[357]]
          name     length      frame     ncbicg 
    "EF628306"      "393"        "2"       "11" 
    
    [[358]]
          name     length      frame     ncbicg 
    "EF628307"      "395"        "2"       "11" 
    
    [[359]]
          name     length      frame     ncbicg 
    "EF628308"      "366"        "2"       "11" 
    
    [[360]]
          name     length      frame     ncbicg 
    "EF628309"      "384"        "2"       "11" 
    
    [[361]]
          name     length      frame     ncbicg 
    "EF628310"      "367"        "2"       "11" 
    
    [[362]]
          name     length      frame     ncbicg 
    "EF628311"      "351"        "0"       "11" 
    
    [[363]]
          name     length      frame     ncbicg 
    "EF628312"      "357"        "0"       "11" 
    
    [[364]]
          name     length      frame     ncbicg 
    "EF628313"      "404"        "2"       "11" 
    
    [[365]]
          name     length      frame     ncbicg 
    "EF628314"      "372"        "0"       "11" 
    
    [[366]]
          name     length      frame     ncbicg 
    "EF628315"      "358"        "2"       "11" 
    
    [[367]]
          name     length      frame     ncbicg 
    "EF628316"      "373"        "0"       "11" 
    
    [[368]]
          name     length      frame     ncbicg 
    "EF628317"      "338"        "2"       "11" 
    
    [[369]]
          name     length      frame     ncbicg 
    "EF628318"      "365"        "1"       "11" 
    
    [[370]]
          name     length      frame     ncbicg 
    "EF628319"      "369"        "0"       "11" 
    
    [[371]]
          name     length      frame     ncbicg 
    "EF628320"      "384"        "2"       "11" 
    
    [[372]]
          name     length      frame     ncbicg 
    "EF628321"      "379"        "1"       "11" 
    
    [[373]]
          name     length      frame     ncbicg 
    "EF628322"      "379"        "1"       "11" 
    
    [[374]]
          name     length      frame     ncbicg 
    "EF628323"      "367"        "2"       "11" 
    
    [[375]]
          name     length      frame     ncbicg 
    "EF628324"      "372"        "2"       "11" 
    
    [[376]]
          name     length      frame     ncbicg 
    "EF628325"      "349"        "1"       "11" 
    
    [[377]]
          name     length      frame     ncbicg 
    "EF628326"      "380"        "2"       "11" 
    
    [[378]]
          name     length      frame     ncbicg 
    "EF628327"      "381"        "0"       "11" 
    
    [[379]]
          name     length      frame     ncbicg 
    "EF628328"      "381"        "1"       "11" 
    
    [[380]]
          name     length      frame     ncbicg 
    "EF628332"      "403"        "2"       "11" 
    
    [[381]]
          name     length      frame     ncbicg 
    "EF628333"      "346"        "0"       "11" 
    
    [[382]]
          name     length      frame     ncbicg 
    "EF628334"      "390"        "0"       "11" 
    
    [[383]]
          name     length      frame     ncbicg 
    "EF628335"      "369"        "1"       "11" 
    
    [[384]]
          name     length      frame     ncbicg 
    "EF628337"      "373"        "1"       "11" 
    
    [[385]]
          name     length      frame     ncbicg 
    "EF628338"      "311"        "1"       "11" 
    
    [[386]]
          name     length      frame     ncbicg 
    "EF628341"      "359"        "2"       "11" 
    
    [[387]]
          name     length      frame     ncbicg 
    "EF628345"      "331"        "0"       "11" 
    
    [[388]]
          name     length      frame     ncbicg 
    "EF628346"      "348"        "1"       "11" 
    
    [[389]]
          name     length      frame     ncbicg 
    "EF628348"      "335"        "1"       "11" 
    
    [[390]]
          name     length      frame     ncbicg 
    "EF628351"      "397"        "0"       "11" 
    
    [[391]]
          name     length      frame     ncbicg 
    "EF628353"      "384"        "2"       "11" 
    
    [[392]]
          name     length      frame     ncbicg 
    "EF628354"      "391"        "2"       "11" 
    
    [[393]]
          name     length      frame     ncbicg 
    "EF628355"      "340"        "2"       "11" 
    
    [[394]]
          name     length      frame     ncbicg 
    "EF628356"      "384"        "2"       "11" 
    
    [[395]]
          name     length      frame     ncbicg 
    "EF628357"      "338"        "2"       "11" 
    
    [[396]]
          name     length      frame     ncbicg 
    "EF628359"      "331"        "2"       "11" 
    
    [[397]]
          name     length      frame     ncbicg 
    "EF628360"      "398"        "2"       "11" 
    
    [[398]]
          name     length      frame     ncbicg 
    "EF628364"      "385"        "2"       "11" 
    
    [[399]]
          name     length      frame     ncbicg 
    "EF628366"      "324"        "2"       "11" 
    
    [[400]]
          name     length      frame     ncbicg 
    "EF628368"      "315"        "0"       "11" 
    
    [[401]]
          name     length      frame     ncbicg 
    "EF661663"      "216"        "0"       "11" 
    
    [[402]]
          name     length      frame     ncbicg 
    "EU325648"      "252"        "0"       "11" 
    
    [[403]]
          name     length      frame     ncbicg 
    "FJ915183"      "184"        "0"       "11" 
    
    [[404]]
          name     length      frame     ncbicg 
    "FJ915184"      "193"        "0"       "11" 
    
    [[405]]
          name     length      frame     ncbicg 
    "FJ915185"      "193"        "0"       "11" 
    
    [[406]]
          name     length      frame     ncbicg 
    "FJ915186"      "184"        "0"       "11" 
    
    [[407]]
          name     length      frame     ncbicg 
    "FJ915187"      "184"        "0"       "11" 
    
    [[408]]
          name     length      frame     ncbicg 
    "FJ915188"      "189"        "0"       "11" 
    
    [[409]]
          name     length      frame     ncbicg 
    "FM172959"       "84"        "0"       "11" 
    
    [[410]]
          name     length      frame     ncbicg 
    "FM172960"       "84"        "0"       "11" 
    
    [[411]]
          name     length      frame     ncbicg 
    "FM172961"       "84"        "0"       "11" 
    
    [[412]]
          name     length      frame     ncbicg 
    "FM172962"       "84"        "0"       "11" 
    
    [[413]]
          name     length      frame     ncbicg 
    "FM172963"       "84"        "0"       "11" 
    
    [[414]]
          name     length      frame     ncbicg 
    "FM172964"       "84"        "0"       "11" 
    
    [[415]]
          name     length      frame     ncbicg 
    "FM172965"       "84"        "0"       "11" 
    
    [[416]]
          name     length      frame     ncbicg 
    "FM172966"       "84"        "0"       "11" 
    
    [[417]]
          name     length      frame     ncbicg 
    "FM172967"       "84"        "0"       "11" 
    
    [[418]]
          name     length      frame     ncbicg 
    "FM172968"       "84"        "0"       "11" 
    
    [[419]]
          name     length      frame     ncbicg 
    "FM172969"       "84"        "0"       "11" 
    
    [[420]]
          name     length      frame     ncbicg 
    "FM172970"       "84"        "0"       "11" 
    
    [[421]]
          name     length      frame     ncbicg 
    "FM172971"       "84"        "0"       "11" 
    
    [[422]]
          name     length      frame     ncbicg 
    "FM172972"       "84"        "0"       "11" 
    
    [[423]]
          name     length      frame     ncbicg 
    "FM172973"       "84"        "0"       "11" 
    
    [[424]]
          name     length      frame     ncbicg 
    "FM172974"       "84"        "0"       "11" 
    
    [[425]]
          name     length      frame     ncbicg 
    "FM172975"       "84"        "0"       "11" 
    
    [[426]]
          name     length      frame     ncbicg 
    "FM172976"       "84"        "0"       "11" 
    
    [[427]]
          name     length      frame     ncbicg 
    "FM172977"      "525"        "0"       "11" 
    
    [[428]]
          name     length      frame     ncbicg 
    "FM172978"      "525"        "0"       "11" 
    
    [[429]]
          name     length      frame     ncbicg 
    "FM172992"      "525"        "0"       "11" 
    
    [[430]]
          name     length      frame     ncbicg 
    "FM172993"      "525"        "0"       "11" 
    
    [[431]]
          name     length      frame     ncbicg 
    "GQ250580"      "513"        "0"       "11" 
    
    [[432]]
          name     length      frame     ncbicg 
    "GQ250581"      "513"        "0"       "11" 
    
    [[433]]
          name     length      frame     ncbicg 
    "GQ293224"      "348"        "0"       "11" 
    
    [[434]]
          name     length      frame     ncbicg 
    "GQ395623"      "345"        "1"       "11" 
    
    [[435]]
          name     length      frame     ncbicg 
    "GQ871909"      "114"        "0"       "11" 
    
    [[436]]
          name     length      frame     ncbicg 
    "GQ871910"      "114"        "0"       "11" 
    
    [[437]]
          name     length      frame     ncbicg 
    "GQ871911"      "114"        "0"       "11" 
    
    [[438]]
          name     length      frame     ncbicg 
    "GQ871912"      "114"        "0"       "11" 
    
    [[439]]
          name     length      frame     ncbicg 
    "GQ871913"      "114"        "0"       "11" 
    
    [[440]]
          name     length      frame     ncbicg 
    "GQ871914"      "123"        "0"       "11" 
    
    [[441]]
          name     length      frame     ncbicg 
    "GQ871915"      "114"        "0"       "11" 
    
    [[442]]
          name     length      frame     ncbicg 
    "GU904010"      "179"        "0"       "11" 
    
    [[443]]
          name     length      frame     ncbicg 
    "GU904011"      "179"        "0"       "11" 
    
    [[444]]
          name     length      frame     ncbicg 
    "GU904012"      "179"        "0"       "11" 
    
    [[445]]
          name     length      frame     ncbicg 
    "GU904013"      "179"        "0"       "11" 
    
    [[446]]
          name     length      frame     ncbicg 
    "GU904014"      "179"        "0"       "11" 
    
    [[447]]
          name     length      frame     ncbicg 
    "GU904015"      "179"        "0"       "11" 
    
    [[448]]
          name     length      frame     ncbicg 
    "GU904016"      "179"        "0"       "11" 
    
    [[449]]
          name     length      frame     ncbicg 
    "GU904017"      "179"        "0"       "11" 
    
    [[450]]
          name     length      frame     ncbicg 
    "GU904018"      "179"        "0"       "11" 
    
    [[451]]
          name     length      frame     ncbicg 
    "GU904019"      "179"        "0"       "11" 
    
    [[452]]
          name     length      frame     ncbicg 
    "GU904020"      "179"        "0"       "11" 
    
    [[453]]
          name     length      frame     ncbicg 
    "GU904021"      "179"        "0"       "11" 
    
    [[454]]
          name     length      frame     ncbicg 
    "GU904022"      "179"        "0"       "11" 
    
    [[455]]
          name     length      frame     ncbicg 
    "GU904023"      "179"        "0"       "11" 
    
    [[456]]
          name     length      frame     ncbicg 
    "GU904024"      "179"        "0"       "11" 
    
    [[457]]
               name          length           frame          ncbicg 
    "HE608151.RPOB"          "3519"             "0"            "11" 
    
    [[458]]
               name          length           frame          ncbicg 
    "HG813240.RPOB"          "3519"             "0"            "11" 
    
    [[459]]
          name     length      frame     ncbicg 
    "HM048901"      "305"        "0"       "11" 
    
    [[460]]
          name     length      frame     ncbicg 
    "HM048902"      "305"        "0"       "11" 
    
    [[461]]
          name     length      frame     ncbicg 
    "HM048903"      "305"        "0"       "11" 
    
    [[462]]
          name     length      frame     ncbicg 
    "HM048904"      "305"        "0"       "11" 
    
    [[463]]
          name     length      frame     ncbicg 
    "HM215248"       "84"        "0"       "11" 
    
    [[464]]
          name     length      frame     ncbicg 
    "HM229777"      "344"        "1"       "11" 
    
    [[465]]
          name     length      frame     ncbicg 
    "HM345980"      "559"        "1"       "11" 
    
    [[466]]
          name     length      frame     ncbicg 
    "HM355826"      "371"        "0"       "11" 
    
    [[467]]
          name     length      frame     ncbicg 
    "HM355827"      "360"        "0"       "11" 
    
    [[468]]
          name     length      frame     ncbicg 
    "HM776946"      "271"        "0"       "11" 
    
    [[469]]
          name     length      frame     ncbicg 
    "HM776950"      "271"        "0"       "11" 
    
    [[470]]
          name     length      frame     ncbicg 
    "HM776953"      "271"        "0"       "11" 
    
    [[471]]
          name     length      frame     ncbicg 
    "HM776954"      "271"        "0"       "11" 
    
    [[472]]
          name     length      frame     ncbicg 
    "HM776955"      "271"        "0"       "11" 
    
    [[473]]
          name     length      frame     ncbicg 
    "HM776956"      "271"        "0"       "11" 
    
    [[474]]
          name     length      frame     ncbicg 
    "HM776958"      "231"        "0"       "11" 
    
    [[475]]
          name     length      frame     ncbicg 
    "HQ286613"      "192"        "0"       "11" 
    
    [[476]]
          name     length      frame     ncbicg 
    "HQ286614"      "192"        "0"       "11" 
    
    [[477]]
          name     length      frame     ncbicg 
    "HQ286615"      "192"        "0"       "11" 
    
    [[478]]
          name     length      frame     ncbicg 
    "HQ286616"      "192"        "0"       "11" 
    
    [[479]]
          name     length      frame     ncbicg 
    "HQ286617"      "192"        "0"       "11" 
    
    [[480]]
          name     length      frame     ncbicg 
    "HQ286618"      "192"        "0"       "11" 
    
    [[481]]
          name     length      frame     ncbicg 
    "HQ286619"      "351"        "0"       "11" 
    
    [[482]]
          name     length      frame     ncbicg 
    "HQ286620"      "351"        "0"       "11" 
    
    [[483]]
          name     length      frame     ncbicg 
    "HQ286621"      "351"        "0"       "11" 
    
    [[484]]
          name     length      frame     ncbicg 
    "HQ286622"      "351"        "0"       "11" 
    
    [[485]]
          name     length      frame     ncbicg 
    "HQ286623"      "351"        "0"       "11" 
    
    [[486]]
          name     length      frame     ncbicg 
    "HQ286624"      "351"        "0"       "11" 
    
    [[487]]
          name     length      frame     ncbicg 
    "HQ286625"      "351"        "0"       "11" 
    
    [[488]]
          name     length      frame     ncbicg 
    "HQ286626"      "351"        "0"       "11" 
    
    [[489]]
          name     length      frame     ncbicg 
    "HQ286627"      "351"        "0"       "11" 
    
    [[490]]
          name     length      frame     ncbicg 
    "HQ286628"      "351"        "0"       "11" 
    
    [[491]]
          name     length      frame     ncbicg 
    "HQ377336"      "332"        "1"       "11" 
    
    [[492]]
          name     length      frame     ncbicg 
    "HQ377337"      "332"        "1"       "11" 
    
    [[493]]
          name     length      frame     ncbicg 
    "HQ377338"      "331"        "0"       "11" 
    
    [[494]]
          name     length      frame     ncbicg 
    "HQ377339"      "331"        "1"       "11" 
    
    [[495]]
          name     length      frame     ncbicg 
    "HQ377340"      "332"        "1"       "11" 
    
    [[496]]
          name     length      frame     ncbicg 
    "HQ377341"      "332"        "1"       "11" 
    
    [[497]]
          name     length      frame     ncbicg 
    "HQ377342"      "332"        "1"       "11" 
    
    [[498]]
          name     length      frame     ncbicg 
    "HQ377343"      "332"        "1"       "11" 
    
    [[499]]
          name     length      frame     ncbicg 
    "HQ377344"      "332"        "1"       "11" 
    
    [[500]]
          name     length      frame     ncbicg 
    "HQ377345"      "332"        "1"       "11" 
    
    [[501]]
          name     length      frame     ncbicg 
    "HQ377346"      "332"        "1"       "11" 
    
    [[502]]
          name     length      frame     ncbicg 
    "HQ377347"      "332"        "1"       "11" 
    
    [[503]]
          name     length      frame     ncbicg 
    "HQ377348"      "332"        "1"       "11" 
    
    [[504]]
          name     length      frame     ncbicg 
    "HQ377349"      "332"        "1"       "11" 
    
    [[505]]
          name     length      frame     ncbicg 
    "HQ377350"      "332"        "1"       "11" 
    
    [[506]]
          name     length      frame     ncbicg 
    "HQ377351"      "332"        "1"       "11" 
    
    [[507]]
          name     length      frame     ncbicg 
    "HQ377352"      "332"        "1"       "11" 
    
    [[508]]
          name     length      frame     ncbicg 
    "HQ377353"      "332"        "1"       "11" 
    
    [[509]]
          name     length      frame     ncbicg 
    "HQ377354"      "332"        "1"       "11" 
    
    [[510]]
          name     length      frame     ncbicg 
    "HQ377355"      "332"        "1"       "11" 
    
    [[511]]
          name     length      frame     ncbicg 
    "HQ540563"      "114"        "0"       "11" 
    
    [[512]]
          name     length      frame     ncbicg 
    "HQ540564"      "114"        "0"       "11" 
    
    [[513]]
          name     length      frame     ncbicg 
    "HQ540565"      "123"        "0"       "11" 
    
    [[514]]
          name     length      frame     ncbicg 
    "HQ540566"      "114"        "0"       "11" 
    
    [[515]]
          name     length      frame     ncbicg 
    "HQ540567"      "114"        "0"       "11" 
    
    [[516]]
          name     length      frame     ncbicg 
    "HQ540568"      "114"        "0"       "11" 
    
    [[517]]
          name     length      frame     ncbicg 
    "HQ540569"      "114"        "0"       "11" 
    
    [[518]]
          name     length      frame     ncbicg 
    "HQ540570"      "114"        "0"       "11" 
    
    [[519]]
          name     length      frame     ncbicg 
    "HQ540571"      "123"        "0"       "11" 
    
    [[520]]
          name     length      frame     ncbicg 
    "HQ540572"      "114"        "0"       "11" 
    
    [[521]]
          name     length      frame     ncbicg 
    "HQ540573"      "114"        "0"       "11" 
    
    [[522]]
          name     length      frame     ncbicg 
    "HQ540574"      "123"        "0"       "11" 
    
    [[523]]
          name     length      frame     ncbicg 
    "HQ540575"      "123"        "0"       "11" 
    
    [[524]]
          name     length      frame     ncbicg 
    "HQ540576"      "123"        "0"       "11" 
    
    [[525]]
          name     length      frame     ncbicg 
    "HQ540577"      "114"        "0"       "11" 
    
    [[526]]
          name     length      frame     ncbicg 
    "HQ540578"      "114"        "0"       "11" 
    
    [[527]]
          name     length      frame     ncbicg 
    "HQ589040"      "114"        "0"       "11" 
    
    [[528]]
          name     length      frame     ncbicg 
    "HQ589041"      "123"        "0"       "11" 
    
    [[529]]
          name     length      frame     ncbicg 
    "HQ589042"      "123"        "0"       "11" 
    
    [[530]]
          name     length      frame     ncbicg 
    "HQ589043"      "123"        "0"       "11" 
    
    [[531]]
          name     length      frame     ncbicg 
    "HQ589044"      "123"        "0"       "11" 
    
    [[532]]
          name     length      frame     ncbicg 
    "HQ589045"      "123"        "0"       "11" 
    
    [[533]]
          name     length      frame     ncbicg 
    "HQ589046"      "123"        "0"       "11" 
    
    [[534]]
          name     length      frame     ncbicg 
    "HQ589047"      "123"        "0"       "11" 
    
    [[535]]
          name     length      frame     ncbicg 
    "HQ589048"      "123"        "0"       "11" 
    
    [[536]]
          name     length      frame     ncbicg 
    "HQ589049"      "114"        "0"       "11" 
    
    [[537]]
          name     length      frame     ncbicg 
    "HQ692873"      "126"        "0"       "11" 
    
    [[538]]
          name     length      frame     ncbicg 
    "HQ692874"      "123"        "0"       "11" 
    
    [[539]]
          name     length      frame     ncbicg 
    "HQ692876"      "120"        "0"       "11" 
    
    [[540]]
          name     length      frame     ncbicg 
    "HQ692877"      "123"        "0"       "11" 
    
    [[541]]
          name     length      frame     ncbicg 
    "HQ692878"      "123"        "0"       "11" 
    
    [[542]]
          name     length      frame     ncbicg 
    "HQ692879"      "126"        "0"       "11" 
    
    [[543]]
          name     length      frame     ncbicg 
    "HQ692880"      "126"        "0"       "11" 
    
    [[544]]
          name     length      frame     ncbicg 
    "HQ692881"      "123"        "0"       "11" 
    
    [[545]]
          name     length      frame     ncbicg 
    "HQ692882"      "123"        "0"       "11" 
    
    [[546]]
          name     length      frame     ncbicg 
    "HQ692883"      "116"        "1"       "11" 
    
    [[547]]
          name     length      frame     ncbicg 
    "HQ709240"       "75"        "0"       "11" 
    
    [[548]]
          name     length      frame     ncbicg 
    "HQ844243"      "179"        "2"       "11" 
    
    [[549]]
          name     length      frame     ncbicg 
    "HQ844244"      "187"        "0"       "11" 
    
    [[550]]
          name     length      frame     ncbicg 
    "HQ844249"      "185"        "2"       "11" 
    
    [[551]]
          name     length      frame     ncbicg 
    "HQ844250"      "185"        "2"       "11" 
    
    [[552]]
          name     length      frame     ncbicg 
    "HQ844251"      "182"        "2"       "11" 
    
    [[553]]
          name     length      frame     ncbicg 
    "HQ844252"      "185"        "2"       "11" 
    
    [[554]]
          name     length      frame     ncbicg 
    "HQ844253"      "180"        "0"       "11" 
    
    [[555]]
          name     length      frame     ncbicg 
    "HQ997366"      "305"        "0"       "11" 
    
    [[556]]
          name     length      frame     ncbicg 
    "HQ997367"      "305"        "0"       "11" 
    
    [[557]]
          name     length      frame     ncbicg 
    "HQ997368"      "308"        "0"       "11" 
    
    [[558]]
          name     length      frame     ncbicg 
    "JF268583"      "350"        "2"       "11" 
    
    [[559]]
          name     length      frame     ncbicg 
    "JF268584"      "350"        "2"       "11" 
    
    [[560]]
          name     length      frame     ncbicg 
    "JF268585"      "350"        "2"       "11" 
    
    [[561]]
          name     length      frame     ncbicg 
    "JF268586"      "350"        "2"       "11" 
    
    [[562]]
          name     length      frame     ncbicg 
    "JF268587"      "350"        "2"       "11" 
    
    [[563]]
          name     length      frame     ncbicg 
    "JF268588"      "350"        "2"       "11" 
    
    [[564]]
          name     length      frame     ncbicg 
    "JF268589"      "350"        "2"       "11" 
    
    [[565]]
          name     length      frame     ncbicg 
    "JF268590"      "350"        "2"       "11" 
    
    [[566]]
          name     length      frame     ncbicg 
    "JF268591"      "350"        "2"       "11" 
    
    [[567]]
          name     length      frame     ncbicg 
    "JF268592"      "350"        "2"       "11" 
    
    [[568]]
          name     length      frame     ncbicg 
    "JF268593"      "350"        "2"       "11" 
    
    [[569]]
          name     length      frame     ncbicg 
    "JF268594"      "350"        "2"       "11" 
    
    [[570]]
          name     length      frame     ncbicg 
    "JF268595"      "350"        "2"       "11" 
    
    [[571]]
          name     length      frame     ncbicg 
    "JF268596"      "350"        "2"       "11" 
    
    [[572]]
          name     length      frame     ncbicg 
    "JF268597"      "350"        "2"       "11" 
    
    [[573]]
          name     length      frame     ncbicg 
    "JF268598"      "350"        "2"       "11" 
    
    [[574]]
          name     length      frame     ncbicg 
    "JF268599"      "350"        "2"       "11" 
    
    [[575]]
          name     length      frame     ncbicg 
    "JF268600"      "350"        "2"       "11" 
    
    [[576]]
          name     length      frame     ncbicg 
    "JF268601"      "350"        "2"       "11" 
    
    [[577]]
          name     length      frame     ncbicg 
    "JF268602"      "350"        "2"       "11" 
    
    [[578]]
          name     length      frame     ncbicg 
    "JF268603"      "350"        "2"       "11" 
    
    [[579]]
          name     length      frame     ncbicg 
    "JF268604"      "350"        "2"       "11" 
    
    [[580]]
          name     length      frame     ncbicg 
    "JF268605"      "350"        "2"       "11" 
    
    [[581]]
          name     length      frame     ncbicg 
    "JF268606"      "350"        "2"       "11" 
    
    [[582]]
          name     length      frame     ncbicg 
    "JF268607"      "353"        "2"       "11" 
    
    [[583]]
          name     length      frame     ncbicg 
    "JF268608"      "350"        "2"       "11" 
    
    [[584]]
          name     length      frame     ncbicg 
    "JF268609"      "350"        "2"       "11" 
    
    [[585]]
          name     length      frame     ncbicg 
    "JF268610"      "350"        "2"       "11" 
    
    [[586]]
          name     length      frame     ncbicg 
    "JF268611"      "350"        "2"       "11" 
    
    [[587]]
          name     length      frame     ncbicg 
    "JF812083"      "333"        "0"       "11" 
    
    [[588]]
          name     length      frame     ncbicg 
    "JF974259"      "483"        "0"       "11" 
    
    [[589]]
          name     length      frame     ncbicg 
    "JN037845"      "486"        "0"       "11" 
    
    [[590]]
          name     length      frame     ncbicg 
    "JN037846"      "489"        "0"       "11" 
    
    [[591]]
          name     length      frame     ncbicg 
    "JN210554"      "489"        "0"       "11" 
    
    [[592]]
          name     length      frame     ncbicg 
    "JN210555"      "489"        "0"       "11" 
    
    [[593]]
          name     length      frame     ncbicg 
    "JN315350"      "364"        "0"       "11" 
    
    [[594]]
          name     length      frame     ncbicg 
    "JN315351"      "364"        "0"       "11" 
    
    [[595]]
          name     length      frame     ncbicg 
    "JN315352"      "364"        "0"       "11" 
    
    [[596]]
          name     length      frame     ncbicg 
    "JN315353"      "364"        "0"       "11" 
    
    [[597]]
          name     length      frame     ncbicg 
    "JN315354"      "652"        "0"       "11" 
    
    [[598]]
          name     length      frame     ncbicg 
    "JN315355"      "649"        "0"       "11" 
    
    [[599]]
          name     length      frame     ncbicg 
    "JN315356"      "649"        "0"       "11" 
    
    [[600]]
          name     length      frame     ncbicg 
    "JN626460"      "650"        "1"       "11" 
    
    [[601]]
          name     length      frame     ncbicg 
    "JN626461"      "650"        "1"       "11" 
    
    [[602]]
          name     length      frame     ncbicg 
    "JN626462"      "650"        "1"       "11" 
    
    [[603]]
          name     length      frame     ncbicg 
    "JN626463"      "650"        "1"       "11" 
    
    [[604]]
          name     length      frame     ncbicg 
    "JN654513"      "435"        "0"       "11" 
    
    [[605]]
          name     length      frame     ncbicg 
    "JN654514"      "435"        "0"       "11" 
    
    [[606]]
          name     length      frame     ncbicg 
    "JN654515"      "435"        "0"       "11" 
    
    [[607]]
          name     length      frame     ncbicg 
    "JN654516"      "426"        "0"       "11" 
    
    [[608]]
          name     length      frame     ncbicg 
    "JN654517"      "435"        "0"       "11" 
    
    [[609]]
          name     length      frame     ncbicg 
    "JN654518"      "435"        "0"       "11" 
    
    [[610]]
          name     length      frame     ncbicg 
    "JN654519"      "435"        "0"       "11" 
    
    [[611]]
          name     length      frame     ncbicg 
    "JN654520"      "435"        "0"       "11" 
    
    [[612]]
          name     length      frame     ncbicg 
    "JN654521"      "426"        "0"       "11" 
    
    [[613]]
          name     length      frame     ncbicg 
    "JN654522"      "435"        "0"       "11" 
    
    [[614]]
          name     length      frame     ncbicg 
    "JN654523"      "435"        "0"       "11" 
    
    [[615]]
          name     length      frame     ncbicg 
    "JN654524"      "435"        "0"       "11" 
    
    [[616]]
          name     length      frame     ncbicg 
    "JN654525"      "435"        "0"       "11" 
    
    [[617]]
          name     length      frame     ncbicg 
    "JN654526"      "300"        "0"       "11" 
    
    [[618]]
          name     length      frame     ncbicg 
    "JN654527"      "300"        "0"       "11" 
    
    [[619]]
          name     length      frame     ncbicg 
    "JN654528"      "300"        "0"       "11" 
    
    [[620]]
          name     length      frame     ncbicg 
    "JN654529"      "300"        "0"       "11" 
    
    [[621]]
          name     length      frame     ncbicg 
    "JN654530"      "300"        "0"       "11" 
    
    [[622]]
          name     length      frame     ncbicg 
    "JN654531"      "300"        "0"       "11" 
    
    [[623]]
          name     length      frame     ncbicg 
    "JN819066"      "493"        "1"       "11" 
    
    [[624]]
          name     length      frame     ncbicg 
    "JN819067"      "493"        "1"       "11" 
    
    [[625]]
          name     length      frame     ncbicg 
    "JN819068"      "472"        "1"       "11" 
    
    [[626]]
          name     length      frame     ncbicg 
    "JN819069"      "499"        "1"       "11" 
    
    [[627]]
          name     length      frame     ncbicg 
    "JQ314433"      "766"        "0"       "11" 
    
    [[628]]
          name     length      frame     ncbicg 
    "JQ314434"      "762"        "0"       "11" 
    
    [[629]]
          name     length      frame     ncbicg 
    "JQ314435"      "761"        "0"       "11" 
    
    [[630]]
          name     length      frame     ncbicg 
    "JQ314436"      "789"        "0"       "11" 
    
    [[631]]
          name     length      frame     ncbicg 
    "JQ314437"      "756"        "0"       "11" 
    
    [[632]]
          name     length      frame     ncbicg 
    "JQ314438"      "793"        "0"       "11" 
    
    [[633]]
          name     length      frame     ncbicg 
    "JQ314439"      "790"        "0"       "11" 
    
    [[634]]
          name     length      frame     ncbicg 
    "JQ314440"      "792"        "0"       "11" 
    
    [[635]]
          name     length      frame     ncbicg 
    "JQ314441"      "794"        "0"       "11" 
    
    [[636]]
          name     length      frame     ncbicg 
    "JQ314442"      "793"        "0"       "11" 
    
    [[637]]
          name     length      frame     ncbicg 
    "JQ314443"      "775"        "0"       "11" 
    
    [[638]]
          name     length      frame     ncbicg 
    "JQ314444"      "500"        "0"       "11" 
    
    [[639]]
          name     length      frame     ncbicg 
    "JQ414012"     "1606"        "2"       "11" 
    
    [[640]]
          name     length      frame     ncbicg 
    "JQ414013"     "1606"        "2"       "11" 
    
    [[641]]
          name     length      frame     ncbicg 
    "JQ414014"     "1606"        "2"       "11" 
    
    [[642]]
          name     length      frame     ncbicg 
    "JQ414015"     "1606"        "2"       "11" 
    
    [[643]]
          name     length      frame     ncbicg 
    "JQ414016"     "1606"        "2"       "11" 
    
    [[644]]
          name     length      frame     ncbicg 
    "JQ414017"     "1606"        "2"       "11" 
    
    [[645]]
          name     length      frame     ncbicg 
    "JQ414018"     "1606"        "2"       "11" 
    
    [[646]]
          name     length      frame     ncbicg 
    "JQ414019"     "1606"        "2"       "11" 
    
    [[647]]
          name     length      frame     ncbicg 
    "JQ425738"      "281"        "2"       "11" 
    
    [[648]]
          name     length      frame     ncbicg 
    "JX294399"      "670"        "1"       "11" 
    
    [[649]]
          name     length      frame     ncbicg 
    "JX303307"     "3537"        "0"       "11" 
    
    [[650]]
          name     length      frame     ncbicg 
    "JX303308"     "3537"        "0"       "11" 
    
    [[651]]
          name     length      frame     ncbicg 
    "JX303309"     "3537"        "0"       "11" 
    
    [[652]]
          name     length      frame     ncbicg 
    "JX303310"     "3537"        "0"       "11" 
    
    [[653]]
          name     length      frame     ncbicg 
    "JX303311"     "3537"        "0"       "11" 
    
    [[654]]
          name     length      frame     ncbicg 
    "JX303312"     "3537"        "0"       "11" 
    
    [[655]]
          name     length      frame     ncbicg 
    "JX303313"     "3537"        "0"       "11" 
    
    [[656]]
          name     length      frame     ncbicg 
    "JX303314"     "3537"        "0"       "11" 
    
    [[657]]
          name     length      frame     ncbicg 
    "JX303315"     "3537"        "0"       "11" 
    
    [[658]]
          name     length      frame     ncbicg 
    "JX303316"     "3537"        "0"       "11" 
    
    [[659]]
          name     length      frame     ncbicg 
    "JX303317"     "3537"        "0"       "11" 
    
    [[660]]
          name     length      frame     ncbicg 
    "JX303318"     "3537"        "0"       "11" 
    
    [[661]]
          name     length      frame     ncbicg 
    "JX303319"     "3537"        "0"       "11" 
    
    [[662]]
          name     length      frame     ncbicg 
    "JX303320"     "3537"        "0"       "11" 
    
    [[663]]
          name     length      frame     ncbicg 
    "JX303321"     "3537"        "0"       "11" 
    
    [[664]]
          name     length      frame     ncbicg 
    "JX303322"     "3537"        "0"       "11" 
    
    [[665]]
          name     length      frame     ncbicg 
    "JX303323"     "3537"        "0"       "11" 
    
    [[666]]
          name     length      frame     ncbicg 
    "JX303324"     "3537"        "0"       "11" 
    
    [[667]]
          name     length      frame     ncbicg 
    "JX303325"     "3537"        "0"       "11" 
    
    [[668]]
          name     length      frame     ncbicg 
    "JX303326"     "3537"        "0"       "11" 
    
    [[669]]
          name     length      frame     ncbicg 
    "JX303327"     "3537"        "0"       "11" 
    
    [[670]]
          name     length      frame     ncbicg 
    "JX303328"     "3537"        "0"       "11" 
    
    [[671]]
          name     length      frame     ncbicg 
    "JX303329"     "3537"        "0"       "11" 
    
    [[672]]
          name     length      frame     ncbicg 
    "JX303330"     "3537"        "0"       "11" 
    
    [[673]]
          name     length      frame     ncbicg 
    "JX303331"     "3537"        "0"       "11" 
    
    [[674]]
          name     length      frame     ncbicg 
    "JX303332"     "3537"        "0"       "11" 
    
    [[675]]
          name     length      frame     ncbicg 
    "KC344736"      "480"        "0"       "11" 
    
    [[676]]
          name     length      frame     ncbicg 
    "KC344737"      "480"        "0"       "11" 
    
    [[677]]
          name     length      frame     ncbicg 
    "KC491372"      "415"        "0"       "11" 
    
    [[678]]
          name     length      frame     ncbicg 
    "KC491373"      "410"        "0"       "11" 
    
    [[679]]
          name     length      frame     ncbicg 
    "KC491374"      "426"        "0"       "11" 
    
    [[680]]
          name     length      frame     ncbicg 
    "KC692347"     "3537"        "0"       "11" 
    
    [[681]]
          name     length      frame     ncbicg 
    "KC692348"     "3537"        "0"       "11" 
    
    [[682]]
          name     length      frame     ncbicg 
    "KC692349"     "3537"        "0"       "11" 
    
    [[683]]
               name          length           frame          ncbicg 
    "KF444078.RPOB"           "216"             "0"            "11" 
    
    [[684]]
               name          length           frame          ncbicg 
    "KF444079.RPOB"           "216"             "0"            "11" 
    
    [[685]]
               name          length           frame          ncbicg 
    "KF444080.RPOB"           "216"             "0"            "11" 
    
    [[686]]
               name          length           frame          ncbicg 
    "KF444081.RPOB"           "216"             "0"            "11" 
    
    [[687]]
               name          length           frame          ncbicg 
    "KF444082.RPOB"           "216"             "0"            "11" 
    
    [[688]]
               name          length           frame          ncbicg 
    "KF444083.RPOB"           "216"             "0"            "11" 
    
    [[689]]
          name     length      frame     ncbicg 
    "KF751683"      "246"        "0"       "11" 
    
    [[690]]
          name     length      frame     ncbicg 
    "KF751684"      "246"        "0"       "11" 
    
    [[691]]
          name     length      frame     ncbicg 
    "KF751685"      "246"        "0"       "11" 
    
    [[692]]
          name     length      frame     ncbicg 
    "KF751686"      "246"        "0"       "11" 
    
    [[693]]
          name     length      frame     ncbicg 
    "KF751687"      "246"        "0"       "11" 
    
    [[694]]
          name     length      frame     ncbicg 
    "KF751689"      "243"        "0"       "11" 
    
    [[695]]
          name     length      frame     ncbicg 
    "KF751690"      "246"        "0"       "11" 
    
    [[696]]
          name     length      frame     ncbicg 
    "KF751691"      "246"        "0"       "11" 
    
    [[697]]
          name     length      frame     ncbicg 
    "KF751692"      "246"        "0"       "11" 
    
    [[698]]
          name     length      frame     ncbicg 
    "KF751693"      "243"        "0"       "11" 
    
    [[699]]
          name     length      frame     ncbicg 
    "KF877732"      "312"        "0"       "11" 
    
    [[700]]
          name     length      frame     ncbicg 
    "KJ095683"     "3519"        "0"       "11" 
    
    [[701]]
          name     length      frame     ncbicg 
    "KJ410758"      "543"        "2"       "11" 
    
    [[702]]
          name     length      frame     ncbicg 
    "KJ410759"      "543"        "2"       "11" 
    
    [[703]]
          name     length      frame     ncbicg 
    "KJ410760"      "543"        "2"       "11" 
    
    [[704]]
          name     length      frame     ncbicg 
    "KJ410761"      "543"        "2"       "11" 
    
    [[705]]
          name     length      frame     ncbicg 
    "KJ499913"      "688"        "1"       "11" 
    
    [[706]]
          name     length      frame     ncbicg 
    "KJ659894"      "358"        "1"       "11" 
    
    [[707]]
          name     length      frame     ncbicg 
    "KJ659895"      "358"        "1"       "11" 
    
    [[708]]
          name     length      frame     ncbicg 
    "KJ659896"      "358"        "1"       "11" 
    
    [[709]]
          name     length      frame     ncbicg 
    "KJ659897"      "358"        "1"       "11" 
    
    [[710]]
          name     length      frame     ncbicg 
    "KJ659898"      "358"        "1"       "11" 
    
    [[711]]
          name     length      frame     ncbicg 
    "KJ659899"      "358"        "1"       "11" 
    
    [[712]]
          name     length      frame     ncbicg 
    "KJ659900"      "358"        "1"       "11" 
    
    [[713]]
          name     length      frame     ncbicg 
    "KJ659901"      "358"        "1"       "11" 
    
    [[714]]
          name     length      frame     ncbicg 
    "KJ659902"      "358"        "1"       "11" 
    
    [[715]]
          name     length      frame     ncbicg 
    "KJ683741"      "688"        "1"       "11" 
    
    [[716]]
          name     length      frame     ncbicg 
    "KJ683742"      "688"        "1"       "11" 
    
    [[717]]
          name     length      frame     ncbicg 
    "KJ683743"      "688"        "1"       "11" 
    
    [[718]]
          name     length      frame     ncbicg 
    "KJ683744"      "688"        "1"       "11" 
    
    [[719]]
          name     length      frame     ncbicg 
    "KJ683745"      "688"        "1"       "11" 
    
    [[720]]
          name     length      frame     ncbicg 
    "KJ683746"      "688"        "1"       "11" 
    
    [[721]]
          name     length      frame     ncbicg 
    "KJ734728"      "472"        "2"       "11" 
    
    [[722]]
          name     length      frame     ncbicg 
    "KJ734729"      "472"        "2"       "11" 
    
    [[723]]
          name     length      frame     ncbicg 
    "KJ734730"      "472"        "2"       "11" 
    
    [[724]]
          name     length      frame     ncbicg 
    "KJ734731"      "472"        "2"       "11" 
    
    [[725]]
          name     length      frame     ncbicg 
    "KJ734732"      "472"        "2"       "11" 
    
    [[726]]
          name     length      frame     ncbicg 
    "KJ734733"      "472"        "2"       "11" 
    
    [[727]]
          name     length      frame     ncbicg 
    "KJ734734"      "472"        "2"       "11" 
    
    [[728]]
          name     length      frame     ncbicg 
    "KJ734735"      "472"        "2"       "11" 
    
    [[729]]
          name     length      frame     ncbicg 
    "KJ734736"      "472"        "2"       "11" 
    
    [[730]]
          name     length      frame     ncbicg 
    "KJ734737"      "472"        "2"       "11" 
    
    [[731]]
          name     length      frame     ncbicg 
    "KJ734738"      "472"        "2"       "11" 
    
    [[732]]
          name     length      frame     ncbicg 
    "KM234063"      "665"        "2"       "11" 
    
    [[733]]
          name     length      frame     ncbicg 
    "KP658669"      "353"        "2"       "11" 
    
    [[734]]
          name     length      frame     ncbicg 
    "KP658670"      "353"        "2"       "11" 
    
    [[735]]
          name     length      frame     ncbicg 
    "KP658671"      "353"        "2"       "11" 
    
    [[736]]
          name     length      frame     ncbicg 
    "KP658672"      "353"        "2"       "11" 
    
    [[737]]
          name     length      frame     ncbicg 
    "KP658673"      "353"        "2"       "11" 
    
    [[738]]
          name     length      frame     ncbicg 
    "KP658674"      "388"        "2"       "11" 
    
    [[739]]
          name     length      frame     ncbicg 
    "KP658675"      "388"        "2"       "11" 
    
    [[740]]
          name     length      frame     ncbicg 
    "KP658676"      "388"        "2"       "11" 
    
    [[741]]
          name     length      frame     ncbicg 
    "KP658677"      "376"        "2"       "11" 
    
    [[742]]
          name     length      frame     ncbicg 
    "KP658678"      "376"        "2"       "11" 
    
    [[743]]
          name     length      frame     ncbicg 
    "KP658679"      "376"        "2"       "11" 
    
    [[744]]
          name     length      frame     ncbicg 
    "KP658680"      "376"        "2"       "11" 
    
    [[745]]
          name     length      frame     ncbicg 
    "KP658681"      "376"        "2"       "11" 
    
    [[746]]
          name     length      frame     ncbicg 
    "KP658682"      "376"        "2"       "11" 
    
    [[747]]
          name     length      frame     ncbicg 
    "KP658683"      "375"        "0"       "11" 
    
    [[748]]
          name     length      frame     ncbicg 
    "KP658684"      "375"        "0"       "11" 
    
    [[749]]
          name     length      frame     ncbicg 
    "KP658685"      "375"        "0"       "11" 
    
    [[750]]
          name     length      frame     ncbicg 
    "KP658686"      "375"        "0"       "11" 
    
    [[751]]
          name     length      frame     ncbicg 
    "KP658687"      "375"        "0"       "11" 
    
    [[752]]
          name     length      frame     ncbicg 
    "KP658688"      "375"        "0"       "11" 
    
    [[753]]
          name     length      frame     ncbicg 
    "KP658689"      "378"        "0"       "11" 
    
    [[754]]
          name     length      frame     ncbicg 
    "KP658690"      "378"        "0"       "11" 
    
    [[755]]
          name     length      frame     ncbicg 
    "KP658691"      "378"        "0"       "11" 
    
    [[756]]
          name     length      frame     ncbicg 
    "KP658692"      "378"        "0"       "11" 
    
    [[757]]
          name     length      frame     ncbicg 
    "KP658693"      "364"        "2"       "11" 
    
    [[758]]
          name     length      frame     ncbicg 
    "KP658694"      "364"        "2"       "11" 
    
    [[759]]
          name     length      frame     ncbicg 
    "KP658695"      "364"        "2"       "11" 
    
    [[760]]
          name     length      frame     ncbicg 
    "KP658696"      "364"        "2"       "11" 
    
    [[761]]
          name     length      frame     ncbicg 
    "KP658697"      "384"        "0"       "11" 
    
    [[762]]
          name     length      frame     ncbicg 
    "KP658698"      "384"        "0"       "11" 
    
    [[763]]
          name     length      frame     ncbicg 
    "KP658699"      "384"        "0"       "11" 
    
    [[764]]
          name     length      frame     ncbicg 
    "KP658700"      "384"        "0"       "11" 
    
    [[765]]
          name     length      frame     ncbicg 
    "KP658701"      "368"        "0"       "11" 
    
    [[766]]
          name     length      frame     ncbicg 
    "KP658702"      "368"        "0"       "11" 
    
    [[767]]
          name     length      frame     ncbicg 
    "KP658703"      "368"        "0"       "11" 
    
    [[768]]
          name     length      frame     ncbicg 
    "KP658704"      "368"        "0"       "11" 
    
    [[769]]
          name     length      frame     ncbicg 
    "KP658705"      "368"        "0"       "11" 
    
    [[770]]
          name     length      frame     ncbicg 
    "KP658706"      "374"        "0"       "11" 
    
    [[771]]
          name     length      frame     ncbicg 
    "KP658707"      "374"        "0"       "11" 
    
    [[772]]
          name     length      frame     ncbicg 
    "KP658708"      "374"        "0"       "11" 
    
    [[773]]
          name     length      frame     ncbicg 
    "KP658709"      "374"        "0"       "11" 
    
    [[774]]
          name     length      frame     ncbicg 
    "KP658710"      "385"        "0"       "11" 
    
    [[775]]
          name     length      frame     ncbicg 
    "KP658711"      "375"        "0"       "11" 
    
    [[776]]
          name     length      frame     ncbicg 
    "KP658712"      "375"        "0"       "11" 
    
    [[777]]
          name     length      frame     ncbicg 
    "KP658713"      "375"        "0"       "11" 
    
    [[778]]
          name     length      frame     ncbicg 
    "KP658714"      "375"        "0"       "11" 
    
    [[779]]
          name     length      frame     ncbicg 
    "KP658715"      "375"        "0"       "11" 
    
    [[780]]
          name     length      frame     ncbicg 
    "KP658716"      "386"        "0"       "11" 
    
    [[781]]
          name     length      frame     ncbicg 
    "KP658717"      "387"        "0"       "11" 
    
    [[782]]
          name     length      frame     ncbicg 
    "KP658718"      "364"        "2"       "11" 
    
    [[783]]
          name     length      frame     ncbicg 
    "KP658719"      "364"        "2"       "11" 
    
    [[784]]
          name     length      frame     ncbicg 
    "KP658720"      "410"        "2"       "11" 
    
    [[785]]
          name     length      frame     ncbicg 
    "KP732539"      "198"        "0"       "11" 
    
    [[786]]
          name     length      frame     ncbicg 
    "KP744369"     "3519"        "0"       "11" 
    
    [[787]]
          name     length      frame     ncbicg 
    "KP744370"     "3519"        "0"       "11" 
    
    [[788]]
          name     length      frame     ncbicg 
    "KP744371"     "3519"        "0"       "11" 
    
    [[789]]
          name     length      frame     ncbicg 
    "KP744372"     "3519"        "0"       "11" 
    
    [[790]]
          name     length      frame     ncbicg 
    "KP744373"     "3519"        "0"       "11" 
    
    [[791]]
          name     length      frame     ncbicg 
    "KR424777"      "353"        "2"       "11" 
    
    [[792]]
          name     length      frame     ncbicg 
    "KT067737"      "211"        "0"       "11" 
    
    [[793]]
          name     length      frame     ncbicg 
    "KT067738"      "211"        "0"       "11" 
    
    [[794]]
          name     length      frame     ncbicg 
    "KX501218"      "267"        "0"       "11" 
    
    [[795]]
          name     length      frame     ncbicg 
    "KX714230"      "258"        "0"       "11" 
    
    [[796]]
          name     length      frame     ncbicg 
    "KX714231"      "258"        "0"       "11" 
    
    [[797]]
               name          length           frame          ncbicg 
    "KY087999.RPOB"          "3380"             "0"            "11" 
    
    [[798]]
          name     length      frame     ncbicg 
    "KY216023"      "267"        "1"       "11" 
    
    [[799]]
          name     length      frame     ncbicg 
    "KY216024"      "267"        "1"       "11" 
    
    [[800]]
          name     length      frame     ncbicg 
    "KY216025"      "267"        "1"       "11" 
    
    [[801]]
          name     length      frame     ncbicg 
    "KY216026"      "267"        "1"       "11" 
    
    [[802]]
          name     length      frame     ncbicg 
    "KY216027"      "267"        "1"       "11" 
    
    [[803]]
          name     length      frame     ncbicg 
    "KY216028"      "267"        "1"       "11" 
    
    [[804]]
          name     length      frame     ncbicg 
    "KY216029"      "267"        "1"       "11" 
    
    [[805]]
          name     length      frame     ncbicg 
    "KY216030"      "267"        "1"       "11" 
    
    [[806]]
          name     length      frame     ncbicg 
    "KY216031"      "267"        "1"       "11" 
    
    [[807]]
          name     length      frame     ncbicg 
    "KY216032"      "267"        "1"       "11" 
    
    [[808]]
          name     length      frame     ncbicg 
    "KY216033"      "267"        "1"       "11" 
    
    [[809]]
          name     length      frame     ncbicg 
    "KY216034"      "267"        "1"       "11" 
    
    [[810]]
          name     length      frame     ncbicg 
    "KY216035"      "267"        "1"       "11" 
    
    [[811]]
          name     length      frame     ncbicg 
    "KY216036"      "267"        "1"       "11" 
    
    [[812]]
          name     length      frame     ncbicg 
    "KY216037"      "267"        "1"       "11" 
    
    [[813]]
          name     length      frame     ncbicg 
    "KY216038"      "267"        "1"       "11" 
    
    [[814]]
          name     length      frame     ncbicg 
    "KY498650"      "235"        "2"       "11" 
    
    [[815]]
          name     length      frame     ncbicg 
    "KY498651"      "235"        "2"       "11" 
    
    [[816]]
          name     length      frame     ncbicg 
    "KY498652"      "235"        "2"       "11" 
    
    [[817]]
          name     length      frame     ncbicg 
    "KY498653"      "235"        "2"       "11" 
    
    [[818]]
          name     length      frame     ncbicg 
    "KY498654"      "235"        "2"       "11" 
    
    [[819]]
          name     length      frame     ncbicg 
    "KY498655"      "235"        "2"       "11" 
    
    [[820]]
          name     length      frame     ncbicg 
    "KY751931"      "215"        "0"       "11" 
    
    [[821]]
          name     length      frame     ncbicg 
    "KY751932"      "215"        "0"       "11" 
    
    [[822]]
          name     length      frame     ncbicg 
    "KY751933"      "215"        "0"       "11" 
    
    [[823]]
          name     length      frame     ncbicg 
    "KY751934"      "215"        "0"       "11" 
    
    [[824]]
          name     length      frame     ncbicg 
    "KY751935"      "215"        "0"       "11" 
    
    [[825]]
          name     length      frame     ncbicg 
    "KY751936"      "215"        "0"       "11" 
    
    [[826]]
          name     length      frame     ncbicg 
    "KY751937"      "215"        "0"       "11" 
    
    [[827]]
          name     length      frame     ncbicg 
    "KY751938"      "215"        "0"       "11" 
    
    [[828]]
          name     length      frame     ncbicg 
    "KY751939"      "215"        "0"       "11" 
    
    [[829]]
          name     length      frame     ncbicg 
    "KY751940"      "215"        "0"       "11" 
    
    [[830]]
          name     length      frame     ncbicg 
    "KY751941"      "215"        "0"       "11" 
    
    [[831]]
          name     length      frame     ncbicg 
    "KY751942"      "215"        "0"       "11" 
    
    [[832]]
          name     length      frame     ncbicg 
    "KY751943"      "215"        "0"       "11" 
    
    [[833]]
          name     length      frame     ncbicg 
    "KY751944"      "215"        "0"       "11" 
    
    [[834]]
          name     length      frame     ncbicg 
    "KY751945"      "215"        "0"       "11" 
    
    [[835]]
          name     length      frame     ncbicg 
    "KY751946"      "215"        "0"       "11" 
    
    [[836]]
          name     length      frame     ncbicg 
    "KY751947"      "215"        "0"       "11" 
    
    [[837]]
          name     length      frame     ncbicg 
    "KY751948"      "215"        "0"       "11" 
    
    [[838]]
          name     length      frame     ncbicg 
    "KY751949"      "215"        "0"       "11" 
    
    [[839]]
          name     length      frame     ncbicg 
    "KY751950"      "215"        "0"       "11" 
    
    [[840]]
          name     length      frame     ncbicg 
    "KY751951"      "215"        "0"       "11" 
    
    [[841]]
          name     length      frame     ncbicg 
    "KY751952"      "215"        "0"       "11" 
    
    [[842]]
          name     length      frame     ncbicg 
    "LN651304"      "144"        "0"       "11" 
    
    [[843]]
          name     length      frame     ncbicg 
    "LN651305"      "144"        "0"       "11" 
    
    [[844]]
          name     length      frame     ncbicg 
    "LN651306"      "225"        "0"       "11" 
    
    [[845]]
          name     length      frame     ncbicg 
    "LN651307"      "225"        "0"       "11" 
    
    [[846]]
          name     length      frame     ncbicg 
    "LN651308"      "393"        "0"       "11" 
    
    [[847]]
          name     length      frame     ncbicg 
    "LN651309"      "393"        "0"       "11" 
    
    [[848]]
          name     length      frame     ncbicg 
    "LT970844"      "316"        "0"       "11" 
    
    [[849]]
               name          length           frame          ncbicg 
    "MF145294.RPOB"          "1644"             "0"            "11" 
    
    [[850]]
               name          length           frame          ncbicg 
    "MF145295.RPOB"          "1644"             "0"            "11" 
    
    [[851]]
               name          length           frame          ncbicg 
    "MF145296.RPOB"          "1644"             "0"            "11" 
    
    [[852]]
               name          length           frame          ncbicg 
    "MF145297.RPOB"          "1644"             "0"            "11" 
    
    [[853]]
               name          length           frame          ncbicg 
    "MF145298.RPOB"          "1644"             "0"            "11" 
    
    [[854]]
               name          length           frame          ncbicg 
    "MF145299.RPOB"          "1644"             "0"            "11" 
    
    [[855]]
               name          length           frame          ncbicg 
    "MF145300.RPOB"          "1644"             "0"            "11" 
    
    [[856]]
               name          length           frame          ncbicg 
    "MF145301.RPOB"          "1644"             "0"            "11" 
    
    [[857]]
               name          length           frame          ncbicg 
    "MF145302.RPOB"          "1644"             "0"            "11" 
    
    [[858]]
          name     length      frame     ncbicg 
    "MF981086"      "262"        "2"       "11" 
    
    [[859]]
          name     length      frame     ncbicg 
    "MF981087"      "258"        "2"       "11" 
    
    [[860]]
          name     length      frame     ncbicg 
    "MG012797"      "263"        "2"       "11" 
    
    [[861]]
          name     length      frame     ncbicg 
    "MG012798"      "279"        "2"       "11" 
    
    [[862]]
          name     length      frame     ncbicg 
    "MG189519"      "340"        "0"       "11" 
    
    [[863]]
          name     length      frame     ncbicg 
    "MG189521"      "313"        "0"       "11" 
    
    [[864]]
          name     length      frame     ncbicg 
    "MG490375"      "918"        "0"       "11" 
    
    [[865]]
          name     length      frame     ncbicg 
    "MG995041"     "1577"        "0"       "11" 
    
    [[866]]
          name     length      frame     ncbicg 
    "MG995042"     "1577"        "0"       "11" 
    
    [[867]]
          name     length      frame     ncbicg 
    "MG995043"     "1577"        "0"       "11" 
    
    [[868]]
          name     length      frame     ncbicg 
    "MG995044"     "1577"        "0"       "11" 
    
    [[869]]
          name     length      frame     ncbicg 
    "MG995045"     "1577"        "0"       "11" 
    
    [[870]]
          name     length      frame     ncbicg 
    "MG995046"     "1577"        "0"       "11" 
    
    [[871]]
          name     length      frame     ncbicg 
    "MG995047"     "1577"        "0"       "11" 
    
    [[872]]
          name     length      frame     ncbicg 
    "MG995048"     "1577"        "0"       "11" 
    
    [[873]]
          name     length      frame     ncbicg 
    "MG995049"     "1577"        "0"       "11" 
    
    [[874]]
          name     length      frame     ncbicg 
    "MG995050"     "1577"        "0"       "11" 
    
    [[875]]
          name     length      frame     ncbicg 
    "MG995051"     "1577"        "0"       "11" 
    
    [[876]]
          name     length      frame     ncbicg 
    "MG995052"     "1577"        "0"       "11" 
    
    [[877]]
          name     length      frame     ncbicg 
    "MG995053"     "1577"        "0"       "11" 
    
    [[878]]
          name     length      frame     ncbicg 
    "MG995054"     "1577"        "0"       "11" 
    
    [[879]]
          name     length      frame     ncbicg 
    "MG995055"     "1577"        "0"       "11" 
    
    [[880]]
          name     length      frame     ncbicg 
    "MG995056"     "1577"        "0"       "11" 
    
    [[881]]
          name     length      frame     ncbicg 
    "MG995057"     "1577"        "0"       "11" 
    
    [[882]]
          name     length      frame     ncbicg 
    "MG995058"     "1577"        "0"       "11" 
    
    [[883]]
          name     length      frame     ncbicg 
    "MG995059"     "1577"        "0"       "11" 
    
    [[884]]
          name     length      frame     ncbicg 
    "MG995060"     "1577"        "0"       "11" 
    
    [[885]]
          name     length      frame     ncbicg 
    "MG995061"     "1577"        "0"       "11" 
    
    [[886]]
          name     length      frame     ncbicg 
    "MG995062"     "1577"        "0"       "11" 
    
    [[887]]
          name     length      frame     ncbicg 
    "MG995063"     "1577"        "0"       "11" 
    
    [[888]]
          name     length      frame     ncbicg 
    "MG995064"     "1577"        "0"       "11" 
    
    [[889]]
          name     length      frame     ncbicg 
    "MG995065"     "1577"        "0"       "11" 
    
    [[890]]
          name     length      frame     ncbicg 
    "MG995066"     "1577"        "0"       "11" 
    
    [[891]]
          name     length      frame     ncbicg 
    "MG995067"     "1577"        "0"       "11" 
    
    [[892]]
          name     length      frame     ncbicg 
    "MG995068"     "1577"        "0"       "11" 
    
    [[893]]
          name     length      frame     ncbicg 
    "MG995069"     "1577"        "0"       "11" 
    
    [[894]]
          name     length      frame     ncbicg 
    "MG995070"     "1577"        "0"       "11" 
    
    [[895]]
          name     length      frame     ncbicg 
    "MG995071"     "1577"        "0"       "11" 
    
    [[896]]
          name     length      frame     ncbicg 
    "MG995072"     "1577"        "0"       "11" 
    
    [[897]]
          name     length      frame     ncbicg 
    "MG995073"     "1577"        "0"       "11" 
    
    [[898]]
          name     length      frame     ncbicg 
    "MG995074"     "1577"        "0"       "11" 
    
    [[899]]
          name     length      frame     ncbicg 
    "MG995075"     "1577"        "0"       "11" 
    
    [[900]]
          name     length      frame     ncbicg 
    "MG995076"     "1577"        "0"       "11" 
    
    [[901]]
          name     length      frame     ncbicg 
    "MG995077"     "1577"        "0"       "11" 
    
    [[902]]
          name     length      frame     ncbicg 
    "MG995078"     "1577"        "0"       "11" 
    
    [[903]]
          name     length      frame     ncbicg 
    "MG995079"     "1577"        "0"       "11" 
    
    [[904]]
          name     length      frame     ncbicg 
    "MG995080"     "1577"        "0"       "11" 
    
    [[905]]
          name     length      frame     ncbicg 
    "MG995081"     "1577"        "0"       "11" 
    
    [[906]]
          name     length      frame     ncbicg 
    "MG995082"     "1577"        "0"       "11" 
    
    [[907]]
          name     length      frame     ncbicg 
    "MG995083"     "1577"        "0"       "11" 
    
    [[908]]
          name     length      frame     ncbicg 
    "MG995084"     "1577"        "0"       "11" 
    
    [[909]]
          name     length      frame     ncbicg 
    "MG995085"     "1577"        "0"       "11" 
    
    [[910]]
          name     length      frame     ncbicg 
    "MG995086"     "1577"        "0"       "11" 
    
    [[911]]
          name     length      frame     ncbicg 
    "MG995087"     "1577"        "0"       "11" 
    
    [[912]]
          name     length      frame     ncbicg 
    "MG995088"     "1577"        "0"       "11" 
    
    [[913]]
          name     length      frame     ncbicg 
    "MG995089"     "1577"        "0"       "11" 
    
    [[914]]
          name     length      frame     ncbicg 
    "MG995090"     "1577"        "0"       "11" 
    
    [[915]]
          name     length      frame     ncbicg 
    "MG995091"     "1577"        "0"       "11" 
    
    [[916]]
          name     length      frame     ncbicg 
    "MG995092"     "1577"        "0"       "11" 
    
    [[917]]
          name     length      frame     ncbicg 
    "MG995093"     "1577"        "0"       "11" 
    
    [[918]]
          name     length      frame     ncbicg 
    "MG995094"     "1577"        "0"       "11" 
    
    [[919]]
          name     length      frame     ncbicg 
    "MG995095"     "1577"        "0"       "11" 
    
    [[920]]
          name     length      frame     ncbicg 
    "MG995096"     "1577"        "0"       "11" 
    
    [[921]]
          name     length      frame     ncbicg 
    "MG995097"     "1577"        "0"       "11" 
    
    [[922]]
          name     length      frame     ncbicg 
    "MG995098"     "1577"        "0"       "11" 
    
    [[923]]
          name     length      frame     ncbicg 
    "MG995099"     "1577"        "0"       "11" 
    
    [[924]]
          name     length      frame     ncbicg 
    "MG995100"     "1577"        "0"       "11" 
    
    [[925]]
          name     length      frame     ncbicg 
    "MG995101"     "1577"        "0"       "11" 
    
    [[926]]
          name     length      frame     ncbicg 
    "MG995102"     "1577"        "0"       "11" 
    
    [[927]]
          name     length      frame     ncbicg 
    "MG995103"     "1577"        "0"       "11" 
    
    [[928]]
          name     length      frame     ncbicg 
    "MG995104"     "1577"        "0"       "11" 
    
    [[929]]
          name     length      frame     ncbicg 
    "MG995105"     "1577"        "0"       "11" 
    
    [[930]]
          name     length      frame     ncbicg 
    "MG995106"     "1577"        "0"       "11" 
    
    [[931]]
          name     length      frame     ncbicg 
    "MG995107"     "1577"        "0"       "11" 
    
    [[932]]
          name     length      frame     ncbicg 
    "MG995108"     "1577"        "0"       "11" 
    
    [[933]]
          name     length      frame     ncbicg 
    "MG995109"     "1577"        "0"       "11" 
    
    [[934]]
          name     length      frame     ncbicg 
    "MG995110"     "1577"        "0"       "11" 
    
    [[935]]
          name     length      frame     ncbicg 
    "MG995111"     "1577"        "0"       "11" 
    
    [[936]]
          name     length      frame     ncbicg 
    "MG995112"     "1577"        "0"       "11" 
    
    [[937]]
          name     length      frame     ncbicg 
    "MG995113"     "1577"        "0"       "11" 
    
    [[938]]
          name     length      frame     ncbicg 
    "MG995114"     "1577"        "0"       "11" 
    
    [[939]]
          name     length      frame     ncbicg 
    "MG995115"     "1577"        "0"       "11" 
    
    [[940]]
              name         length          frame         ncbicg 
    "MSGRPOB.RPOB"         "3534"            "0"           "11" 
    
    [[941]]
               name          length           frame          ncbicg 
    "MTU12205.RPOB"          "3278"             "0"            "11" 
    
    [[942]]
          name     length      frame     ncbicg 
    "MTU70422"      "174"        "0"       "11" 
    
    [[943]]
        name   length    frame   ncbicg 
    "S71246"     "69"      "0"     "11" 
    



```R
q2$req    # 463 Sequences as on Jan 29th 2018
```


    [[1]]
               name          length           frame          ncbicg 
    "AE005174.RPOB"          "4029"             "0"            "11" 
    
    [[2]]
               name          length           frame          ncbicg 
    "AE014075.RPOB"          "4062"             "0"            "11" 
    
    [[3]]
          name     length      frame     ncbicg 
    "AJ854258"      "715"        "0"       "11" 
    
    [[4]]
               name          length           frame          ncbicg 
    "AM946981.RPOB"          "4029"             "0"            "11" 
    
    [[5]]
               name          length           frame          ncbicg 
    "AP009048.RPOB"          "4029"             "0"            "11" 
    
    [[6]]
               name          length           frame          ncbicg 
    "AP010953.RPOB"          "4029"             "0"            "11" 
    
    [[7]]
               name          length           frame          ncbicg 
    "AP010958.RPOB"          "4029"             "0"            "11" 
    
    [[8]]
               name          length           frame          ncbicg 
    "AP010960.RPOB"          "4029"             "0"            "11" 
    
    [[9]]
               name          length           frame          ncbicg 
    "AP012030.RPOB"          "4029"             "0"            "11" 
    
    [[10]]
               name          length           frame          ncbicg 
    "AP012306.RPOB"          "4029"             "0"            "11" 
    
    [[11]]
               name          length           frame          ncbicg 
    "AP017610.RPOB"          "4029"             "0"            "11" 
    
    [[12]]
               name          length           frame          ncbicg 
    "AP017617.RPOB"          "4029"             "0"            "11" 
    
    [[13]]
               name          length           frame          ncbicg 
    "AP017620.RPOB"          "4029"             "0"            "11" 
    
    [[14]]
          name     length      frame     ncbicg 
    "AY426236"      "330"        "0"       "11" 
    
    [[15]]
               name          length           frame          ncbicg 
    "CP000243.RPOB"          "4062"             "0"            "11" 
    
    [[16]]
               name          length           frame          ncbicg 
    "CP000468.RPOB"          "4062"             "0"            "11" 
    
    [[17]]
               name          length           frame          ncbicg 
    "CP000800.RPOB"          "4029"             "0"            "11" 
    
    [[18]]
               name          length           frame          ncbicg 
    "CP000802.RPOB"          "4029"             "0"            "11" 
    
    [[19]]
               name          length           frame          ncbicg 
    "CP000819.RPOB"          "4029"             "0"            "11" 
    
    [[20]]
               name          length           frame          ncbicg 
    "CP000948.RPOB"          "4029"             "0"            "11" 
    
    [[21]]
               name          length           frame          ncbicg 
    "CP000970.RPOB"          "4029"             "0"            "11" 
    
    [[22]]
               name          length           frame          ncbicg 
    "CP001164.RPOB"          "4029"             "0"            "11" 
    
    [[23]]
               name          length           frame          ncbicg 
    "CP001368.RPOB"          "4029"             "0"            "11" 
    
    [[24]]
               name          length           frame          ncbicg 
    "CP001396.RPOB"          "4029"             "0"            "11" 
    
    [[25]]
               name          length           frame          ncbicg 
    "CP001509.RPOB"          "4029"             "0"            "11" 
    
    [[26]]
               name          length           frame          ncbicg 
    "CP001671.RPOB"          "4029"             "0"            "11" 
    
    [[27]]
               name          length           frame          ncbicg 
    "CP001846.RPOB"          "4029"             "0"            "11" 
    
    [[28]]
               name          length           frame          ncbicg 
    "CP001855.RPOB"          "4029"             "0"            "11" 
    
    [[29]]
               name          length           frame          ncbicg 
    "CP001925.RPOB"          "4029"             "0"            "11" 
    
    [[30]]
               name          length           frame          ncbicg 
    "CP001969.RPOB"          "4029"             "0"            "11" 
    
    [[31]]
               name          length           frame          ncbicg 
    "CP002167.RPOB"          "4029"             "0"            "11" 
    
    [[32]]
               name          length           frame          ncbicg 
    "CP002185.RPOB"          "4029"             "0"            "11" 
    
    [[33]]
               name          length           frame          ncbicg 
    "CP002211.RPOB"          "4062"             "0"            "11" 
    
    [[34]]
               name          length           frame          ncbicg 
    "CP002212.RPOB"          "4062"             "0"            "11" 
    
    [[35]]
               name          length           frame          ncbicg 
    "CP002291.RPOB"          "4029"             "0"            "11" 
    
    [[36]]
               name          length           frame          ncbicg 
    "CP002729.RPOB"          "4029"             "0"            "11" 
    
    [[37]]
               name          length           frame          ncbicg 
    "CP002967.RPOB"          "4029"             "0"            "11" 
    
    [[38]]
               name          length           frame          ncbicg 
    "CP002970.RPOB"          "4029"             "0"            "11" 
    
    [[39]]
               name          length           frame          ncbicg 
    "CP003034.RPOB"          "4029"             "0"            "11" 
    
    [[40]]
               name          length           frame          ncbicg 
    "CP003109.RPOB"          "4029"             "0"            "11" 
    
    [[41]]
               name          length           frame          ncbicg 
    "CP003289.RPOB"          "4029"             "0"            "11" 
    
    [[42]]
               name          length           frame          ncbicg 
    "CP003297.RPOB"          "4029"             "0"            "11" 
    
    [[43]]
               name          length           frame          ncbicg 
    "CP003301.RPOB"          "4029"             "0"            "11" 
    
    [[44]]
               name          length           frame          ncbicg 
    "CP006027.RPOB"          "4029"             "0"            "11" 
    
    [[45]]
               name          length           frame          ncbicg 
    "CP006262.RPOB"          "4029"             "0"            "11" 
    
    [[46]]
               name          length           frame          ncbicg 
    "CP006584.RPOB"          "4029"             "0"            "11" 
    
    [[47]]
               name          length           frame          ncbicg 
    "CP006632.RPOB"          "4029"             "0"            "11" 
    
    [[48]]
               name          length           frame          ncbicg 
    "CP006636.RPOB"          "4029"             "0"            "11" 
    
    [[49]]
               name          length           frame          ncbicg 
    "CP006698.RPOB"          "4029"             "0"            "11" 
    
    [[50]]
               name          length           frame          ncbicg 
    "CP006784.RPOB"          "4029"             "0"            "11" 
    
    [[51]]
               name          length           frame          ncbicg 
    "CP006830.RPOB"          "4029"             "0"            "11" 
    
    [[52]]
               name          length           frame          ncbicg 
    "CP006834.RPOB"          "4029"             "0"            "11" 
    
    [[53]]
               name          length           frame          ncbicg 
    "CP007133.RPOB"          "4029"             "0"            "11" 
    
    [[54]]
               name          length           frame          ncbicg 
    "CP007136.RPOB"          "4029"             "0"            "11" 
    
    [[55]]
               name          length           frame          ncbicg 
    "CP007149.RPOB"          "4029"             "0"            "11" 
    
    [[56]]
               name          length           frame          ncbicg 
    "CP007265.RPOB"          "4029"             "0"            "11" 
    
    [[57]]
               name          length           frame          ncbicg 
    "CP007275.RPOB"          "4029"             "0"            "11" 
    
    [[58]]
               name          length           frame          ncbicg 
    "CP007390.RPOB"          "4029"             "0"            "11" 
    
    [[59]]
               name          length           frame          ncbicg 
    "CP007391.RPOB"          "4029"             "0"            "11" 
    
    [[60]]
               name          length           frame          ncbicg 
    "CP007392.RPOB"          "4029"             "0"            "11" 
    
    [[61]]
               name          length           frame          ncbicg 
    "CP007393.RPOB"          "4029"             "0"            "11" 
    
    [[62]]
               name          length           frame          ncbicg 
    "CP007394.RPOB"          "4029"             "0"            "11" 
    
    [[63]]
               name          length           frame          ncbicg 
    "CP007442.RPOB"          "4029"             "0"            "11" 
    
    [[64]]
               name          length           frame          ncbicg 
    "CP007491.RPOB"          "4062"             "0"            "11" 
    
    [[65]]
               name          length           frame          ncbicg 
    "CP007594.RPOB"          "4029"             "0"            "11" 
    
    [[66]]
               name          length           frame          ncbicg 
    "CP007799.RPOB"          "4029"             "0"            "11" 
    
    [[67]]
               name          length           frame          ncbicg 
    "CP008697.RPOB"          "4029"             "0"            "11" 
    
    [[68]]
               name          length           frame          ncbicg 
    "CP008801.RPOB"          "4029"             "0"            "11" 
    
    [[69]]
               name          length           frame          ncbicg 
    "CP008805.RPOB"          "4029"             "0"            "11" 
    
    [[70]]
               name          length           frame          ncbicg 
    "CP009072.RPOB"          "4029"             "0"            "11" 
    
    [[71]]
               name          length           frame          ncbicg 
    "CP009104.RPOB"          "4029"             "0"            "11" 
    
    [[72]]
               name          length           frame          ncbicg 
    "CP009106.RPOB"          "4029"             "0"            "11" 
    
    [[73]]
               name          length           frame          ncbicg 
    "CP009166.RPOB"          "4029"             "0"            "11" 
    
    [[74]]
               name          length           frame          ncbicg 
    "CP009273.RPOB"          "4029"             "0"            "11" 
    
    [[75]]
               name          length           frame          ncbicg 
    "CP009578.RPOB"          "4029"             "0"            "11" 
    
    [[76]]
               name          length           frame          ncbicg 
    "CP009644.RPOB"          "4029"             "0"            "11" 
    
    [[77]]
               name          length           frame          ncbicg 
    "CP009685.RPOB"          "4029"             "0"            "11" 
    
    [[78]]
               name          length           frame          ncbicg 
    "CP009789.RPOB"          "4029"             "0"            "11" 
    
    [[79]]
               name          length           frame          ncbicg 
    "CP009859.RPOB"          "4029"             "0"            "11" 
    
    [[80]]
               name          length           frame          ncbicg 
    "CP010116.RPOB"          "4029"             "0"            "11" 
    
    [[81]]
               name          length           frame          ncbicg 
    "CP010117.RPOB"          "4029"             "0"            "11" 
    
    [[82]]
               name          length           frame          ncbicg 
    "CP010119.RPOB"          "4029"             "0"            "11" 
    
    [[83]]
               name          length           frame          ncbicg 
    "CP010121.RPOB"          "4029"             "0"            "11" 
    
    [[84]]
               name          length           frame          ncbicg 
    "CP010122.RPOB"          "4029"             "0"            "11" 
    
    [[85]]
               name          length           frame          ncbicg 
    "CP010125.RPOB"          "4029"             "0"            "11" 
    
    [[86]]
               name          length           frame          ncbicg 
    "CP010129.RPOB"          "4029"             "0"            "11" 
    
    [[87]]
               name          length           frame          ncbicg 
    "CP010132.RPOB"          "4029"             "0"            "11" 
    
    [[88]]
               name          length           frame          ncbicg 
    "CP010133.RPOB"          "4029"             "0"            "11" 
    
    [[89]]
               name          length           frame          ncbicg 
    "CP010134.RPOB"          "4029"             "0"            "11" 
    
    [[90]]
               name          length           frame          ncbicg 
    "CP010137.RPOB"          "4029"             "0"            "11" 
    
    [[91]]
               name          length           frame          ncbicg 
    "CP010140.RPOB"          "4029"             "0"            "11" 
    
    [[92]]
               name          length           frame          ncbicg 
    "CP010143.RPOB"          "4029"             "0"            "11" 
    
    [[93]]
               name          length           frame          ncbicg 
    "CP010145.RPOB"          "4029"             "0"            "11" 
    
    [[94]]
               name          length           frame          ncbicg 
    "CP010148.RPOB"          "4029"             "0"            "11" 
    
    [[95]]
               name          length           frame          ncbicg 
    "CP010150.RPOB"          "4029"             "0"            "11" 
    
    [[96]]
               name          length           frame          ncbicg 
    "CP010151.RPOB"          "4029"             "0"            "11" 
    
    [[97]]
               name          length           frame          ncbicg 
    "CP010152.RPOB"          "4029"             "0"            "11" 
    
    [[98]]
               name          length           frame          ncbicg 
    "CP010157.RPOB"          "4029"             "0"            "11" 
    
    [[99]]
               name          length           frame          ncbicg 
    "CP010160.RPOB"          "4029"             "0"            "11" 
    
    [[100]]
               name          length           frame          ncbicg 
    "CP010163.RPOB"          "4029"             "0"            "11" 
    
    [[101]]
               name          length           frame          ncbicg 
    "CP010167.RPOB"          "4029"             "0"            "11" 
    
    [[102]]
               name          length           frame          ncbicg 
    "CP010169.RPOB"          "4029"             "0"            "11" 
    
    [[103]]
               name          length           frame          ncbicg 
    "CP010170.RPOB"          "4029"             "0"            "11" 
    
    [[104]]
               name          length           frame          ncbicg 
    "CP010171.RPOB"          "4029"             "0"            "11" 
    
    [[105]]
               name          length           frame          ncbicg 
    "CP010172.RPOB"          "4029"             "0"            "11" 
    
    [[106]]
               name          length           frame          ncbicg 
    "CP010176.RPOB"          "4029"             "0"            "11" 
    
    [[107]]
               name          length           frame          ncbicg 
    "CP010177.RPOB"          "4029"             "0"            "11" 
    
    [[108]]
               name          length           frame          ncbicg 
    "CP010178.RPOB"          "4029"             "0"            "11" 
    
    [[109]]
               name          length           frame          ncbicg 
    "CP010180.RPOB"          "4029"             "0"            "11" 
    
    [[110]]
               name          length           frame          ncbicg 
    "CP010183.RPOB"          "4029"             "0"            "11" 
    
    [[111]]
               name          length           frame          ncbicg 
    "CP010186.RPOB"          "4029"             "0"            "11" 
    
    [[112]]
               name          length           frame          ncbicg 
    "CP010191.RPOB"          "4029"             "0"            "11" 
    
    [[113]]
               name          length           frame          ncbicg 
    "CP010196.RPOB"          "4029"             "0"            "11" 
    
    [[114]]
               name          length           frame          ncbicg 
    "CP010200.RPOB"          "4029"             "0"            "11" 
    
    [[115]]
               name          length           frame          ncbicg 
    "CP010206.RPOB"          "4029"             "0"            "11" 
    
    [[116]]
               name          length           frame          ncbicg 
    "CP010213.RPOB"          "4029"             "0"            "11" 
    
    [[117]]
               name          length           frame          ncbicg 
    "CP010219.RPOB"          "4029"             "0"            "11" 
    
    [[118]]
               name          length           frame          ncbicg 
    "CP010221.RPOB"          "4029"             "0"            "11" 
    
    [[119]]
               name          length           frame          ncbicg 
    "CP010226.RPOB"          "4029"             "0"            "11" 
    
    [[120]]
               name          length           frame          ncbicg 
    "CP010228.RPOB"          "4029"             "0"            "11" 
    
    [[121]]
               name          length           frame          ncbicg 
    "CP010229.RPOB"          "4029"             "0"            "11" 
    
    [[122]]
               name          length           frame          ncbicg 
    "CP010230.RPOB"          "4029"             "0"            "11" 
    
    [[123]]
               name          length           frame          ncbicg 
    "CP010231.RPOB"          "4029"             "0"            "11" 
    
    [[124]]
               name          length           frame          ncbicg 
    "CP010235.RPOB"          "4029"             "0"            "11" 
    
    [[125]]
               name          length           frame          ncbicg 
    "CP010236.RPOB"          "4029"             "0"            "11" 
    
    [[126]]
               name          length           frame          ncbicg 
    "CP010237.RPOB"          "4029"             "0"            "11" 
    
    [[127]]
               name          length           frame          ncbicg 
    "CP010238.RPOB"          "4029"             "0"            "11" 
    
    [[128]]
               name          length           frame          ncbicg 
    "CP010240.RPOB"          "4029"             "0"            "11" 
    
    [[129]]
               name          length           frame          ncbicg 
    "CP010242.RPOB"          "4029"             "0"            "11" 
    
    [[130]]
               name          length           frame          ncbicg 
    "CP010304.RPOB"          "4029"             "0"            "11" 
    
    [[131]]
               name          length           frame          ncbicg 
    "CP010315.RPOB"          "4029"             "0"            "11" 
    
    [[132]]
               name          length           frame          ncbicg 
    "CP010344.RPOB"          "4029"             "0"            "11" 
    
    [[133]]
               name          length           frame          ncbicg 
    "CP010371.RPOB"          "4029"             "0"            "11" 
    
    [[134]]
               name          length           frame          ncbicg 
    "CP010438.RPOB"          "4029"             "0"            "11" 
    
    [[135]]
               name          length           frame          ncbicg 
    "CP010439.RPOB"          "4029"             "0"            "11" 
    
    [[136]]
               name          length           frame          ncbicg 
    "CP010440.RPOB"          "4029"             "0"            "11" 
    
    [[137]]
               name          length           frame          ncbicg 
    "CP010441.RPOB"          "4029"             "0"            "11" 
    
    [[138]]
               name          length           frame          ncbicg 
    "CP010442.RPOB"          "4029"             "0"            "11" 
    
    [[139]]
               name          length           frame          ncbicg 
    "CP010443.RPOB"          "4029"             "0"            "11" 
    
    [[140]]
               name          length           frame          ncbicg 
    "CP010444.RPOB"          "4029"             "0"            "11" 
    
    [[141]]
               name          length           frame          ncbicg 
    "CP010445.RPOB"          "4029"             "0"            "11" 
    
    [[142]]
               name          length           frame          ncbicg 
    "CP010585.RPOB"          "4029"             "0"            "11" 
    
    [[143]]
               name          length           frame          ncbicg 
    "CP010816.RPOB"          "4029"             "0"            "11" 
    
    [[144]]
               name          length           frame          ncbicg 
    "CP010876.RPOB"          "4029"             "0"            "11" 
    
    [[145]]
               name          length           frame          ncbicg 
    "CP011018.RPOB"          "4029"             "0"            "11" 
    
    [[146]]
               name          length           frame          ncbicg 
    "CP011061.RPOB"          "4029"             "0"            "11" 
    
    [[147]]
               name          length           frame          ncbicg 
    "CP011113.RPOB"          "4029"             "0"            "11" 
    
    [[148]]
               name          length           frame          ncbicg 
    "CP011124.RPOB"          "4029"             "0"            "11" 
    
    [[149]]
               name          length           frame          ncbicg 
    "CP011134.RPOB"          "4029"             "0"            "11" 
    
    [[150]]
               name          length           frame          ncbicg 
    "CP011320.RPOB"          "4029"             "0"            "11" 
    
    [[151]]
               name          length           frame          ncbicg 
    "CP011321.RPOB"          "4029"             "0"            "11" 
    
    [[152]]
               name          length           frame          ncbicg 
    "CP011322.RPOB"          "4029"             "0"            "11" 
    
    [[153]]
               name          length           frame          ncbicg 
    "CP011323.RPOB"          "4029"             "0"            "11" 
    
    [[154]]
               name          length           frame          ncbicg 
    "CP011324.RPOB"          "4029"             "0"            "11" 
    
    [[155]]
               name          length           frame          ncbicg 
    "CP011331.RPOB"          "4029"             "0"            "11" 
    
    [[156]]
               name          length           frame          ncbicg 
    "CP011342.RPOB"          "4029"             "0"            "11" 
    
    [[157]]
               name          length           frame          ncbicg 
    "CP011343.RPOB"          "4029"             "0"            "11" 
    
    [[158]]
               name          length           frame          ncbicg 
    "CP011416.RPOB"          "4029"             "0"            "11" 
    
    [[159]]
               name          length           frame          ncbicg 
    "CP011495.RPOB"          "4029"             "0"            "11" 
    
    [[160]]
               name          length           frame          ncbicg 
    "CP011915.RPOB"          "4029"             "0"            "11" 
    
    [[161]]
               name          length           frame          ncbicg 
    "CP011938.RPOB"          "4029"             "0"            "11" 
    
    [[162]]
               name          length           frame          ncbicg 
    "CP012112.RPOB"          "4029"             "0"            "11" 
    
    [[163]]
               name          length           frame          ncbicg 
    "CP012125.RPOB"          "4029"             "0"            "11" 
    
    [[164]]
               name          length           frame          ncbicg 
    "CP012126.RPOB"          "4029"             "0"            "11" 
    
    [[165]]
               name          length           frame          ncbicg 
    "CP012127.RPOB"          "4029"             "0"            "11" 
    
    [[166]]
               name          length           frame          ncbicg 
    "CP012625.RPOB"          "4029"             "0"            "11" 
    
    [[167]]
               name          length           frame          ncbicg 
    "CP012631.RPOB"          "4029"             "0"            "11" 
    
    [[168]]
               name          length           frame          ncbicg 
    "CP012633.RPOB"          "4029"             "0"            "11" 
    
    [[169]]
               name          length           frame          ncbicg 
    "CP012635.RPOB"          "4029"             "0"            "11" 
    
    [[170]]
               name          length           frame          ncbicg 
    "CP012802.RPOB"          "4029"             "0"            "11" 
    
    [[171]]
               name          length           frame          ncbicg 
    "CP012868.RPOB"          "4029"             "0"            "11" 
    
    [[172]]
               name          length           frame          ncbicg 
    "CP012869.RPOB"          "4029"             "0"            "11" 
    
    [[173]]
               name          length           frame          ncbicg 
    "CP012870.RPOB"          "4029"             "0"            "11" 
    
    [[174]]
               name          length           frame          ncbicg 
    "CP013025.RPOB"          "4029"             "0"            "11" 
    
    [[175]]
               name          length           frame          ncbicg 
    "CP013029.RPOB"          "4029"             "0"            "11" 
    
    [[176]]
               name          length           frame          ncbicg 
    "CP013031.RPOB"          "4029"             "0"            "11" 
    
    [[177]]
               name          length           frame          ncbicg 
    "CP013048.RPOB"          "4029"             "0"            "11" 
    
    [[178]]
               name          length           frame          ncbicg 
    "CP013112.RPOB"          "4029"             "0"            "11" 
    
    [[179]]
               name          length           frame          ncbicg 
    "CP013253.RPOB"          "4029"             "0"            "11" 
    
    [[180]]
               name          length           frame          ncbicg 
    "CP013483.RPOB"          "4029"             "0"            "11" 
    
    [[181]]
               name          length           frame          ncbicg 
    "CP013658.RPOB"          "4029"             "0"            "11" 
    
    [[182]]
               name          length           frame          ncbicg 
    "CP013662.RPOB"          "4029"             "0"            "11" 
    
    [[183]]
               name          length           frame          ncbicg 
    "CP013663.RPOB"          "4029"             "0"            "11" 
    
    [[184]]
               name          length           frame          ncbicg 
    "CP013831.RPOB"          "4029"             "0"            "11" 
    
    [[185]]
               name          length           frame          ncbicg 
    "CP013835.RPOB"          "4029"             "0"            "11" 
    
    [[186]]
               name          length           frame          ncbicg 
    "CP013837.RPOB"          "4029"             "0"            "11" 
    
    [[187]]
               name          length           frame          ncbicg 
    "CP013952.RPOB"          "4029"             "0"            "11" 
    
    [[188]]
               name          length           frame          ncbicg 
    "CP014092.RPOB"          "4029"             "0"            "11" 
    
    [[189]]
               name          length           frame          ncbicg 
    "CP014111.RPOB"          "4029"             "0"            "11" 
    
    [[190]]
               name          length           frame          ncbicg 
    "CP014197.RPOB"          "4029"             "0"            "11" 
    
    [[191]]
               name          length           frame          ncbicg 
    "CP014225.RPOB"          "4029"             "0"            "11" 
    
    [[192]]
               name          length           frame          ncbicg 
    "CP014268.RPOB"          "4029"             "0"            "11" 
    
    [[193]]
               name          length           frame          ncbicg 
    "CP014269.RPOB"          "4029"             "0"            "11" 
    
    [[194]]
               name          length           frame          ncbicg 
    "CP014270.RPOB"          "4029"             "0"            "11" 
    
    [[195]]
               name          length           frame          ncbicg 
    "CP014272.RPOB"          "4029"             "0"            "11" 
    
    [[196]]
               name          length           frame          ncbicg 
    "CP014316.RPOB"          "4029"             "0"            "11" 
    
    [[197]]
               name          length           frame          ncbicg 
    "CP014348.RPOB"          "4029"             "0"            "11" 
    
    [[198]]
               name          length           frame          ncbicg 
    "CP014488.RPOB"          "4029"             "0"            "11" 
    
    [[199]]
               name          length           frame          ncbicg 
    "CP014492.RPOB"          "4029"             "0"            "11" 
    
    [[200]]
               name          length           frame          ncbicg 
    "CP014495.RPOB"          "4029"             "0"            "11" 
    
    [[201]]
               name          length           frame          ncbicg 
    "CP014497.RPOB"          "4029"             "0"            "11" 
    
    [[202]]
               name          length           frame          ncbicg 
    "CP014522.RPOB"          "4029"             "0"            "11" 
    
    [[203]]
               name          length           frame          ncbicg 
    "CP014641.RPOB"          "4029"             "0"            "11" 
    
    [[204]]
               name          length           frame          ncbicg 
    "CP014642.RPOB"          "4029"             "0"            "11" 
    
    [[205]]
               name          length           frame          ncbicg 
    "CP014667.RPOB"          "4029"             "0"            "11" 
    
    [[206]]
               name          length           frame          ncbicg 
    "CP014752.RPOB"          "4029"             "0"            "11" 
    
    [[207]]
               name          length           frame          ncbicg 
    "CP015020.RPOB"          "4029"             "0"            "11" 
    
    [[208]]
               name          length           frame          ncbicg 
    "CP015023.RPOB"          "4029"             "0"            "11" 
    
    [[209]]
               name          length           frame          ncbicg 
    "CP015069.RPOB"          "4029"             "0"            "11" 
    
    [[210]]
               name          length           frame          ncbicg 
    "CP015074.RPOB"          "4029"             "0"            "11" 
    
    [[211]]
               name          length           frame          ncbicg 
    "CP015076.RPOB"          "4029"             "0"            "11" 
    
    [[212]]
                 name            length             frame            ncbicg 
    "CP015085.PE3285"            "4062"               "0"              "11" 
    
    [[213]]
               name          length           frame          ncbicg 
    "CP015138.RPOB"          "4029"             "0"            "11" 
    
    [[214]]
               name          length           frame          ncbicg 
    "CP015159.RPOB"          "4029"             "0"            "11" 
    
    [[215]]
               name          length           frame          ncbicg 
    "CP015228.RPOB"          "4029"             "0"            "11" 
    
    [[216]]
               name          length           frame          ncbicg 
    "CP015229.RPOB"          "4029"             "0"            "11" 
    
    [[217]]
               name          length           frame          ncbicg 
    "CP015240.RPOB"          "4029"             "0"            "11" 
    
    [[218]]
               name          length           frame          ncbicg 
    "CP015241.RPOB"          "4029"             "0"            "11" 
    
    [[219]]
               name          length           frame          ncbicg 
    "CP015244.RPOB"          "4029"             "0"            "11" 
    
    [[220]]
               name          length           frame          ncbicg 
    "CP015834.RPOB"          "4029"             "0"            "11" 
    
    [[221]]
               name          length           frame          ncbicg 
    "CP015842.RPOB"          "4029"             "0"            "11" 
    
    [[222]]
               name          length           frame          ncbicg 
    "CP015843.RPOB"          "4029"             "0"            "11" 
    
    [[223]]
               name          length           frame          ncbicg 
    "CP015846.RPOB"          "4029"             "0"            "11" 
    
    [[224]]
               name          length           frame          ncbicg 
    "CP018801.RPOB"          "4029"             "0"            "11" 
    
    [[225]]
               name          length           frame          ncbicg 
    "CP019271.RPOB"          "4029"             "0"            "11" 
    
    [[226]]
               name          length           frame          ncbicg 
    "CP019944.RPOB"          "4029"             "0"            "11" 
    
    [[227]]
               name          length           frame          ncbicg 
    "CP020368.RPOB"          "4029"             "0"            "11" 
    
    [[228]]
               name          length           frame          ncbicg 
    "CP020509.RPOB"          "4029"             "0"            "11" 
    
    [[229]]
               name          length           frame          ncbicg 
    "CP020514.RPOB"          "4029"             "0"            "11" 
    
    [[230]]
               name          length           frame          ncbicg 
    "CP020516.RPOB"          "4029"             "0"            "11" 
    
    [[231]]
               name          length           frame          ncbicg 
    "CP020520.RPOB"          "4029"             "0"            "11" 
    
    [[232]]
               name          length           frame          ncbicg 
    "CP021179.RPOB"          "4029"             "0"            "11" 
    
    [[233]]
               name          length           frame          ncbicg 
    "CP021288.RPOB"          "4029"             "0"            "11" 
    
    [[234]]
               name          length           frame          ncbicg 
    "CP021335.RPOB"          "4029"             "0"            "11" 
    
    [[235]]
               name          length           frame          ncbicg 
    "CP021339.RPOB"          "4029"             "0"            "11" 
    
    [[236]]
               name          length           frame          ncbicg 
    "CP021454.RPOB"          "4029"             "0"            "11" 
    
    [[237]]
               name          length           frame          ncbicg 
    "CP021840.RPOB"          "4029"             "0"            "11" 
    
    [[238]]
               name          length           frame          ncbicg 
    "CP021844.RPOB"          "4029"             "0"            "11" 
    
    [[239]]
               name          length           frame          ncbicg 
    "CP022050.RPOB"          "4029"             "0"            "11" 
    
    [[240]]
               name          length           frame          ncbicg 
    "CP022086.RPOB"          "4029"             "0"            "11" 
    
    [[241]]
               name          length           frame          ncbicg 
    "CP022154.RPOB"          "4029"             "0"            "11" 
    
    [[242]]
               name          length           frame          ncbicg 
    "CP022164.RPOB"          "4029"             "0"            "11" 
    
    [[243]]
               name          length           frame          ncbicg 
    "CP022279.RPOB"          "4029"             "0"            "11" 
    
    [[244]]
               name          length           frame          ncbicg 
    "CP022393.RPOB"          "4029"             "0"            "11" 
    
    [[245]]
               name          length           frame          ncbicg 
    "CP022407.RPOB"          "4029"             "0"            "11" 
    
    [[246]]
               name          length           frame          ncbicg 
    "CP022466.RPOB"          "4029"             "0"            "11" 
    
    [[247]]
               name          length           frame          ncbicg 
    "CP022689.RPOB"          "4029"             "0"            "11" 
    
    [[248]]
               name          length           frame          ncbicg 
    "CP023142.RPOB"          "4029"             "0"            "11" 
    
    [[249]]
               name          length           frame          ncbicg 
    "CP023200.RPOB"          "4029"             "0"            "11" 
    
    [[250]]
               name          length           frame          ncbicg 
    "CP023201.RPOB"          "4029"             "0"            "11" 
    
    [[251]]
               name          length           frame          ncbicg 
    "CP023346.RPOB"          "4029"             "0"            "11" 
    
    [[252]]
               name          length           frame          ncbicg 
    "CP023349.RPOB"          "4029"             "0"            "11" 
    
    [[253]]
               name          length           frame          ncbicg 
    "CP023353.RPOB"          "4029"             "0"            "11" 
    
    [[254]]
               name          length           frame          ncbicg 
    "CP023357.RPOB"          "4029"             "0"            "11" 
    
    [[255]]
               name          length           frame          ncbicg 
    "CP023359.RPOB"          "4029"             "0"            "11" 
    
    [[256]]
               name          length           frame          ncbicg 
    "CP023364.RPOB"          "4029"             "0"            "11" 
    
    [[257]]
               name          length           frame          ncbicg 
    "CP023366.RPOB"          "4029"             "0"            "11" 
    
    [[258]]
               name          length           frame          ncbicg 
    "CP023371.RPOB"          "4029"             "0"            "11" 
    
    [[259]]
               name          length           frame          ncbicg 
    "CP023377.RPOB"          "4029"             "0"            "11" 
    
    [[260]]
               name          length           frame          ncbicg 
    "CP023383.RPOB"          "4029"             "0"            "11" 
    
    [[261]]
               name          length           frame          ncbicg 
    "CP023386.RPOB"          "4029"             "0"            "11" 
    
    [[262]]
               name          length           frame          ncbicg 
    "CP023388.RPOB"          "4029"             "0"            "11" 
    
    [[263]]
               name          length           frame          ncbicg 
    "CP023531.RPOB"          "4029"             "0"            "11" 
    
    [[264]]
               name          length           frame          ncbicg 
    "CP023535.RPOB"          "4029"             "0"            "11" 
    
    [[265]]
               name          length           frame          ncbicg 
    "CP023541.RPOB"          "4029"             "0"            "11" 
    
    [[266]]
               name          length           frame          ncbicg 
    "CP023644.RPOB"          "4029"             "0"            "11" 
    
    [[267]]
               name          length           frame          ncbicg 
    "CP023673.RPOB"          "4029"             "0"            "11" 
    
    [[268]]
               name          length           frame          ncbicg 
    "CP023870.RPOB"          "4029"             "0"            "11" 
    
    [[269]]
               name          length           frame          ncbicg 
    "CP023899.RPOB"          "4029"             "0"            "11" 
    
    [[270]]
               name          length           frame          ncbicg 
    "CP023960.RPOB"          "4029"             "0"            "11" 
    
    [[271]]
               name          length           frame          ncbicg 
    "CP024056.RPOB"          "4029"             "0"            "11" 
    
    [[272]]
               name          length           frame          ncbicg 
    "CP024090.RPOB"          "4029"             "0"            "11" 
    
    [[273]]
               name          length           frame          ncbicg 
    "CP024092.RPOB"          "4029"             "0"            "11" 
    
    [[274]]
               name          length           frame          ncbicg 
    "CP024127.RPOB"          "4029"             "0"            "11" 
    
    [[275]]
               name          length           frame          ncbicg 
    "CP024131.RPOB"          "4029"             "0"            "11" 
    
    [[276]]
               name          length           frame          ncbicg 
    "CP024134.RPOB"          "4029"             "0"            "11" 
    
    [[277]]
               name          length           frame          ncbicg 
    "CP024138.RPOB"          "4029"             "0"            "11" 
    
    [[278]]
               name          length           frame          ncbicg 
    "CP024141.RPOB"          "4029"             "0"            "11" 
    
    [[279]]
               name          length           frame          ncbicg 
    "CP024147.RPOB"          "4029"             "0"            "11" 
    
    [[280]]
               name          length           frame          ncbicg 
    "CP024155.RPOB"          "4029"             "0"            "11" 
    
    [[281]]
               name          length           frame          ncbicg 
    "CP024618.RPOB"          "4029"             "0"            "11" 
    
    [[282]]
               name          length           frame          ncbicg 
    "CP024650.RPOB"          "4029"             "0"            "11" 
    
    [[283]]
               name          length           frame          ncbicg 
    "CP024717.RPOB"          "4029"             "0"            "11" 
    
    [[284]]
               name          length           frame          ncbicg 
    "CP024720.RPOB"          "4029"             "0"            "11" 
    
    [[285]]
               name          length           frame          ncbicg 
    "CP024801.RPOB"          "4029"             "0"            "11" 
    
    [[286]]
               name          length           frame          ncbicg 
    "CP024815.RPOB"          "4029"             "0"            "11" 
    
    [[287]]
               name          length           frame          ncbicg 
    "CP024821.RPOB"          "4029"             "0"            "11" 
    
    [[288]]
               name          length           frame          ncbicg 
    "CP024826.RPOB"          "4029"             "0"            "11" 
    
    [[289]]
               name          length           frame          ncbicg 
    "CP024830.RPOB"          "4029"             "0"            "11" 
    
    [[290]]
               name          length           frame          ncbicg 
    "CP024851.RPOB"          "4029"             "0"            "11" 
    
    [[291]]
               name          length           frame          ncbicg 
    "CP024855.RPOB"          "4029"             "0"            "11" 
    
    [[292]]
               name          length           frame          ncbicg 
    "CP024859.RPOB"          "4029"             "0"            "11" 
    
    [[293]]
               name          length           frame          ncbicg 
    "CP024862.RPOB"          "4029"             "0"            "11" 
    
    [[294]]
               name          length           frame          ncbicg 
    "CP024886.RPOB"          "4029"             "0"            "11" 
    
    [[295]]
               name          length           frame          ncbicg 
    "CP024889.RPOB"          "4029"             "0"            "11" 
    
    [[296]]
               name          length           frame          ncbicg 
    "CP024978.RPOB"          "4029"             "0"            "11" 
    
    [[297]]
               name          length           frame          ncbicg 
    "CP024997.RPOB"          "4029"             "0"            "11" 
    
    [[298]]
               name          length           frame          ncbicg 
    "CP025036.RPOB"          "4029"             "0"            "11" 
    
    [[299]]
               name          length           frame          ncbicg 
    "CP025048.RPOB"          "4029"             "0"            "11" 
    
    [[300]]
               name          length           frame          ncbicg 
    "CP025251.RPOB"          "4029"             "0"            "11" 
    
    [[301]]
               name          length           frame          ncbicg 
    "CP025268.RPOB"          "4029"             "0"            "11" 
    
    [[302]]
               name          length           frame          ncbicg 
    "CP025328.RPOB"          "4029"             "0"            "11" 
    
    [[303]]
               name          length           frame          ncbicg 
    "CP025401.RPOB"          "4029"             "0"            "11" 
    
    [[304]]
               name          length           frame          ncbicg 
    "CP025627.RPOB"          "4029"             "0"            "11" 
    
    [[305]]
               name          length           frame          ncbicg 
    "CP025703.RPOB"          "4029"             "0"            "11" 
    
    [[306]]
               name          length           frame          ncbicg 
    "CP025707.RPOB"          "4029"             "0"            "11" 
    
    [[307]]
               name          length           frame          ncbicg 
    "CP025716.RPOB"          "4029"             "0"            "11" 
    
    [[308]]
               name          length           frame          ncbicg 
    "CP025739.RPOB"          "4029"             "0"            "11" 
    
    [[309]]
               name          length           frame          ncbicg 
    "CP025747.RPOB"          "4029"             "0"            "11" 
    
    [[310]]
               name          length           frame          ncbicg 
    "CP025753.RPOB"          "4029"             "0"            "11" 
    
    [[311]]
               name          length           frame          ncbicg 
    "CP025950.RPOB"          "4029"             "0"            "11" 
    
    [[312]]
               name          length           frame          ncbicg 
    "CP026085.RPOB"          "4029"             "0"            "11" 
    
    [[313]]
               name          length           frame          ncbicg 
    "CP026199.RPOB"          "4029"             "0"            "11" 
    
    [[314]]
               name          length           frame          ncbicg 
    "CP026202.RPOB"          "4029"             "0"            "11" 
    
    [[315]]
               name          length           frame          ncbicg 
    "CP026399.RPOB"          "4029"             "0"            "11" 
    
    [[316]]
               name          length           frame          ncbicg 
    "CP026491.RPOB"          "4029"             "0"            "11" 
    
    [[317]]
               name          length           frame          ncbicg 
    "CP026580.RPOB"          "4029"             "0"            "11" 
    
    [[318]]
               name          length           frame          ncbicg 
    "CP026612.RPOB"          "4029"             "0"            "11" 
    
    [[319]]
               name          length           frame          ncbicg 
    "CP026755.RPOB"          "4029"             "0"            "11" 
    
    [[320]]
               name          length           frame          ncbicg 
    "CP026853.RPOB"          "4029"             "0"            "11" 
    
    [[321]]
               name          length           frame          ncbicg 
    "CP027060.RPOB"          "4029"             "0"            "11" 
    
    [[322]]
               name          length           frame          ncbicg 
    "CP027105.RPOB"          "4029"             "0"            "11" 
    
    [[323]]
               name          length           frame          ncbicg 
    "CP027118.RPOB"          "4029"             "0"            "11" 
    
    [[324]]
               name          length           frame          ncbicg 
    "CP027126.RPOB"          "4029"             "0"            "11" 
    
    [[325]]
               name          length           frame          ncbicg 
    "CP027134.RPOB"          "4029"             "0"            "11" 
    
    [[326]]
               name          length           frame          ncbicg 
    "CP027140.RPOB"          "4029"             "0"            "11" 
    
    [[327]]
               name          length           frame          ncbicg 
    "CP027205.RPOB"          "4029"             "0"            "11" 
    
    [[328]]
               name          length           frame          ncbicg 
    "CP027255.RPOB"          "4029"             "0"            "11" 
    
    [[329]]
               name          length           frame          ncbicg 
    "CP027394.RPOB"          "4029"             "0"            "11" 
    
    [[330]]
               name          length           frame          ncbicg 
    "CP027430.RPOB"          "4029"             "0"            "11" 
    
    [[331]]
               name          length           frame          ncbicg 
    "CP027534.RPOB"          "4029"             "0"            "11" 
    
    [[332]]
               name          length           frame          ncbicg 
    "CP027701.RPOB"          "4029"             "0"            "11" 
    
    [[333]]
               name          length           frame          ncbicg 
    "CP028166.RPOB"          "4029"             "0"            "11" 
    
    [[334]]
               name          length           frame          ncbicg 
    "CP028192.RPOB"          "4029"             "0"            "11" 
    
    [[335]]
               name          length           frame          ncbicg 
    "CU651637.RPOB"          "4029"             "0"            "11" 
    
    [[336]]
               name          length           frame          ncbicg 
    "CU928145.RPOB"          "4029"             "0"            "11" 
    
    [[337]]
               name          length           frame          ncbicg 
    "CU928160.RPOB"          "4029"             "0"            "11" 
    
    [[338]]
               name          length           frame          ncbicg 
    "CU928161.RPOB"          "4029"             "0"            "11" 
    
    [[339]]
               name          length           frame          ncbicg 
    "CU928162.RPOB"          "4029"             "0"            "11" 
    
    [[340]]
               name          length           frame          ncbicg 
    "CU928163.RPOB"          "4029"             "0"            "11" 
    
    [[341]]
               name          length           frame          ncbicg 
    "CU928164.RPOB"          "4029"             "0"            "11" 
    
    [[342]]
              name         length          frame         ncbicg 
    "ECOUW89.RPOB"         "4029"            "0"           "11" 
    
    [[343]]
               name          length           frame          ncbicg 
    "ECU76222.RPOB"          "4029"             "0"            "11" 
    
    [[344]]
          name     length      frame     ncbicg 
    "ECU77436"      "512"        "0"       "11" 
    
    [[345]]
          name     length      frame     ncbicg 
    "EU010107"      "501"        "0"       "11" 
    
    [[346]]
               name          length           frame          ncbicg 
    "FM180568.RPOB"          "4029"             "0"            "11" 
    
    [[347]]
               name          length           frame          ncbicg 
    "FN554766.RPOB"          "4029"             "0"            "11" 
    
    [[348]]
          name     length      frame     ncbicg 
    "HE979716"      "550"        "0"       "11" 
    
    [[349]]
          name     length      frame     ncbicg 
    "HE979717"      "550"        "0"       "11" 
    
    [[350]]
          name     length      frame     ncbicg 
    "HE979721"      "550"        "0"       "11" 
    
    [[351]]
          name     length      frame     ncbicg 
    "HE979722"      "550"        "0"       "11" 
    
    [[352]]
          name     length      frame     ncbicg 
    "HE979723"      "550"        "0"       "11" 
    
    [[353]]
          name     length      frame     ncbicg 
    "HE979724"      "550"        "0"       "11" 
    
    [[354]]
          name     length      frame     ncbicg 
    "HE979725"      "550"        "0"       "11" 
    
    [[355]]
          name     length      frame     ncbicg 
    "HE979729"      "550"        "0"       "11" 
    
    [[356]]
          name     length      frame     ncbicg 
    "HE979730"      "550"        "0"       "11" 
    
    [[357]]
          name     length      frame     ncbicg 
    "HE979731"      "550"        "0"       "11" 
    
    [[358]]
          name     length      frame     ncbicg 
    "HE979732"      "550"        "0"       "11" 
    
    [[359]]
          name     length      frame     ncbicg 
    "HE979733"      "550"        "0"       "11" 
    
    [[360]]
          name     length      frame     ncbicg 
    "HE979734"      "550"        "0"       "11" 
    
    [[361]]
          name     length      frame     ncbicg 
    "HE979735"      "550"        "0"       "11" 
    
    [[362]]
          name     length      frame     ncbicg 
    "HE979736"      "550"        "0"       "11" 
    
    [[363]]
          name     length      frame     ncbicg 
    "HE979737"      "550"        "0"       "11" 
    
    [[364]]
               name          length           frame          ncbicg 
    "HF572917.RPOB"          "4029"             "0"            "11" 
    
    [[365]]
               name          length           frame          ncbicg 
    "HG428755.RPOB"          "4029"             "0"            "11" 
    
    [[366]]
               name          length           frame          ncbicg 
    "HG941718.RPOB"          "4062"             "0"            "11" 
    
    [[367]]
          name     length      frame     ncbicg 
    "JN707606"     "4029"        "0"       "11" 
    
    [[368]]
          name     length      frame     ncbicg 
    "JN707607"     "4029"        "0"       "11" 
    
    [[369]]
          name     length      frame     ncbicg 
    "JN707608"     "4029"        "0"       "11" 
    
    [[370]]
          name     length      frame     ncbicg 
    "JN707609"     "4029"        "0"       "11" 
    
    [[371]]
          name     length      frame     ncbicg 
    "JN707610"     "4029"        "0"       "11" 
    
    [[372]]
          name     length      frame     ncbicg 
    "JN707611"     "4029"        "0"       "11" 
    
    [[373]]
          name     length      frame     ncbicg 
    "JN707612"     "4029"        "0"       "11" 
    
    [[374]]
          name     length      frame     ncbicg 
    "JN707613"     "4029"        "0"       "11" 
    
    [[375]]
          name     length      frame     ncbicg 
    "JN707614"     "4029"        "0"       "11" 
    
    [[376]]
          name     length      frame     ncbicg 
    "JN707615"     "4029"        "0"       "11" 
    
    [[377]]
          name     length      frame     ncbicg 
    "JN707616"     "4029"        "0"       "11" 
    
    [[378]]
          name     length      frame     ncbicg 
    "JN707617"     "4029"        "0"       "11" 
    
    [[379]]
          name     length      frame     ncbicg 
    "JN707618"     "4029"        "0"       "11" 
    
    [[380]]
          name     length      frame     ncbicg 
    "JN707619"     "4029"        "0"       "11" 
    
    [[381]]
          name     length      frame     ncbicg 
    "JN707620"     "4029"        "0"       "11" 
    
    [[382]]
          name     length      frame     ncbicg 
    "JN707621"     "4029"        "0"       "11" 
    
    [[383]]
          name     length      frame     ncbicg 
    "JN707622"     "4029"        "0"       "11" 
    
    [[384]]
          name     length      frame     ncbicg 
    "JN707623"     "4029"        "0"       "11" 
    
    [[385]]
          name     length      frame     ncbicg 
    "JN707624"     "4029"        "0"       "11" 
    
    [[386]]
          name     length      frame     ncbicg 
    "JN707625"     "4029"        "0"       "11" 
    
    [[387]]
          name     length      frame     ncbicg 
    "JN707626"     "4029"        "0"       "11" 
    
    [[388]]
          name     length      frame     ncbicg 
    "JN707627"     "4029"        "0"       "11" 
    
    [[389]]
          name     length      frame     ncbicg 
    "JN707628"     "4029"        "0"       "11" 
    
    [[390]]
          name     length      frame     ncbicg 
    "JN707629"     "4029"        "0"       "11" 
    
    [[391]]
          name     length      frame     ncbicg 
    "JN707630"     "4029"        "0"       "11" 
    
    [[392]]
          name     length      frame     ncbicg 
    "JN707631"     "4029"        "0"       "11" 
    
    [[393]]
          name     length      frame     ncbicg 
    "JN707632"     "4029"        "0"       "11" 
    
    [[394]]
          name     length      frame     ncbicg 
    "JN707633"     "4029"        "0"       "11" 
    
    [[395]]
          name     length      frame     ncbicg 
    "JN707634"     "4029"        "0"       "11" 
    
    [[396]]
          name     length      frame     ncbicg 
    "JN707635"     "4029"        "0"       "11" 
    
    [[397]]
          name     length      frame     ncbicg 
    "JN707636"     "4029"        "0"       "11" 
    
    [[398]]
          name     length      frame     ncbicg 
    "JN707637"     "4029"        "0"       "11" 
    
    [[399]]
          name     length      frame     ncbicg 
    "JN707638"     "4029"        "0"       "11" 
    
    [[400]]
          name     length      frame     ncbicg 
    "JN707639"     "4029"        "0"       "11" 
    
    [[401]]
          name     length      frame     ncbicg 
    "JN707640"     "4029"        "0"       "11" 
    
    [[402]]
          name     length      frame     ncbicg 
    "JN707641"     "4029"        "0"       "11" 
    
    [[403]]
          name     length      frame     ncbicg 
    "JN707642"     "4029"        "0"       "11" 
    
    [[404]]
          name     length      frame     ncbicg 
    "JN707643"     "4029"        "0"       "11" 
    
    [[405]]
          name     length      frame     ncbicg 
    "JN707644"     "4029"        "0"       "11" 
    
    [[406]]
          name     length      frame     ncbicg 
    "JN707645"     "4029"        "0"       "11" 
    
    [[407]]
          name     length      frame     ncbicg 
    "JN707646"     "4029"        "0"       "11" 
    
    [[408]]
          name     length      frame     ncbicg 
    "JN707647"     "4029"        "0"       "11" 
    
    [[409]]
          name     length      frame     ncbicg 
    "JN707648"     "4029"        "0"       "11" 
    
    [[410]]
          name     length      frame     ncbicg 
    "JN707649"     "4029"        "0"       "11" 
    
    [[411]]
          name     length      frame     ncbicg 
    "JN707650"     "4029"        "0"       "11" 
    
    [[412]]
          name     length      frame     ncbicg 
    "JN707651"     "4029"        "0"       "11" 
    
    [[413]]
          name     length      frame     ncbicg 
    "JN707652"     "4029"        "0"       "11" 
    
    [[414]]
          name     length      frame     ncbicg 
    "JN707653"     "4029"        "0"       "11" 
    
    [[415]]
          name     length      frame     ncbicg 
    "JN707654"     "4029"        "0"       "11" 
    
    [[416]]
          name     length      frame     ncbicg 
    "JN707655"     "4029"        "0"       "11" 
    
    [[417]]
          name     length      frame     ncbicg 
    "JN707656"     "4029"        "0"       "11" 
    
    [[418]]
          name     length      frame     ncbicg 
    "JN707657"     "4029"        "0"       "11" 
    
    [[419]]
          name     length      frame     ncbicg 
    "JN707658"     "4029"        "0"       "11" 
    
    [[420]]
          name     length      frame     ncbicg 
    "JN707659"     "4029"        "0"       "11" 
    
    [[421]]
          name     length      frame     ncbicg 
    "JN707660"     "4029"        "0"       "11" 
    
    [[422]]
          name     length      frame     ncbicg 
    "JN707661"     "4029"        "0"       "11" 
    
    [[423]]
          name     length      frame     ncbicg 
    "JN707662"     "4029"        "0"       "11" 
    
    [[424]]
          name     length      frame     ncbicg 
    "JN707663"     "4029"        "0"       "11" 
    
    [[425]]
          name     length      frame     ncbicg 
    "JN707664"     "4029"        "0"       "11" 
    
    [[426]]
          name     length      frame     ncbicg 
    "JN707665"     "4029"        "0"       "11" 
    
    [[427]]
          name     length      frame     ncbicg 
    "JN707666"     "4029"        "0"       "11" 
    
    [[428]]
          name     length      frame     ncbicg 
    "JN707667"     "4029"        "0"       "11" 
    
    [[429]]
          name     length      frame     ncbicg 
    "JN707668"     "4029"        "0"       "11" 
    
    [[430]]
          name     length      frame     ncbicg 
    "JN707669"     "4029"        "0"       "11" 
    
    [[431]]
          name     length      frame     ncbicg 
    "JN707670"     "4029"        "0"       "11" 
    
    [[432]]
          name     length      frame     ncbicg 
    "JN707671"     "4029"        "0"       "11" 
    
    [[433]]
          name     length      frame     ncbicg 
    "JN707672"     "4029"        "0"       "11" 
    
    [[434]]
          name     length      frame     ncbicg 
    "JN707673"     "4029"        "0"       "11" 
    
    [[435]]
          name     length      frame     ncbicg 
    "JN707674"     "4029"        "0"       "11" 
    
    [[436]]
          name     length      frame     ncbicg 
    "JN707675"     "4029"        "0"       "11" 
    
    [[437]]
          name     length      frame     ncbicg 
    "JN707676"     "4029"        "0"       "11" 
    
    [[438]]
          name     length      frame     ncbicg 
    "JN707677"     "4029"        "0"       "11" 
    
    [[439]]
          name     length      frame     ncbicg 
    "JN707678"     "4029"        "0"       "11" 
    
    [[440]]
          name     length      frame     ncbicg 
    "JN707679"     "4029"        "0"       "11" 
    
    [[441]]
          name     length      frame     ncbicg 
    "JN707680"     "4029"        "0"       "11" 
    
    [[442]]
          name     length      frame     ncbicg 
    "JN707681"     "4029"        "0"       "11" 
    
    [[443]]
          name     length      frame     ncbicg 
    "JN707682"     "4029"        "0"       "11" 
    
    [[444]]
          name     length      frame     ncbicg 
    "JN707683"     "4029"        "0"       "11" 
    
    [[445]]
          name     length      frame     ncbicg 
    "JN811645"      "225"        "1"       "11" 
    
    [[446]]
          name     length      frame     ncbicg 
    "JN811646"      "201"        "1"       "11" 
    
    [[447]]
          name     length      frame     ncbicg 
    "JN811647"      "224"        "0"       "11" 
    
    [[448]]
          name     length      frame     ncbicg 
    "JN811648"      "217"        "2"       "11" 
    
    [[449]]
          name     length      frame     ncbicg 
    "JN811649"      "221"        "2"       "11" 
    
    [[450]]
          name     length      frame     ncbicg 
    "JN811650"      "218"        "0"       "11" 
    
    [[451]]
          name     length      frame     ncbicg 
    "JX471606"      "527"        "2"       "11" 
    
    [[452]]
          name     length      frame     ncbicg 
    "KJ186119"      "832"        "0"       "11" 
    
    [[453]]
          name     length      frame     ncbicg 
    "KP670778"     "4029"        "0"       "11" 
    
    [[454]]
          name     length      frame     ncbicg 
    "KP670779"     "4029"        "0"       "11" 
    
    [[455]]
          name     length      frame     ncbicg 
    "KP670780"     "4029"        "0"       "11" 
    
    [[456]]
          name     length      frame     ncbicg 
    "KP670781"     "4029"        "0"       "11" 
    
    [[457]]
          name     length      frame     ncbicg 
    "KP670782"     "4029"        "0"       "11" 
    
    [[458]]
          name     length      frame     ncbicg 
    "KP670783"     "4029"        "0"       "11" 
    
    [[459]]
          name     length      frame     ncbicg 
    "KP670784"     "4029"        "0"       "11" 
    
    [[460]]
          name     length      frame     ncbicg 
    "KP670785"     "4029"        "0"       "11" 
    
    [[461]]
          name     length      frame     ncbicg 
    "KP670786"     "4029"        "0"       "11" 
    
    [[462]]
          name     length      frame     ncbicg 
    "KP670787"     "4029"        "0"       "11" 
    
    [[463]]
          name     length      frame     ncbicg 
    "KP670788"     "4029"        "0"       "11" 
    
    [[464]]
          name     length      frame     ncbicg 
    "KP670789"     "4029"        "0"       "11" 
    
    [[465]]
          name     length      frame     ncbicg 
    "LC272544"      "971"        "2"       "11" 
    
    [[466]]
          name     length      frame     ncbicg 
    "LC272545"      "971"        "2"       "11" 
    
    [[467]]
          name     length      frame     ncbicg 
    "LC272546"      "971"        "2"       "11" 
    
    [[468]]
          name     length      frame     ncbicg 
    "LC272547"      "971"        "2"       "11" 
    
    [[469]]
          name     length      frame     ncbicg 
    "LC272548"      "971"        "2"       "11" 
    
    [[470]]
          name     length      frame     ncbicg 
    "LC272549"      "971"        "2"       "11" 
    
    [[471]]
          name     length      frame     ncbicg 
    "LC272550"      "971"        "2"       "11" 
    
    [[472]]
          name     length      frame     ncbicg 
    "LC272551"      "971"        "2"       "11" 
    
    [[473]]
          name     length      frame     ncbicg 
    "LC272552"      "971"        "2"       "11" 
    
    [[474]]
          name     length      frame     ncbicg 
    "LC272553"      "971"        "2"       "11" 
    
    [[475]]
          name     length      frame     ncbicg 
    "LC272554"      "971"        "2"       "11" 
    
    [[476]]
          name     length      frame     ncbicg 
    "LC272555"      "971"        "2"       "11" 
    
    [[477]]
          name     length      frame     ncbicg 
    "LC272556"      "971"        "2"       "11" 
    
    [[478]]
          name     length      frame     ncbicg 
    "LC272557"      "971"        "2"       "11" 
    
    [[479]]
          name     length      frame     ncbicg 
    "LC272558"      "971"        "2"       "11" 
    
    [[480]]
          name     length      frame     ncbicg 
    "LC272559"      "971"        "2"       "11" 
    
    [[481]]
          name     length      frame     ncbicg 
    "LC272560"      "971"        "2"       "11" 
    
    [[482]]
               name          length           frame          ncbicg 
    "LM651922.RPOB"           "953"             "2"            "11" 
    
    [[483]]
               name          length           frame          ncbicg 
    "LM651924.RPOB"           "878"             "2"            "11" 
    
    [[484]]
               name          length           frame          ncbicg 
    "LM993812.RPOB"          "4029"             "0"            "11" 
    
    [[485]]
               name          length           frame          ncbicg 
    "LM995446.RPOB"          "4029"             "0"            "11" 
    
    [[486]]
               name          length           frame          ncbicg 
    "LN832404.RPOB"          "4029"             "0"            "11" 
    
    [[487]]
               name          length           frame          ncbicg 
    "LT601384.RPOB"          "4029"             "0"            "11" 
    
    [[488]]
               name          length           frame          ncbicg 
    "LT827011.RPOB"          "4029"             "0"            "11" 
    
    [[489]]
               name          length           frame          ncbicg 
    "LT903847.RPOB"          "4029"             "0"            "11" 
    
    [[490]]
               name          length           frame          ncbicg 
    "LT906474.RPOB"          "4029"             "0"            "11" 
    
    [[491]]
             name        length         frame        ncbicg 
    "U00096.RPOB"        "4029"           "0"          "11" 
    
    [[492]]
             name        length         frame        ncbicg 
    "V00339.RPOB"        "4029"           "0"          "11" 
    



```R
myActino <- getSequence(q1$req[[658]])
myActino
# From these sequences, choose one sequence from each species (M. tuberculosis and E. coli) that will be compared. 
# In our case, we go with 'JX303316' in the former and 'AE005174.RPOB' in the latter species 
# (the sequence index numbers 658 and 1, respectively, in the corresponding qwa objects)
```


<ol class=list-inline>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
</ol>




```R
myProteo <- getSequence(q2$req[[1]])
myProteo
```


<ol class=list-inline>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'y'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
</ol>




```R
table(myActino)
# To compute the number of each base in the sequence, use the generic 'table' function
```


    myActino
       a    c    g    t 
     679 1084 1190  584 



```R
table(myProteo)
# To compute the number of each base in the sequence, use the generic 'table' function
length(myProteo)
```


    myProteo
       a    c    g    t    y 
    1013 1026 1090  899    1 



4029



```R
table(myProteo)/length(myProteo)
# ‘length’ function 
# To get the fraction for each base, do it in on a line by dividing the outcome of the 'table' function by the length of the sequence 
# (multiply it with 100 to get the result in a percentage)
```


    myProteo
               a            c            g            t            y 
    0.2514271531 0.2546537602 0.2705385952 0.2231322909 0.0002482005 



```R
myseq <- "AAAATGCAGTAACCCATGCCAAAATGCAGTAA"
myseq <- strsplit(myseq, "")
myseq <- unlist(myseq)
myseq

```


<ol class=list-inline>
	<li>'A'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'T'</li>
	<li>'G'</li>
	<li>'C'</li>
	<li>'A'</li>
	<li>'G'</li>
	<li>'T'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'C'</li>
	<li>'C'</li>
	<li>'C'</li>
	<li>'A'</li>
	<li>'T'</li>
	<li>'G'</li>
	<li>'C'</li>
	<li>'C'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'T'</li>
	<li>'G'</li>
	<li>'C'</li>
	<li>'A'</li>
	<li>'G'</li>
	<li>'T'</li>
	<li>'A'</li>
	<li>'A'</li>
</ol>




```R
table(myseq)
```


    myseq
     A  C  G  T 
    15  7  5  5 



```R
myseq <-
"MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWR
VISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFYLKMK
GDYFRYLSEVASGDNKQTTVSNQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNS
PEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTWTSENQGDEGENLYFQ"
myseq <- strsplit(myseq, "")
myseq <- unlist(myseq)
myseq
```


<ol class=list-inline>
	<li>'M'</li>
	<li>'T'</li>
	<li>'M'</li>
	<li>'D'</li>
	<li>'K'</li>
	<li>'S'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'V'</li>
	<li>'Q'</li>
	<li>'K'</li>
	<li>'A'</li>
	<li>'K'</li>
	<li>'L'</li>
	<li>'A'</li>
	<li>'E'</li>
	<li>'Q'</li>
	<li>'A'</li>
	<li>'E'</li>
	<li>'R'</li>
	<li>'Y'</li>
	<li>'D'</li>
	<li>'D'</li>
	<li>'M'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'M'</li>
	<li>'K'</li>
	<li>'A'</li>
	<li>'V'</li>
	<li>'T'</li>
	<li>'E'</li>
	<li>'Q'</li>
	<li>'G'</li>
	<li>'H'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'S'</li>
	<li>'N'</li>
	<li>'E'</li>
	<li>'E'</li>
	<li>'R'</li>
	<li>'N'</li>
	<li>'L'</li>
	<li>'L'</li>
	<li>'S'</li>
	<li>'V'</li>
	<li>'A'</li>
	<li>'Y'</li>
	<li>'K'</li>
	<li>'N'</li>
	<li>'V'</li>
	<li>'V'</li>
	<li>'G'</li>
	<li>'A'</li>
	<li>'R'</li>
	<li>'R'</li>
	<li>'S'</li>
	<li>'S'</li>
	<li>'W'</li>
	<li>'R'</li>
	<li>'\n'</li>
	<li>'V'</li>
	<li>'I'</li>
	<li>'S'</li>
	<li>'S'</li>
	<li>'I'</li>
	<li>'E'</li>
	<li>'Q'</li>
	<li>'K'</li>
	<li>'T'</li>
	<li>'E'</li>
	<li>'R'</li>
	<li>'N'</li>
	<li>'E'</li>
	<li>'K'</li>
	<li>'K'</li>
	<li>'Q'</li>
	<li>'Q'</li>
	<li>'M'</li>
	<li>'G'</li>
	<li>'K'</li>
	<li>'E'</li>
	<li>'Y'</li>
	<li>'R'</li>
	<li>'E'</li>
	<li>'K'</li>
	<li>'I'</li>
	<li>'E'</li>
	<li>'A'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'Q'</li>
	<li>'D'</li>
	<li>'I'</li>
	<li>'C'</li>
	<li>'N'</li>
	<li>'D'</li>
	<li>'V'</li>
	<li>'L'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'L'</li>
	<li>'D'</li>
	<li>'K'</li>
	<li>'Y'</li>
	<li>'L'</li>
	<li>'I'</li>
	<li>'P'</li>
	<li>'N'</li>
	<li>'A'</li>
	<li>'T'</li>
	<li>'Q'</li>
	<li>'P'</li>
	<li>'E'</li>
	<li>'S'</li>
	<li>'K'</li>
	<li>'V'</li>
	<li>'F'</li>
	<li>'Y'</li>
	<li>'L'</li>
	<li>'K'</li>
	<li>'M'</li>
	<li>'K'</li>
	<li>'\n'</li>
	<li>'G'</li>
	<li>'D'</li>
	<li>'Y'</li>
	<li>'F'</li>
	<li>'R'</li>
	<li>'Y'</li>
	<li>'L'</li>
	<li>'S'</li>
	<li>'E'</li>
	<li>'V'</li>
	<li>'A'</li>
	<li>'S'</li>
	<li>'G'</li>
	<li>'D'</li>
	<li>'N'</li>
	<li>'K'</li>
	<li>'Q'</li>
	<li>'T'</li>
	<li>'T'</li>
	<li>'V'</li>
	<li>'S'</li>
	<li>'N'</li>
	<li>'Q'</li>
	<li>'Q'</li>
	<li>'A'</li>
	<li>'Y'</li>
	<li>'Q'</li>
	<li>'E'</li>
	<li>'A'</li>
	<li>'F'</li>
	<li>'E'</li>
	<li>'I'</li>
	<li>'S'</li>
	<li>'K'</li>
	<li>'K'</li>
	<li>'E'</li>
	<li>'M'</li>
	<li>'Q'</li>
	<li>'P'</li>
	<li>'T'</li>
	<li>'H'</li>
	<li>'P'</li>
	<li>'I'</li>
	<li>'R'</li>
	<li>'L'</li>
	<li>'G'</li>
	<li>'L'</li>
	<li>'A'</li>
	<li>'L'</li>
	<li>'N'</li>
	<li>'F'</li>
	<li>'S'</li>
	<li>'V'</li>
	<li>'F'</li>
	<li>'Y'</li>
	<li>'Y'</li>
	<li>'E'</li>
	<li>'I'</li>
	<li>'L'</li>
	<li>'N'</li>
	<li>'S'</li>
	<li>'\n'</li>
	<li>'P'</li>
	<li>'E'</li>
	<li>'K'</li>
	<li>'A'</li>
	<li>'C'</li>
	<li>'S'</li>
	<li>'L'</li>
	<li>'A'</li>
	<li>'K'</li>
	<li>'T'</li>
	<li>'A'</li>
	<li>'F'</li>
	<li>'D'</li>
	<li>'E'</li>
	<li>'A'</li>
	<li>'I'</li>
	<li>'A'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'D'</li>
	<li>'T'</li>
	<li>'L'</li>
	<li>'N'</li>
	<li>'E'</li>
	<li>'E'</li>
	<li>'S'</li>
	<li>'Y'</li>
	<li>'K'</li>
	<li>'D'</li>
	<li>'S'</li>
	<li>'T'</li>
	<li>'L'</li>
	<li>'I'</li>
	<li>'M'</li>
	<li>'Q'</li>
	<li>'L'</li>
	<li>'L'</li>
	<li>'R'</li>
	<li>'D'</li>
	<li>'N'</li>
	<li>'L'</li>
	<li>'T'</li>
	<li>'W'</li>
	<li>'T'</li>
	<li>'S'</li>
	<li>'E'</li>
	<li>'N'</li>
	<li>'Q'</li>
	<li>'G'</li>
	<li>'D'</li>
	<li>'E'</li>
	<li>'G'</li>
	<li>'E'</li>
	<li>'N'</li>
	<li>'L'</li>
	<li>'Y'</li>
	<li>'F'</li>
	<li>'Q'</li>
</ol>




```R
myseq <-
"MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWR
VISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFYLKMK
GDYFRYLSEVASGDNKQTTVSNQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNS
PEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTWTSENQGDEGENLYFQ"
myseq <- s2c(myseq)
myseq
# 's2c' converses a string into a vector of chars.
# such as "BigBang" into a vector of chars such as c("B", "i", "g", "B", "a", "n", "g").
```


<ol class=list-inline>
	<li>'M'</li>
	<li>'T'</li>
	<li>'M'</li>
	<li>'D'</li>
	<li>'K'</li>
	<li>'S'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'V'</li>
	<li>'Q'</li>
	<li>'K'</li>
	<li>'A'</li>
	<li>'K'</li>
	<li>'L'</li>
	<li>'A'</li>
	<li>'E'</li>
	<li>'Q'</li>
	<li>'A'</li>
	<li>'E'</li>
	<li>'R'</li>
	<li>'Y'</li>
	<li>'D'</li>
	<li>'D'</li>
	<li>'M'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'A'</li>
	<li>'M'</li>
	<li>'K'</li>
	<li>'A'</li>
	<li>'V'</li>
	<li>'T'</li>
	<li>'E'</li>
	<li>'Q'</li>
	<li>'G'</li>
	<li>'H'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'S'</li>
	<li>'N'</li>
	<li>'E'</li>
	<li>'E'</li>
	<li>'R'</li>
	<li>'N'</li>
	<li>'L'</li>
	<li>'L'</li>
	<li>'S'</li>
	<li>'V'</li>
	<li>'A'</li>
	<li>'Y'</li>
	<li>'K'</li>
	<li>'N'</li>
	<li>'V'</li>
	<li>'V'</li>
	<li>'G'</li>
	<li>'A'</li>
	<li>'R'</li>
	<li>'R'</li>
	<li>'S'</li>
	<li>'S'</li>
	<li>'W'</li>
	<li>'R'</li>
	<li>'\n'</li>
	<li>'V'</li>
	<li>'I'</li>
	<li>'S'</li>
	<li>'S'</li>
	<li>'I'</li>
	<li>'E'</li>
	<li>'Q'</li>
	<li>'K'</li>
	<li>'T'</li>
	<li>'E'</li>
	<li>'R'</li>
	<li>'N'</li>
	<li>'E'</li>
	<li>'K'</li>
	<li>'K'</li>
	<li>'Q'</li>
	<li>'Q'</li>
	<li>'M'</li>
	<li>'G'</li>
	<li>'K'</li>
	<li>'E'</li>
	<li>'Y'</li>
	<li>'R'</li>
	<li>'E'</li>
	<li>'K'</li>
	<li>'I'</li>
	<li>'E'</li>
	<li>'A'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'Q'</li>
	<li>'D'</li>
	<li>'I'</li>
	<li>'C'</li>
	<li>'N'</li>
	<li>'D'</li>
	<li>'V'</li>
	<li>'L'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'L'</li>
	<li>'D'</li>
	<li>'K'</li>
	<li>'Y'</li>
	<li>'L'</li>
	<li>'I'</li>
	<li>'P'</li>
	<li>'N'</li>
	<li>'A'</li>
	<li>'T'</li>
	<li>'Q'</li>
	<li>'P'</li>
	<li>'E'</li>
	<li>'S'</li>
	<li>'K'</li>
	<li>'V'</li>
	<li>'F'</li>
	<li>'Y'</li>
	<li>'L'</li>
	<li>'K'</li>
	<li>'M'</li>
	<li>'K'</li>
	<li>'\n'</li>
	<li>'G'</li>
	<li>'D'</li>
	<li>'Y'</li>
	<li>'F'</li>
	<li>'R'</li>
	<li>'Y'</li>
	<li>'L'</li>
	<li>'S'</li>
	<li>'E'</li>
	<li>'V'</li>
	<li>'A'</li>
	<li>'S'</li>
	<li>'G'</li>
	<li>'D'</li>
	<li>'N'</li>
	<li>'K'</li>
	<li>'Q'</li>
	<li>'T'</li>
	<li>'T'</li>
	<li>'V'</li>
	<li>'S'</li>
	<li>'N'</li>
	<li>'Q'</li>
	<li>'Q'</li>
	<li>'A'</li>
	<li>'Y'</li>
	<li>'Q'</li>
	<li>'E'</li>
	<li>'A'</li>
	<li>'F'</li>
	<li>'E'</li>
	<li>'I'</li>
	<li>'S'</li>
	<li>'K'</li>
	<li>'K'</li>
	<li>'E'</li>
	<li>'M'</li>
	<li>'Q'</li>
	<li>'P'</li>
	<li>'T'</li>
	<li>'H'</li>
	<li>'P'</li>
	<li>'I'</li>
	<li>'R'</li>
	<li>'L'</li>
	<li>'G'</li>
	<li>'L'</li>
	<li>'A'</li>
	<li>'L'</li>
	<li>'N'</li>
	<li>'F'</li>
	<li>'S'</li>
	<li>'V'</li>
	<li>'F'</li>
	<li>'Y'</li>
	<li>'Y'</li>
	<li>'E'</li>
	<li>'I'</li>
	<li>'L'</li>
	<li>'N'</li>
	<li>'S'</li>
	<li>'\n'</li>
	<li>'P'</li>
	<li>'E'</li>
	<li>'K'</li>
	<li>'A'</li>
	<li>'C'</li>
	<li>'S'</li>
	<li>'L'</li>
	<li>'A'</li>
	<li>'K'</li>
	<li>'T'</li>
	<li>'A'</li>
	<li>'F'</li>
	<li>'D'</li>
	<li>'E'</li>
	<li>'A'</li>
	<li>'I'</li>
	<li>'A'</li>
	<li>'E'</li>
	<li>'L'</li>
	<li>'D'</li>
	<li>'T'</li>
	<li>'L'</li>
	<li>'N'</li>
	<li>'E'</li>
	<li>'E'</li>
	<li>'S'</li>
	<li>'Y'</li>
	<li>'K'</li>
	<li>'D'</li>
	<li>'S'</li>
	<li>'T'</li>
	<li>'L'</li>
	<li>'I'</li>
	<li>'M'</li>
	<li>'Q'</li>
	<li>'L'</li>
	<li>'L'</li>
	<li>'R'</li>
	<li>'D'</li>
	<li>'N'</li>
	<li>'L'</li>
	<li>'T'</li>
	<li>'W'</li>
	<li>'T'</li>
	<li>'S'</li>
	<li>'E'</li>
	<li>'N'</li>
	<li>'Q'</li>
	<li>'G'</li>
	<li>'D'</li>
	<li>'E'</li>
	<li>'G'</li>
	<li>'E'</li>
	<li>'N'</li>
	<li>'L'</li>
	<li>'Y'</li>
	<li>'F'</li>
	<li>'Q'</li>
</ol>



* There are two ways to make a string of protein or DNA(RNA) sequence into a list of single vectors:<br>
    * one is via `strsplit()`--->`unlist()`<br>
    * another is via `s2c()` function<br>
    * `c2s()` functions converting a vector of chars such as c("a", "d", "d") into a single string such as "add". 


```R
table(myseq)
```


    myseq
    \n  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y 
     3 20  2 13 29  7  8  2 10 20 24  8 14  5 16 10 18 12 11  2 12 



```R
library(graphics)
```


```R
Actinobacteria= GC(myActino)
Actinobacteria
```


0.642917726887193



```R
Proteobacteria=GC(myProteo)
Proteobacteria
```


0.5253227408143


To know the frequency of every possible pair of nucleotides in the sequence, use the `count` function, for every character pair (we can also do this for triples and so on by choosing the right value for the wordsize argument):


```R
seqinr::count(myActino, wordsize=2)
```


    
     aa  ac  ag  at  ca  cc  cg  ct  ga  gc  gg  gt  ta  tc  tg  tt 
    131 232 194 121 219 287 431 147 291 340 325 234  38 225 239  82 


The GC function of `seqinr` uses a similar method, but as an extension, it computes the fraction of `G` and `C` in the nucleotide. We can do this manually using the values of `C (index 2)` and `G (index 3)` as follows:


```R
myGC <- sum(table(myseq)[2], table(myseq)[3])/sum(table(myseq))
myGC
```


0.0894308943089431



```R

```
