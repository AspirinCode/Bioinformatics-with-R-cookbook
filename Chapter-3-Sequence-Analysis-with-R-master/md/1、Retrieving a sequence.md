
# Retrieving a sequence
Before we do sequence analysis, the first thing we need is a sequence of DNA or protein. We can retrieve such sequence data either by visiting the database hosting page via a browser or by accessing the data mart form within R via corresponding commands/functions.


```R
source("http://www.bioconductor.org/biocLite.R")
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
biocLite("seqinr")
```

    BioC_mirror: https://bioconductor.org
    Using Bioconductor 3.6 (BiocInstaller 1.28.0), R 3.4.2 (2017-09-28).
    Installing package(s) 'seqinr'
    also installing the dependencies 'MASS', 'ade4', 'segmented'
    
    

    package 'MASS' successfully unpacked and MD5 sums checked
    package 'ade4' successfully unpacked and MD5 sums checked
    package 'segmented' successfully unpacked and MD5 sums checked
    package 'seqinr' successfully unpacked and MD5 sums checked
    
    The downloaded binary packages are in
    	C:\Users\Master1\AppData\Local\Temp\RtmpspRJxw\downloaded_packages
    

    Old packages: 'GenomicRanges', 'matrixStats', 'tibble', 'digest', 'pbdZMQ',
      'stringi'
    


```R
library(seqinr)
```


```R
choosebank()
```


     [1] "genbank"         "embl"            "emblwgs"         "swissprot"      
     [5] "ensembl"         "hogenom7"        "hogenom"         "hogenomdna"     
     [9] "hovergendna"     "hovergen"        "hogenom5"        "hogenom5dna"    
    [13] "hogenom4"        "hogenom4dna"     "homolens"        "homolensdna"    
    [17] "hobacnucl"       "hobacprot"       "phever2"         "phever2dna"     
    [21] "refseq"          "greviews"        "bacterial"       "archaeal"       
    [25] "protozoan"       "ensprotists"     "ensfungi"        "ensmetazoa"     
    [29] "ensplants"       "ensemblbacteria" "mito"            "polymorphix"    
    [33] "emglib"          "refseqViruses"   "ribodb"          "taxodb"         



```R
choosebank("genbank")
```


```R
q1 <- query("BRCA1", "SP=Homo sapiens AND K=BRCA1")
q1
```


    175 SQ for SP=Homo sapiens AND K=BRCA1


* When we have the data bank open, we always query the bank by using the `query` command.


```R
?query
```


```R
attributes(q1)
# Have a look at the returned value of the `query` command by looking at the various components of the object
```


<dl>
	<dt>$names</dt>
		<dd><ol class=list-inline>
	<li>'call'</li>
	<li>'name'</li>
	<li>'nelem'</li>
	<li>'typelist'</li>
	<li>'req'</li>
	<li>'socket'</li>
</ol>
</dd>
	<dt>$class</dt>
		<dd>'qaw'</dd>
</dl>




```R
q1$req
# To check the set of all of the sequences that were retrieved, choose the category and get the name, length, and 
# other attributes for every available sequence in the query
```


    [[1]]
                name           length            frame           ncbicg 
    "AB621825.BRCA1"             "71"              "0"              "1" 
    
    [[2]]
                name           length            frame           ncbicg 
    "AF005068.BRCA1"           "5379"              "0"              "1" 
    
    [[3]]
                name           length            frame           ncbicg 
    "AF284812.BRCA1"             "84"              "0"              "1" 
    
    [[4]]
                name           length            frame           ncbicg 
    "AF507076.BRCA1"             "84"              "0"              "1" 
    
    [[5]]
                name           length            frame           ncbicg 
    "AF507077.BRCA1"             "84"              "0"              "1" 
    
    [[6]]
                name           length            frame           ncbicg 
    "AF507078.BRCA1"             "84"              "0"              "1" 
    
    [[7]]
                name           length            frame           ncbicg 
    "AY093484.BRCA1"             "46"              "1"              "1" 
    
    [[8]]
                name           length            frame           ncbicg 
    "AY093485.BRCA1"             "84"              "0"              "1" 
    
    [[9]]
                name           length            frame           ncbicg 
    "AY093486.BRCA1"             "84"              "0"              "1" 
    
    [[10]]
                name           length            frame           ncbicg 
    "AY093487.BRCA1"             "84"              "0"              "1" 
    
    [[11]]
                name           length            frame           ncbicg 
    "AY093489.BRCA1"             "84"              "0"              "1" 
    
    [[12]]
                name           length            frame           ncbicg 
    "AY093490.BRCA1"             "84"              "0"              "1" 
    
    [[13]]
                name           length            frame           ncbicg 
    "AY093491.BRCA1"             "84"              "0"              "1" 
    
    [[14]]
                name           length            frame           ncbicg 
    "AY093492.BRCA1"             "84"              "0"              "1" 
    
    [[15]]
                name           length            frame           ncbicg 
    "AY093493.BRCA1"             "84"              "0"              "1" 
    
    [[16]]
          name     length      frame     ncbicg 
    "AY144588"       "68"        "2"        "1" 
    
    [[17]]
          name     length      frame     ncbicg 
    "AY150865"       "68"        "2"        "1" 
    
    [[18]]
                name           length            frame           ncbicg 
    "AY273801.BRCA1"           "5592"              "0"              "1" 
    
    [[19]]
                name           length            frame           ncbicg 
    "AY304547.BRCA1"           "3426"              "2"              "1" 
    
    [[20]]
          name     length      frame     ncbicg 
    "AY438030"      "951"        "0"        "1" 
    
    [[21]]
          name     length      frame     ncbicg 
    "AY438031"     "1152"        "0"        "1" 
    
    [[22]]
          name     length      frame     ncbicg 
    "AY706911"      "272"        "0"        "1" 
    
    [[23]]
                name           length            frame           ncbicg 
    "AY706912.BRCA1"             "44"              "2"              "1" 
    
    [[24]]
          name     length      frame     ncbicg 
    "AY706913"       "50"        "2"        "1" 
    
    [[25]]
          name     length      frame     ncbicg 
    "AY751490"     "5524"        "0"        "1" 
    
    [[26]]
                name           length            frame           ncbicg 
    "BC030969.BRCA1"           "1971"              "0"              "1" 
    
    [[27]]
          name     length      frame     ncbicg 
    "BC062429"     "1322"        "2"        "1" 
    
    [[28]]
                name           length            frame           ncbicg 
    "BC072418.BRCA1"           "2100"              "0"              "1" 
    
    [[29]]
                name           length            frame           ncbicg 
    "BC085615.BRCA1"           "1872"              "0"              "1" 
    
    [[30]]
                name           length            frame           ncbicg 
    "BC106745.BRCA1"           "1419"              "0"              "1" 
    
    [[31]]
          name     length      frame     ncbicg 
    "BC106746"      "779"        "2"        "1" 
    
    [[32]]
                name           length            frame           ncbicg 
    "BC115037.BRCA1"           "4065"              "0"              "1" 
    
    [[33]]
          name     length      frame     ncbicg 
    "DQ075361"       "50"        "2"        "1" 
    
    [[34]]
          name     length      frame     ncbicg 
    "DQ116737"      "498"        "0"        "1" 
    
    [[35]]
          name     length      frame     ncbicg 
    "DQ145822"      "120"        "1"        "1" 
    
    [[36]]
          name     length      frame     ncbicg 
    "DQ145823"      "101"        "0"        "1" 
    
    [[37]]
          name     length      frame     ncbicg 
    "DQ145824"      "114"        "1"        "1" 
    
    [[38]]
                name           length            frame           ncbicg 
    "DQ145825.BRCA1"            "106"              "1"              "1" 
    
    [[39]]
          name     length      frame     ncbicg 
    "DQ145826"      "120"        "2"        "1" 
    
    [[40]]
                name           length            frame           ncbicg 
    "DQ190450.BRCA1"           "5592"              "0"              "1" 
    
    [[41]]
                name           length            frame           ncbicg 
    "DQ190451.BRCA1"           "5592"              "0"              "1" 
    
    [[42]]
                name           length            frame           ncbicg 
    "DQ190452.BRCA1"           "5592"              "0"              "1" 
    
    [[43]]
                name           length            frame           ncbicg 
    "DQ190453.BRCA1"           "5592"              "0"              "1" 
    
    [[44]]
                name           length            frame           ncbicg 
    "DQ190454.BRCA1"           "5592"              "0"              "1" 
    
    [[45]]
                name           length            frame           ncbicg 
    "DQ190455.BRCA1"           "5592"              "0"              "1" 
    
    [[46]]
                name           length            frame           ncbicg 
    "DQ190456.BRCA1"           "5592"              "0"              "1" 
    
    [[47]]
                name           length            frame           ncbicg 
    "DQ190457.BRCA1"           "5467"              "0"              "1" 
    
    [[48]]
          name     length      frame     ncbicg 
    "DQ299305"      "120"        "1"        "1" 
    
    [[49]]
          name     length      frame     ncbicg 
    "DQ299306"      "120"        "1"        "1" 
    
    [[50]]
          name     length      frame     ncbicg 
    "DQ299307"      "120"        "1"        "1" 
    
    [[51]]
          name     length      frame     ncbicg 
    "DQ299308"      "120"        "1"        "1" 
    
    [[52]]
          name     length      frame     ncbicg 
    "DQ299309"      "120"        "1"        "1" 
    
    [[53]]
          name     length      frame     ncbicg 
    "DQ299310"      "120"        "1"        "1" 
    
    [[54]]
          name     length      frame     ncbicg 
    "DQ299311"      "120"        "1"        "1" 
    
    [[55]]
          name     length      frame     ncbicg 
    "DQ299312"      "120"        "1"        "1" 
    
    [[56]]
          name     length      frame     ncbicg 
    "DQ299313"      "120"        "1"        "1" 
    
    [[57]]
          name     length      frame     ncbicg 
    "DQ299314"      "120"        "1"        "1" 
    
    [[58]]
          name     length      frame     ncbicg 
    "DQ299315"      "120"        "1"        "1" 
    
    [[59]]
          name     length      frame     ncbicg 
    "DQ299316"      "101"        "0"        "1" 
    
    [[60]]
          name     length      frame     ncbicg 
    "DQ299317"      "101"        "0"        "1" 
    
    [[61]]
          name     length      frame     ncbicg 
    "DQ299318"      "101"        "0"        "1" 
    
    [[62]]
          name     length      frame     ncbicg 
    "DQ299319"      "101"        "0"        "1" 
    
    [[63]]
          name     length      frame     ncbicg 
    "DQ299320"      "114"        "1"        "1" 
    
    [[64]]
          name     length      frame     ncbicg 
    "DQ299321"      "120"        "2"        "1" 
    
    [[65]]
          name     length      frame     ncbicg 
    "DQ299322"      "120"        "2"        "1" 
    
    [[66]]
          name     length      frame     ncbicg 
    "DQ299323"      "120"        "2"        "1" 
    
    [[67]]
          name     length      frame     ncbicg 
    "DQ299324"      "120"        "2"        "1" 
    
    [[68]]
          name     length      frame     ncbicg 
    "DQ299325"      "120"        "2"        "1" 
    
    [[69]]
          name     length      frame     ncbicg 
    "DQ299326"      "120"        "2"        "1" 
    
    [[70]]
          name     length      frame     ncbicg 
    "DQ299327"      "120"        "2"        "1" 
    
    [[71]]
          name     length      frame     ncbicg 
    "DQ299328"      "120"        "2"        "1" 
    
    [[72]]
          name     length      frame     ncbicg 
    "DQ299330"      "120"        "2"        "1" 
    
    [[73]]
          name     length      frame     ncbicg 
    "DQ299331"      "120"        "2"        "1" 
    
    [[74]]
                name           length            frame           ncbicg 
    "DQ363751.BRCA1"           "1065"              "0"              "1" 
    
    [[75]]
                name           length            frame           ncbicg 
    "DQ478408.BRCA1"           "5467"              "0"              "1" 
    
    [[76]]
                name           length            frame           ncbicg 
    "FJ940752.BRCA1"             "75"              "2"              "1" 
    
    [[77]]
          name     length      frame     ncbicg 
    "HE600032"      "137"        "0"        "1" 
    
    [[78]]
                name           length            frame           ncbicg 
    "HE600033.BRCA1"            "150"              "0"              "1" 
    
    [[79]]
          name     length      frame     ncbicg 
    "HE600034"      "139"        "2"        "1" 
    
    [[80]]
          name     length      frame     ncbicg 
    "HE600035"      "139"        "2"        "1" 
    
    [[81]]
          name     length      frame     ncbicg 
    "HE600036"      "139"        "2"        "1" 
    
    [[82]]
          name     length      frame     ncbicg 
    "HE600037"      "132"        "1"        "1" 
    
    [[83]]
          name     length      frame     ncbicg 
    "HE600038"      "140"        "0"        "1" 
    
    [[84]]
                name           length            frame           ncbicg 
    "HSU14680.BRCA1"           "5592"              "0"              "1" 
    
    [[85]]
          name     length      frame     ncbicg 
    "HSU18009"     "2381"        "0"        "1" 
    
    [[86]]
          name     length      frame     ncbicg 
    "HSU18018"     "2333"        "0"        "1" 
    
    [[87]]
                name           length            frame           ncbicg 
    "HSU37574.BRCA1"             "80"              "0"              "1" 
    
    [[88]]
                name           length            frame           ncbicg 
    "HSU61268.BRCA1"             "80"              "0"              "1" 
    
    [[89]]
          name     length      frame     ncbicg 
    "HSU64805"     "2280"        "0"        "1" 
    
    [[90]]
                name           length            frame           ncbicg 
    "HSU68041.BRCA1"            "702"              "0"              "1" 
    
    [[91]]
                name           length            frame           ncbicg 
    "JN384124.BRCA1"            "259"              "0"              "1" 
    
    [[92]]
          name     length      frame     ncbicg 
    "JN686490"     "5592"        "0"        "1" 
    
    [[93]]
                name           length            frame           ncbicg 
    "JX480460.BRCA1"             "78"              "0"              "1" 
    
    [[94]]
          name     length      frame     ncbicg 
    "JX480461"      "112"        "0"        "1" 
    
    [[95]]
                name           length            frame           ncbicg 
    "JX480462.BRCA1"            "228"              "0"              "1" 
    
    [[96]]
                name           length            frame           ncbicg 
    "JX480463.BRCA1"             "63"              "0"              "1" 
    
    [[97]]
                name           length            frame           ncbicg 
    "JX480464.BRCA1"            "201"              "0"              "1" 
    
    [[98]]
                name           length            frame           ncbicg 
    "JX480465.BRCA1"             "42"              "0"              "1" 
    
    [[99]]
                name           length            frame           ncbicg 
    "JX480466.BRCA1"            "180"              "0"              "1" 
    
    [[100]]
          name     length      frame     ncbicg 
    "JX480467"       "90"        "0"        "1" 
    
    [[101]]
                name           length            frame           ncbicg 
    "KJ625149.BRCA1"             "36"              "0"              "1" 
    
    [[102]]
                name           length            frame           ncbicg 
    "KJ625150.BRCA1"             "81"              "0"              "1" 
    
    [[103]]
                name           length            frame           ncbicg 
    "KJ625151.BRCA1"             "70"              "1"              "1" 
    
    [[104]]
                name           length            frame           ncbicg 
    "KJ625152.BRCA1"             "30"              "0"              "1" 
    
    [[105]]
                name           length            frame           ncbicg 
    "KJ625153.BRCA1"             "63"              "0"              "1" 
    
    [[106]]
                name           length            frame           ncbicg 
    "KJ625154.BRCA1"             "75"              "0"              "1" 
    
    [[107]]
                name           length            frame           ncbicg 
    "KJ625155.BRCA1"            "281"              "2"              "1" 
    
    [[108]]
                name           length            frame           ncbicg 
    "KJ625156.BRCA1"            "150"              "0"              "1" 
    
    [[109]]
                name           length            frame           ncbicg 
    "KJ625157.BRCA1"            "114"              "0"              "1" 
    
    [[110]]
                name           length            frame           ncbicg 
    "KJ625158.BRCA1"            "150"              "0"              "1" 
    
    [[111]]
                name           length            frame           ncbicg 
    "KJ625159.BRCA1"            "248"              "2"              "1" 
    
    [[112]]
                name           length            frame           ncbicg 
    "KJ625160.BRCA1"            "322"              "1"              "1" 
    
    [[113]]
                name           length            frame           ncbicg 
    "KJ625161.BRCA1"            "280"              "1"              "1" 
    
    [[114]]
          name     length      frame     ncbicg 
    "KJ625162"      "422"        "1"        "1" 
    
    [[115]]
                name           length            frame           ncbicg 
    "KJ625163.BRCA1"            "345"              "0"              "1" 
    
    [[116]]
                name           length            frame           ncbicg 
    "KJ625164.BRCA1"            "116"              "2"              "1" 
    
    [[117]]
                name           length            frame           ncbicg 
    "KJ625165.BRCA1"            "254"              "2"              "1" 
    
    [[118]]
                name           length            frame           ncbicg 
    "KJ625166.BRCA1"            "263"              "2"              "1" 
    
    [[119]]
                name           length            frame           ncbicg 
    "KJ625167.BRCA1"            "383"              "2"              "1" 
    
    [[120]]
                name           length            frame           ncbicg 
    "KJ625168.BRCA1"            "347"              "2"              "1" 
    
    [[121]]
                name           length            frame           ncbicg 
    "KJ625169.BRCA1"            "273"              "0"              "1" 
    
    [[122]]
                name           length            frame           ncbicg 
    "KJ625170.BRCA1"             "53"              "2"              "1" 
    
    [[123]]
                name           length            frame           ncbicg 
    "KJ625171.BRCA1"             "17"              "2"              "1" 
    
    [[124]]
                name           length            frame           ncbicg 
    "KJ625172.BRCA1"            "173"              "2"              "1" 
    
    [[125]]
                name           length            frame           ncbicg 
    "KJ625173.BRCA1"             "23"              "2"              "1" 
    
    [[126]]
                name           length            frame           ncbicg 
    "KJ625174.BRCA1"             "38"              "2"              "1" 
    
    [[127]]
                name           length            frame           ncbicg 
    "KJ625175.BRCA1"            "112"              "1"              "1" 
    
    [[128]]
                name           length            frame           ncbicg 
    "KJ625176.BRCA1"            "143"              "1"              "1" 
    
    [[129]]
              name         length          frame         ncbicg 
    "KJ625176.PE2"          "112"            "1"            "1" 
    
    [[130]]
                name           length            frame           ncbicg 
    "KJ625177.BRCA1"             "83"              "2"              "1" 
    
    [[131]]
                name           length            frame           ncbicg 
    "KJ625178.BRCA1"            "116"              "0"              "1" 
    
    [[132]]
                name           length            frame           ncbicg 
    "KJ625179.BRCA1"            "116"              "0"              "1" 
    
    [[133]]
          name     length      frame     ncbicg 
    "KM434065"      "577"        "0"        "1" 
    
    [[134]]
                name           length            frame           ncbicg 
    "KP255396.BRCA1"            "102"              "0"              "1" 
    
    [[135]]
                name           length            frame           ncbicg 
    "KP255397.BRCA1"            "122"              "2"              "1" 
    
    [[136]]
                name           length            frame           ncbicg 
    "KP255398.BRCA1"            "221"              "2"              "1" 
    
    [[137]]
                name           length            frame           ncbicg 
    "KP255399.BRCA1"            "356"              "2"              "1" 
    
    [[138]]
                name           length            frame           ncbicg 
    "KP255400.BRCA1"             "91"              "1"              "1" 
    
    [[139]]
                name           length            frame           ncbicg 
    "KP255401.BRCA1"            "125"              "2"              "1" 
    
    [[140]]
                name           length            frame           ncbicg 
    "KP255402.BRCA1"            "296"              "2"              "1" 
    
    [[141]]
                name           length            frame           ncbicg 
    "KP255403.BRCA1"            "175"              "1"              "1" 
    
    [[142]]
                name           length            frame           ncbicg 
    "KP272102.BRCA1"            "143"              "1"              "1" 
    
    [[143]]
                name           length            frame           ncbicg 
    "KP272103.BRCA1"             "88"              "1"              "1" 
    
    [[144]]
                name           length            frame           ncbicg 
    "KP272104.BRCA1"            "170"              "2"              "1" 
    
    [[145]]
                name           length            frame           ncbicg 
    "KP272105.BRCA1"            "421"              "1"              "1" 
    
    [[146]]
                name           length            frame           ncbicg 
    "KP272106.BRCA1"            "227"              "2"              "1" 
    
    [[147]]
          name     length      frame     ncbicg 
    "KP404097"      "630"        "0"        "1" 
    
    [[148]]
                name           length            frame           ncbicg 
    "KP455327.BRCA1"            "172"              "0"              "1" 
    
    [[149]]
          name     length      frame     ncbicg 
    "KP701015"      "327"        "0"        "1" 
    
    [[150]]
                name           length            frame           ncbicg 
    "KP701016.BRCA1"            "126"              "0"              "1" 
    
    [[151]]
          name     length      frame     ncbicg 
    "KP729136"      "717"        "0"        "1" 
    
    [[152]]
          name     length      frame     ncbicg 
    "KP729137"      "311"        "2"        "1" 
    
    [[153]]
          name     length      frame     ncbicg 
    "KP744861"      "777"        "0"        "1" 
    
    [[154]]
                name           length            frame           ncbicg 
    "KP753383.BRCA1"             "89"              "2"              "1" 
    
    [[155]]
                name           length            frame           ncbicg 
    "KT120061.BRCA1"             "80"              "0"              "1" 
    
    [[156]]
          name     length      frame     ncbicg 
    "KT152888"      "324"        "0"        "1" 
    
    [[157]]
                name           length            frame           ncbicg 
    "KT152889.BRCA1"            "176"              "2"              "1" 
    
    [[158]]
                name           length            frame           ncbicg 
    "KT152890.BRCA1"           "2223"              "0"              "1" 
    
    [[159]]
                name           length            frame           ncbicg 
    "KT844468.BRCA1"             "66"              "0"              "1" 
    
    [[160]]
                name           length            frame           ncbicg 
    "KT844469.BRCA1"             "89"              "1"              "1" 
    
    [[161]]
          name     length      frame     ncbicg 
    "KU359055"       "88"        "0"        "1" 
    
    [[162]]
          name     length      frame     ncbicg 
    "KU359056"       "88"        "0"        "1" 
    
    [[163]]
          name     length      frame     ncbicg 
    "KU359057"       "88"        "0"        "1" 
    
    [[164]]
          name     length      frame     ncbicg 
    "KU359058"       "88"        "0"        "1" 
    
    [[165]]
          name     length      frame     ncbicg 
    "KU359059"       "88"        "0"        "1" 
    
    [[166]]
          name     length      frame     ncbicg 
    "KU359060"       "88"        "0"        "1" 
    
    [[167]]
          name     length      frame     ncbicg 
    "KU359061"       "88"        "0"        "1" 
    
    [[168]]
          name     length      frame     ncbicg 
    "KU359062"       "88"        "0"        "1" 
    
    [[169]]
          name     length      frame     ncbicg 
    "KU359063"       "88"        "0"        "1" 
    
    [[170]]
          name     length      frame     ncbicg 
    "KU359064"       "88"        "0"        "1" 
    
    [[171]]
                name           length            frame           ncbicg 
    "KX580312.BRCA1"            "127"              "1"              "1" 
    
    [[172]]
                name           length            frame           ncbicg 
    "KX944478.BRCA1"           "1025"              "2"              "1" 
    
    [[173]]
              name         length          frame         ncbicg 
    "L78833.BRCA1"         "5592"            "0"            "1" 
    
    [[174]]
              name         length          frame         ncbicg 
    "S78558.BRCA1"           "66"            "0"            "1" 
    
    [[175]]
              name         length          frame         ncbicg 
    "Y08757.BRCA1"           "88"            "0"            "1" 
    



```R
q2 <- query("BRCA1", "SP=Homo sapiens AND AC=U61268")
q2
# To fetch a specific sequence by its accession number, use the `AC` attribute in the `query` command
```


    1 SQ for SP=Homo sapiens AND AC=U61268



```R
myseq <- getSequence(q2$req[[1]])
myseq
# To fetch a specific sequence from the `query` object, use the `getsequence` command
```


<ol class=list-inline>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
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
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
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
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
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
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'m'</li>
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
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
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
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
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
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'g'</li>
	<li>'t'</li>
	<li>'g'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'c'</li>
	<li>'c'</li>
	<li>'a'</li>
	<li>'c'</li>
	<li>'t'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'a'</li>
	<li>'t'</li>
</ol>




```R
annots <- getAnnot(q2$req[[1]])
annots
# fetch the annotation (the other attributes of a sequence) 
```


<ol class=list-inline>
	<li><span style=white-space:pre-wrap>'LOCUS       HSU61268                1338 bp    DNA     linear   PRI 16-JAN-1997'</span></li>
	<li><span style=white-space:pre-wrap>'DEFINITION  Human breast and ovarian cancer susceptibility (BRCA1) gene, exon'</span></li>
	<li><span style=white-space:pre-wrap>'            2, partial flanking introns, and partial cds.'</span></li>
	<li><span style=white-space:pre-wrap>'ACCESSION   U61268'</span></li>
	<li><span style=white-space:pre-wrap>'VERSION     U61268.1'</span></li>
	<li><span style=white-space:pre-wrap>'KEYWORDS    .'</span></li>
	<li><span style=white-space:pre-wrap>'SOURCE      Homo sapiens (human)'</span></li>
	<li><span style=white-space:pre-wrap>'  ORGANISM  Homo sapiens'</span></li>
	<li><span style=white-space:pre-wrap>'            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;'</span></li>
	<li><span style=white-space:pre-wrap>'            Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;'</span></li>
	<li><span style=white-space:pre-wrap>'            Catarrhini; Hominidae; Homo.'</span></li>
	<li><span style=white-space:pre-wrap>'REFERENCE   1  (bases 1 to 1338)'</span></li>
	<li><span style=white-space:pre-wrap>'  AUTHORS   Struewing,J.P., Robbins,C. and Brody,L.C.'</span></li>
	<li><span style=white-space:pre-wrap>'  TITLE     BRCA1: genomic sequence flanking exon 2'</span></li>
	<li><span style=white-space:pre-wrap>'  JOURNAL   Unpublished'</span></li>
	<li><span style=white-space:pre-wrap>'REFERENCE   2  (bases 1 to 1338)'</span></li>
	<li><span style=white-space:pre-wrap>'  AUTHORS   Struewing,J.P., Robbins,C. and Brody,L.C.'</span></li>
	<li><span style=white-space:pre-wrap>'  TITLE     Direct Submission'</span></li>
	<li><span style=white-space:pre-wrap>'  JOURNAL   Submitted (19-JUN-1996) NCI/NCHGR, NIH, 6130 Executive Blvd MSC'</span></li>
	<li><span style=white-space:pre-wrap>'            7372, Bethesda, MD 20892-7372, USA'</span></li>
	<li><span style=white-space:pre-wrap>'FEATURES             Location/Qualifiers'</span></li>
	<li><span style=white-space:pre-wrap>'     source          1..1338'</span></li>
	<li><span style=white-space:pre-wrap>'                     /organism="Homo sapiens"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /mol_type="genomic DNA"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /db_xref="taxon:9606"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /chromosome="17"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /map="17q21"'</span></li>
	<li><span style=white-space:pre-wrap>'     repeat_region   1..84'</span></li>
	<li><span style=white-space:pre-wrap>'                     /note="Alu"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /rpt_family="Alu"'</span></li>
	<li><span style=white-space:pre-wrap>'     misc_feature    632'</span></li>
	<li><span style=white-space:pre-wrap>'                     /note="probable polymorphism (c/t)"'</span></li>
	<li><span style=white-space:pre-wrap>'     gene            747..845'</span></li>
	<li><span style=white-space:pre-wrap>'                     /gene="BRCA1"'</span></li>
	<li><span style=white-space:pre-wrap>'     exon            747..845'</span></li>
	<li><span style=white-space:pre-wrap>'                     /gene="BRCA1"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /number=2'</span></li>
	<li><span style=white-space:pre-wrap>'     CDS             766..&gt;845'</span></li>
	<li><span style=white-space:pre-wrap>'                     /gene="BRCA1"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /codon_start=1'</span></li>
	<li><span style=white-space:pre-wrap>'                     /product="BRCA1"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /protein_id="AAB40910.1"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /translation="MDLSALRVEEVQNVINAMQKILECPI"'</span></li>
	<li><span style=white-space:pre-wrap>'     repeat_region   1201..1338'</span></li>
	<li><span style=white-space:pre-wrap>'                     /note="Alu"'</span></li>
	<li><span style=white-space:pre-wrap>'                     /rpt_family="Alu"'</span></li>
	<li><span style=white-space:pre-wrap>'ORIGIN      '</span></li>
</ol>




```R
q2$req[[1]]
```


          name     length      frame     ncbicg 
    "HSU61268"     "1338"        "0"        "1" 



```R
cat(annots, sep = "\n")
cat(annots, file = "D:/Try-practice/Chapter 3/3.1Retrieving annots.txt", sep = "\n")
# like `print`
```

    LOCUS       HSU61268                1338 bp    DNA     linear   PRI 16-JAN-1997
    DEFINITION  Human breast and ovarian cancer susceptibility (BRCA1) gene, exon
                2, partial flanking introns, and partial cds.
    ACCESSION   U61268
    VERSION     U61268.1
    KEYWORDS    .
    SOURCE      Homo sapiens (human)
      ORGANISM  Homo sapiens
                Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
                Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
                Catarrhini; Hominidae; Homo.
    REFERENCE   1  (bases 1 to 1338)
      AUTHORS   Struewing,J.P., Robbins,C. and Brody,L.C.
      TITLE     BRCA1: genomic sequence flanking exon 2
      JOURNAL   Unpublished
    REFERENCE   2  (bases 1 to 1338)
      AUTHORS   Struewing,J.P., Robbins,C. and Brody,L.C.
      TITLE     Direct Submission
      JOURNAL   Submitted (19-JUN-1996) NCI/NCHGR, NIH, 6130 Executive Blvd MSC
                7372, Bethesda, MD 20892-7372, USA
    FEATURES             Location/Qualifiers
         source          1..1338
                         /organism="Homo sapiens"
                         /mol_type="genomic DNA"
                         /db_xref="taxon:9606"
                         /chromosome="17"
                         /map="17q21"
         repeat_region   1..84
                         /note="Alu"
                         /rpt_family="Alu"
         misc_feature    632
                         /note="probable polymorphism (c/t)"
         gene            747..845
                         /gene="BRCA1"
         exon            747..845
                         /gene="BRCA1"
                         /number=2
         CDS             766..>845
                         /gene="BRCA1"
                         /codon_start=1
                         /product="BRCA1"
                         /protein_id="AAB40910.1"
                         /translation="MDLSALRVEEVQNVINAMQKILECPI"
         repeat_region   1201..1338
                         /note="Alu"
                         /rpt_family="Alu"
    ORIGIN      
    


```R
?cat
```


```R
closebank()
```


```R

```
