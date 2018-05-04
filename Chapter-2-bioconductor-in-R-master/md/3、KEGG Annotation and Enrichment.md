
# KEGG Annotation


```R
source("http://bioconductor.org/biocLite.R") 
biocLite("KEGG.db")
biocLite("KEGGREST")
biocLite("clusterProfiler")  # for enrichment
```

    Bioconductor version 3.6 (BiocInstaller 1.28.0), ?biocLite for help
    


```R
library(KEGG.db) # loads the  library
myEIDs <- c("1109", "6718") # Create vecotor of input Entrez IDs
kegg <- as.character(unlist(mget(as.character(myEIDs), KEGGEXTID2PATHID, ifnotfound=NA)))
kegg <- sapply(strsplit(kegg, "hsa"), "[[", 2)
myPath <- unlist(mget(kegg, KEGGPATHID2NAME, ifnotfound=list(NA)))
myPath
KEGGPATHID2EXTID$hsa00140
KEGGPATHID2EXTID$sce00100
```


<dl class=dl-horizontal>
	<dt>00120</dt>
		<dd>'Primary bile acid biosynthesis'</dd>
	<dt>00140</dt>
		<dd>'Steroid hormone biosynthesis'</dd>
	<dt>00980</dt>
		<dd>'Metabolism of xenobiotics by cytochrome P450'</dd>
	<dt>01100</dt>
		<dd>'Metabolic pathways'</dd>
	<dt>00120</dt>
		<dd>'Primary bile acid biosynthesis'</dd>
	<dt>00140</dt>
		<dd>'Steroid hormone biosynthesis'</dd>
	<dt>01100</dt>
		<dd>'Metabolic pathways'</dd>
</dl>




<ol class=list-inline>
	<li>'100510686'</li>
	<li>'10720'</li>
	<li>'10941'</li>
	<li>'1109'</li>
	<li>'1312'</li>
	<li>'1543'</li>
	<li>'1545'</li>
	<li>'1551'</li>
	<li>'1576'</li>
	<li>'1577'</li>
	<li>'1581'</li>
	<li>'1583'</li>
	<li>'1584'</li>
	<li>'1585'</li>
	<li>'1586'</li>
	<li>'1588'</li>
	<li>'1589'</li>
	<li>'1645'</li>
	<li>'1646'</li>
	<li>'3283'</li>
	<li>'3284'</li>
	<li>'3290'</li>
	<li>'3291'</li>
	<li>'3292'</li>
	<li>'3293'</li>
	<li>'3294'</li>
	<li>'412'</li>
	<li>'51144'</li>
	<li>'51478'</li>
	<li>'54490'</li>
	<li>'54575'</li>
	<li>'54576'</li>
	<li>'54577'</li>
	<li>'54578'</li>
	<li>'54579'</li>
	<li>'54600'</li>
	<li>'54657'</li>
	<li>'54658'</li>
	<li>'54659'</li>
	<li>'574537'</li>
	<li>'64816'</li>
	<li>'6715'</li>
	<li>'6716'</li>
	<li>'6718'</li>
	<li>'6783'</li>
	<li>'6820'</li>
	<li>'7363'</li>
	<li>'7364'</li>
	<li>'7365'</li>
	<li>'7366'</li>
	<li>'7367'</li>
	<li>'7923'</li>
	<li>'79644'</li>
	<li>'79799'</li>
	<li>'8630'</li>
	<li>'8644'</li>
	<li>'9420'</li>
</ol>




<ol class=list-inline>
	<li>'YCR048W'</li>
	<li>'YGL001C'</li>
	<li>'YGL012W'</li>
	<li>'YGR060W'</li>
	<li>'YGR175C'</li>
	<li>'YHR007C'</li>
	<li>'YHR072W'</li>
	<li>'YHR190W'</li>
	<li>'YLR056W'</li>
	<li>'YLR100W'</li>
	<li>'YML008C'</li>
	<li>'YMR015C'</li>
	<li>'YMR202W'</li>
	<li>'YNL280C'</li>
	<li>'YNR019W'</li>
</ol>



The `KEGGREST` library can query the KEGG pathway database via the KEGG REST server for all the genes, enzymes, compounds, and reactions that are involved in the interactions in the desired pathway.


```R
library(KEGGREST)
genes <- keggGet("hsa00140")
genes
```


<ol>
	<li><dl>
	<dt>$ENTRY</dt>
		<dd><strong>Pathway:</strong> 'hsa00140'</dd>
	<dt>$NAME</dt>
		<dd>'Steroid hormone biosynthesis - Homo sapiens (human)'</dd>
	<dt>$DESCRIPTION</dt>
		<dd>'Steroid hormones derived from cholesterol are a class of biologically active compounds in vertebrates. The cholesterol side-chain cleavage enzyme CYP11A1 catalyzes conversion of cholesterol, a C27 compound, to the first C21 steroid, pregnenolone, which is converted by a bifunctional enzyme complex to the gestagen hormone, progesterone [MD:M00107]. Pregnenolone and progesterone are the starting materials for the three groups of steroids: C21 steroids of glucocorticoids and mineralocorticoids, C19 steroids of androgens, and C18 steroids of estrogens. (i) Progesterone is converted by hydroxylations at carbons 21 and 11 to corticosterone, which is further modified by hydroxylation and oxydoreduction at carbon 18 to yield aldosterone, a mineralcorticoid [MD:M00108]. Cortisol, the main glucocorticoid, is formed from 17alpha-hydroxyprogesterone with 11-deoxycortisol as an intermediate [MD:M00109]. (ii) Male hormone testosterone is formed from pregnenolone by two pathways, delta5 pathway via dehydroepiandrosterone and delta4 pathway via androstenedione [MD:M00110]. The enzyme CYP17A1 is responsible for the 17,20 lyase and 17alpha-hydroxylase activities in respective pathways. (iii) Female hormones estrone and estradiol are formed from testosterone and 4-androstene-3,17-dione by oxidative removal of the C19 methyl group and subsequent aromatization of ring A [MD:M00111]. In addition to these three groups, recent studies show that there is another group, termed neurosteroids, synthesized in the brain rather than the peripheral endocrine gland.'</dd>
	<dt>$CLASS</dt>
		<dd>'Metabolism; Lipid metabolism'</dd>
	<dt>$PATHWAY_MAP</dt>
		<dd><strong>hsa00140:</strong> 'Steroid hormone biosynthesis'</dd>
	<dt>$MODULE</dt>
		<dd><dl class=dl-horizontal>
	<dt>hsa_M00107</dt>
		<dd>'Steroid hormone biosynthesis, cholesterol =&gt; prognenolone =&gt; progesterone [PATH:hsa00140]'</dd>
	<dt>hsa_M00108</dt>
		<dd>'C21-Steroid hormone biosynthesis, progesterone =&gt; corticosterone/aldosterone [PATH:hsa00140]'</dd>
	<dt>hsa_M00109</dt>
		<dd>'C21-Steroid hormone biosynthesis, progesterone =&gt; cortisol/cortisone [PATH:hsa00140]'</dd>
	<dt>hsa_M00110</dt>
		<dd>'C19/C18-Steroid hormone biosynthesis, pregnenolone =&gt; androstenedione =&gt; estrone [PATH:hsa00140]'</dd>
</dl>
</dd>
	<dt>$DISEASE</dt>
		<dd><dl class=dl-horizontal>
	<dt>H00134</dt>
		<dd>'X-linked ichthyosis (XLI)'</dd>
	<dt>H00216</dt>
		<dd>'Congenital adrenal hyperplasia (CAH)'</dd>
	<dt>H00258</dt>
		<dd>'Aldosterone synthase deficiency'</dd>
	<dt>H00259</dt>
		<dd>'Apparent mineralocorticoid excess syndrome'</dd>
	<dt>H00599</dt>
		<dd>'46,XX disorders of sex development (Disorders related to androgen excess)'</dd>
	<dt>H00602</dt>
		<dd>'Glucocorticoid-remediable aldosteronism (GRA)'</dd>
	<dt>H00608</dt>
		<dd>'46,XY disorders of sex development (Disorders in androgen synthesis or action)'</dd>
	<dt>H00628</dt>
		<dd>'Congenital bile acid synthesis defect (CBAS)'</dd>
	<dt>H00794</dt>
		<dd>'Aromatase excess syndrome'</dd>
	<dt>H01075</dt>
		<dd>'Peters anomaly'</dd>
	<dt>H01111</dt>
		<dd>'Cortisone reductase deficiency (CRD)'</dd>
	<dt>H01203</dt>
		<dd>'Primary congenital glaucoma (PCG)'</dd>
	<dt>H01709</dt>
		<dd>'Glucocorticoid-induced osteonecrosis'</dd>
	<dt>H02020</dt>
		<dd>'Aromatase deficiency'</dd>
</dl>
</dd>
	<dt>$DRUG</dt>
		<dd><ol class=list-inline>
	<li>'D00152'</li>
	<li>'Phloroglucinol (JAN)'</li>
	<li>'D00153'</li>
	<li>'Testolactone (USP/INN)'</li>
	<li>'D00156'</li>
	<li>'Glycyrrhetinic acid (JAN)'</li>
	<li>'D00321'</li>
	<li>'Finasteride (JAN/USP/INN)'</li>
	<li>'D00351'</li>
	<li>'Ketoconazole (JP17/USP)'</li>
	<li>'D00410'</li>
	<li>'Metyrapone (JP17/USP/INN)'</li>
	<li>'D00420'</li>
	<li>'Mitotane (JAN/USP/INN)'</li>
	<li>'D00574'</li>
	<li>'Aminoglutethimide (USP/INN)'</li>
	<li>'D00960'</li>
	<li>'Anastrozole (JAN/USAN/INN)'</li>
	<li>'D00963'</li>
	<li>'Exemestane (JAN/USP/INN)'</li>
	<li>'D00964'</li>
	<li>'Letrozole (JAN/USP/INN)'</li>
	<li>'D01134'</li>
	<li>'Epristeride (JAN/USAN/INN)'</li>
	<li>'D01180'</li>
	<li>'Trilostane (JAN/USAN)'</li>
	<li>'D01259'</li>
	<li>'Flopropione (JP17/INN)'</li>
	<li>'D01899'</li>
	<li>'Carbenoxolone sodium (JAN/USAN)'</li>
	<li>'D02451'</li>
	<li>'Fadrozole hydrochloride (USAN)'</li>
	<li>'D03107'</li>
	<li>'Bexlosteride (USAN/INN)'</li>
	<li>'D03749'</li>
	<li>'Plomestane (USAN/INN)'</li>
	<li>'D03778'</li>
	<li>'Fadrozole hydrochloride hydrate (JAN)'</li>
	<li>'D03781'</li>
	<li>'Liarozole fumarate (USAN)'</li>
	<li>'D03784'</li>
	<li>'Liarozole hydrochloride (USAN)'</li>
	<li>'D03786'</li>
	<li>'Vorozole (USAN/INN)'</li>
	<li>'D03820'</li>
	<li>'Dutasteride (JAN/USAN/INN)'</li>
	<li>'D04035'</li>
	<li>'Epostane (USAN/INN)'</li>
	<li>'D04498'</li>
	<li>'Idronoxil (USAN/INN)'</li>
	<li>'D04646'</li>
	<li>'Izonsteride (USAN/INN)'</li>
	<li>'D05019'</li>
	<li>'Metyrapone tartrate (USAN)'</li>
	<li>'D07260'</li>
	<li>'Formestane (INN)'</li>
	<li>'D07615'</li>
	<li>'Carbenoxolone (INN)'</li>
	<li>'D07940'</li>
	<li>'Fadrozole (INN)'</li>
	<li>'D08369'</li>
	<li>'Polyphloroglucinol phosphate'</li>
	<li>'D09701'</li>
	<li>'Abiraterone acetate (JAN/USAN)'</li>
	<li>'D09915'</li>
	<li>'Irosustat (USAN/INN)'</li>
	<li>'D10125'</li>
	<li>'Galeterone (USAN)'</li>
	<li>'D10146'</li>
	<li>'Orteronel (JAN/USAN)'</li>
	<li>'D11061'</li>
	<li>'Osilodrostat (USAN/INN)'</li>
</ol>
</dd>
	<dt>$DBLINKS</dt>
		<dd><ol class=list-inline>
	<li>'BSID: 82940'</li>
	<li>'GO: 0034754'</li>
</ol>
</dd>
	<dt>$ORGANISM</dt>
		<dd><strong>Homo sapiens (human) [GN:hsa]:</strong> <span style=white-space:pre-wrap>'NA  Homo sapiens (human) [GN:hsa]'</span></dd>
	<dt>$GENE</dt>
		<dd><ol class=list-inline>
	<li>'1583'</li>
	<li>'CYP11A1; cytochrome P450 family 11 subfamily A member 1 [KO:K00498] [EC:1.14.15.6]'</li>
	<li>'1586'</li>
	<li>'CYP17A1; cytochrome P450 family 17 subfamily A member 1 [KO:K00512] [EC:1.14.14.32 1.14.14.19]'</li>
	<li>'412'</li>
	<li>'STS; steroid sulfatase [KO:K01131] [EC:3.1.6.2]'</li>
	<li>'6820'</li>
	<li>'SULT2B1; sulfotransferase family 2B member 1 [KO:K01015] [EC:2.8.2.2]'</li>
	<li>'1589'</li>
	<li>'CYP21A2; cytochrome P450 family 21 subfamily A member 2 [KO:K00513] [EC:1.14.14.16]'</li>
	<li>'3283'</li>
	<li>'HSD3B1; hydroxy-delta-5-steroid dehydrogenase, 3 beta- and steroid delta-isomerase 1 [KO:K00070] [EC:5.3.3.1 1.1.1.145]'</li>
	<li>'3284'</li>
	<li>'HSD3B2; hydroxy-delta-5-steroid dehydrogenase, 3 beta- and steroid delta-isomerase 2 [KO:K00070] [EC:5.3.3.1 1.1.1.145]'</li>
	<li>'6715'</li>
	<li>'SRD5A1; steroid 5 alpha-reductase 1 [KO:K12343] [EC:1.3.1.22]'</li>
	<li>'6716'</li>
	<li>'SRD5A2; steroid 5 alpha-reductase 2 [KO:K12344] [EC:1.3.1.22]'</li>
	<li>'79644'</li>
	<li>'SRD5A3; steroid 5 alpha-reductase 3 [KO:K12345] [EC:1.3.1.94 1.3.1.22]'</li>
	<li>'1646'</li>
	<li>'AKR1C2; aldo-keto reductase family 1 member C2 [KO:K00089] [EC:1.1.1.357 1.1.1.213]'</li>
	<li>'8644'</li>
	<li>'AKR1C3; aldo-keto reductase family 1 member C3 [KO:K04119] [EC:1.1.1.357 1.1.1.213 1.1.1.188 1.1.1.51]'</li>
	<li>'1584'</li>
	<li>'CYP11B1; cytochrome P450 family 11 subfamily B member 1 [KO:K00497] [EC:1.14.15.4]'</li>
	<li>'1585'</li>
	<li>'CYP11B2; cytochrome P450 family 11 subfamily B member 2 [KO:K07433] [EC:1.14.15.5 1.14.15.4]'</li>
	<li>'6718'</li>
	<li>'AKR1D1; aldo-keto reductase family 1 member D1 [KO:K00251] [EC:1.3.1.3]'</li>
	<li>'1109'</li>
	<li>'AKR1C4; aldo-keto reductase family 1 member C4 [KO:K00037] [EC:1.1.1.225 1.1.1.357 1.1.1.50]'</li>
	<li>'3290'</li>
	<li>'HSD11B1; hydroxysteroid 11-beta dehydrogenase 1 [KO:K15680] [EC:1.1.1.146]'</li>
	<li>'3291'</li>
	<li>'HSD11B2; hydroxysteroid 11-beta dehydrogenase 2 [KO:K00071] [EC:1.1.1.-]'</li>
	<li>'1645'</li>
	<li>'AKR1C1; aldo-keto reductase family 1 member C1 [KO:K00212] [EC:1.3.1.20 1.1.1.357 1.1.1.149]'</li>
	<li>'9420'</li>
	<li>'CYP7B1; cytochrome P450 family 7 subfamily B member 1 [KO:K07430] [EC:1.14.14.29]'</li>
	<li>'6783'</li>
	<li>'SULT1E1; sulfotransferase family 1E member 1 [KO:K01016] [EC:2.8.2.4]'</li>
	<li>'3292'</li>
	<li>'HSD17B1; hydroxysteroid 17-beta dehydrogenase 1 [KO:K00044] [EC:1.1.1.62]'</li>
	<li>'3294'</li>
	<li>'HSD17B2; hydroxysteroid 17-beta dehydrogenase 2 [KO:K13368] [EC:1.1.1.239 1.1.1.62]'</li>
	<li>'8630'</li>
	<li>'HSD17B6; hydroxysteroid 17-beta dehydrogenase 6 [KO:K13369] [EC:1.1.1.239 1.1.1.105 1.1.1.62]'</li>
	<li>'51478'</li>
	<li>'HSD17B7; hydroxysteroid 17-beta dehydrogenase 7 [KO:K13373] [EC:1.1.1.270 1.1.1.62]'</li>
	<li>'7923'</li>
	<li>'HSD17B8; hydroxysteroid 17-beta dehydrogenase 8 [KO:K13370] [EC:1.1.1.239 1.1.1.62]'</li>
	<li>'51144'</li>
	<li>'HSD17B12; hydroxysteroid 17-beta dehydrogenase 12 [KO:K10251] [EC:1.1.1.330 1.1.1.62]'</li>
	<li>'1543'</li>
	<li>'CYP1A1; cytochrome P450 family 1 subfamily A member 1 [KO:K07408] [EC:1.14.14.1]'</li>
	<li>'1544'</li>
	<li>'CYP1A2; cytochrome P450 family 1 subfamily A member 2 [KO:K07409] [EC:1.14.14.1]'</li>
	<li>'1577'</li>
	<li>'CYP3A5; cytochrome P450 family 3 subfamily A member 5 [KO:K17690] [EC:1.14.14.1]'</li>
	<li>'1551'</li>
	<li>'CYP3A7; cytochrome P450 family 3 subfamily A member 7 [KO:K17691] [EC:1.14.14.1]'</li>
	<li>'100861540'</li>
	<li>'CYP3A7-CYP3A51P; CYP3A7-CYP3A51P readthrough [KO:K17691] [EC:1.14.14.1]'</li>
	<li>'1571'</li>
	<li>'CYP2E1; cytochrome P450 family 2 subfamily E member 1 [KO:K07415] [EC:1.14.13.-]'</li>
	<li>'1576'</li>
	<li>'CYP3A4; cytochrome P450 family 3 subfamily A member 4 [KO:K17689] [EC:1.14.13.- 1.14.14.57 1.14.14.56 1.14.14.55 1.14.13.32]'</li>
	<li>'1545'</li>
	<li>'CYP1B1; cytochrome P450 family 1 subfamily B member 1 [KO:K07410] [EC:1.14.14.1]'</li>
	<li>'1588'</li>
	<li>'CYP19A1; cytochrome P450 family 19 subfamily A member 1 [KO:K07434] [EC:1.14.14.14]'</li>
	<li>'1581'</li>
	<li>'CYP7A1; cytochrome P450 family 7 subfamily A member 1 [KO:K00489] [EC:1.14.14.23]'</li>
	<li>'10941'</li>
	<li>'UGT2A1; UDP glucuronosyltransferase family 2 member A1 complex locus [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'79799'</li>
	<li>'UGT2A3; UDP glucuronosyltransferase family 2 member A3 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'7367'</li>
	<li>'UGT2B17; UDP glucuronosyltransferase family 2 member B17 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'10720'</li>
	<li>'UGT2B11; UDP glucuronosyltransferase family 2 member B11 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54490'</li>
	<li>'UGT2B28; UDP glucuronosyltransferase family 2 member B28 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54578'</li>
	<li>'UGT1A6; UDP glucuronosyltransferase family 1 member A6 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54657'</li>
	<li>'UGT1A4; UDP glucuronosyltransferase family 1 member A4 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54658'</li>
	<li>'UGT1A1; UDP glucuronosyltransferase family 1 member A1 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54659'</li>
	<li>'UGT1A3; UDP glucuronosyltransferase family 1 member A3 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'7365'</li>
	<li>'UGT2B10; UDP glucuronosyltransferase family 2 member B10 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54600'</li>
	<li>'UGT1A9; UDP glucuronosyltransferase family 1 member A9 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'7364'</li>
	<li>'UGT2B7; UDP glucuronosyltransferase family 2 member B7 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54575'</li>
	<li>'UGT1A10; UDP glucuronosyltransferase family 1 member A10 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54576'</li>
	<li>'UGT1A8; UDP glucuronosyltransferase family 1 member A8 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54579'</li>
	<li>'UGT1A5; UDP glucuronosyltransferase family 1 member A5 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'7366'</li>
	<li>'UGT2B15; UDP glucuronosyltransferase family 2 member B15 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'54577'</li>
	<li>'UGT1A7; UDP glucuronosyltransferase family 1 member A7 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'7363'</li>
	<li>'UGT2B4; UDP glucuronosyltransferase family 2 member B4 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'574537'</li>
	<li>'UGT2A2; UDP glucuronosyltransferase family 2 member A2 [KO:K00699] [EC:2.4.1.17]'</li>
	<li>'1312'</li>
	<li>'COMT; catechol-O-methyltransferase [KO:K00545] [EC:2.1.1.6]'</li>
	<li>'220074'</li>
	<li>'LRTOMT; leucine rich transmembrane and O-methyltransferase domain containing [KO:K00545] [EC:2.1.1.6]'</li>
	<li>'3293'</li>
	<li>'HSD17B3; hydroxysteroid 17-beta dehydrogenase 3 [KO:K10207] [EC:1.1.1.64]'</li>
</ol>
</dd>
	<dt>$COMPOUND</dt>
		<dd><dl class=dl-horizontal>
	<dt>C00187</dt>
		<dd>'Cholesterol'</dd>
	<dt>C00280</dt>
		<dd>'Androstenedione'</dd>
	<dt>C00410</dt>
		<dd>'Progesterone'</dd>
	<dt>C00468</dt>
		<dd>'Estrone'</dd>
	<dt>C00523</dt>
		<dd>'Androsterone'</dd>
	<dt>C00535</dt>
		<dd>'Testosterone'</dd>
	<dt>C00674</dt>
		<dd>'5alpha-Androstane-3,17-dione'</dd>
	<dt>C00735</dt>
		<dd>'Cortisol'</dd>
	<dt>C00762</dt>
		<dd>'Cortisone'</dd>
	<dt>C00951</dt>
		<dd>'Estradiol-17beta'</dd>
	<dt>C01124</dt>
		<dd>'18-Hydroxycorticosterone'</dd>
	<dt>C01176</dt>
		<dd>'17alpha-Hydroxyprogesterone'</dd>
	<dt>C01227</dt>
		<dd>'Dehydroepiandrosterone'</dd>
	<dt>C01780</dt>
		<dd>'Aldosterone'</dd>
	<dt>C01953</dt>
		<dd>'Pregnenolone'</dd>
	<dt>C02140</dt>
		<dd>'Corticosterone'</dd>
	<dt>C02373</dt>
		<dd>'4-Methylpentanal'</dd>
	<dt>C02537</dt>
		<dd>'Estradiol-17alpha'</dd>
	<dt>C02538</dt>
		<dd>'Estrone 3-sulfate'</dd>
	<dt>C03205</dt>
		<dd>'11-Deoxycorticosterone'</dd>
	<dt>C03681</dt>
		<dd>'5alpha-Pregnane-3,20-dione'</dd>
	<dt>C03747</dt>
		<dd>'11alpha-Hydroxyprogesterone'</dd>
	<dt>C03772</dt>
		<dd>'5beta-Androstane-3,17-dione'</dd>
	<dt>C03852</dt>
		<dd>'Androstan-3alpha,17beta-diol'</dd>
	<dt>C03917</dt>
		<dd>'Dihydrotestosterone'</dd>
	<dt>C03935</dt>
		<dd>'6beta-Hydroxy-17beta-estradiol'</dd>
	<dt>C04042</dt>
		<dd>'20alpha-Hydroxy-4-pregnen-3-one'</dd>
	<dt>C04295</dt>
		<dd>'Androstenediol'</dd>
	<dt>C04373</dt>
		<dd>'Etiocholanolone'</dd>
	<dt>C04518</dt>
		<dd>'17alpha,20alpha-Dihydroxypregn-4-en-3-one'</dd>
	<dt>C04555</dt>
		<dd>'Dehydroepiandrosterone sulfate'</dd>
	<dt>C04676</dt>
		<dd>'Testololactone'</dd>
	<dt>C05138</dt>
		<dd>'17alpha-Hydroxypregnenolone'</dd>
	<dt>C05139</dt>
		<dd>'16alpha-Hydroxydehydroepiandrosterone'</dd>
	<dt>C05140</dt>
		<dd>'16alpha-Hydroxyandrost-4-ene-3,17-dione'</dd>
	<dt>C05141</dt>
		<dd>'Estriol'</dd>
	<dt>C05284</dt>
		<dd>'11beta-Hydroxyandrost-4-ene-3,17-dione'</dd>
	<dt>C05285</dt>
		<dd>'Adrenosterone'</dd>
	<dt>C05290</dt>
		<dd>'19-Hydroxyandrost-4-ene-3,17-dione'</dd>
	<dt>C05291</dt>
		<dd>'7alpha-Hydroxytestosterone'</dd>
	<dt>C05293</dt>
		<dd>'5beta-Dihydrotestosterone'</dd>
	<dt>C05294</dt>
		<dd>'19-Hydroxytestosterone'</dd>
	<dt>C05295</dt>
		<dd>'19-Oxotestosterone'</dd>
	<dt>C05296</dt>
		<dd>'7alpha-Hydroxyandrost-4-ene-3,17-dione'</dd>
	<dt>C05297</dt>
		<dd>'19-Oxoandrost-4-ene-3,17-dione'</dd>
	<dt>C05298</dt>
		<dd>'2-Hydroxyestrone'</dd>
	<dt>C05299</dt>
		<dd>'2-Methoxyestrone'</dd>
	<dt>C05300</dt>
		<dd>'16alpha-Hydroxyestrone'</dd>
	<dt>C05301</dt>
		<dd>'2-Hydroxyestradiol'</dd>
	<dt>C05302</dt>
		<dd>'2-Methoxy-17beta-estradiol'</dd>
	<dt>C05469</dt>
		<dd>'17alpha,21-Dihydroxy-5beta-pregnane-3,11,20-trione'</dd>
	<dt>C05470</dt>
		<dd>'Tetrahydrocortisone'</dd>
	<dt>C05471</dt>
		<dd>'11beta,17alpha,21-Trihydroxy-5beta-pregnane-3,20-dione'</dd>
	<dt>C05472</dt>
		<dd>'Urocortisol'</dd>
	<dt>C05473</dt>
		<dd>'11beta,21-Dihydroxy-3,20-oxo-5beta-pregnan-18-al'</dd>
	<dt>C05474</dt>
		<dd>'3alpha,11beta,21-Trihydroxy-20-oxo-5beta-pregnan-18-al'</dd>
	<dt>C05475</dt>
		<dd>'11beta,21-Dihydroxy-5beta-pregnane-3,20-dione'</dd>
	<dt>C05476</dt>
		<dd>'Tetrahydrocorticosterone'</dd>
	<dt>C05477</dt>
		<dd>'21-Hydroxy-5beta-pregnane-3,11,20-trione'</dd>
	<dt>C05478</dt>
		<dd>'3alpha,21-Dihydroxy-5beta-pregnane-11,20-dione'</dd>
	<dt>C05479</dt>
		<dd>'5beta-Pregnane-3,20-dione'</dd>
	<dt>C05480</dt>
		<dd>'Pregnanolone'</dd>
	<dt>C05481</dt>
		<dd>'Cortolone'</dd>
	<dt>C05482</dt>
		<dd>'Cortol'</dd>
	<dt>C05483</dt>
		<dd>'3alpha,20alpha,21-Trihydroxy-5beta-pregnan-11-one'</dd>
	<dt>C05484</dt>
		<dd>'Pregnanediol'</dd>
	<dt>C05485</dt>
		<dd>'21-Hydroxypregnenolone'</dd>
	<dt>C05487</dt>
		<dd>'17alpha,21-Dihydroxypregnenolone'</dd>
	<dt>C05488</dt>
		<dd>'11-Deoxycortisol'</dd>
	<dt>C05489</dt>
		<dd>'11beta,17alpha,21-Trihydroxypregnenolone'</dd>
	<dt>C05490</dt>
		<dd>'11-Dehydrocorticosterone'</dd>
	<dt>C05497</dt>
		<dd>'21-Deoxycortisol'</dd>
	<dt>C05498</dt>
		<dd>'11beta-Hydroxyprogesterone'</dd>
	<dt>C05499</dt>
		<dd>'17alpha,20alpha-Dihydroxycholesterol'</dd>
	<dt>C05500</dt>
		<dd>'20alpha-Hydroxycholesterol'</dd>
	<dt>C05501</dt>
		<dd>'20alpha,22beta-Dihydroxycholesterol'</dd>
	<dt>C05502</dt>
		<dd>'22(R)-Hydroxycholesterol'</dd>
	<dt>C05503</dt>
		<dd>'Estradiol-17beta 3-glucuronide'</dd>
	<dt>C05504</dt>
		<dd>'16-Glucuronide-estriol'</dd>
	<dt>C08357</dt>
		<dd>'Estradiol-17beta 3-sulfate'</dd>
	<dt>C08358</dt>
		<dd>'2-Methoxyestrone 3-sulfate'</dd>
	<dt>C08359</dt>
		<dd>'2-Methoxyestradiol-17beta 3-sulfate'</dd>
	<dt>C11131</dt>
		<dd>'2-Methoxy-estradiol-17beta 3-glucuronide'</dd>
	<dt>C11132</dt>
		<dd>'2-Methoxyestrone 3-glucuronide'</dd>
	<dt>C11133</dt>
		<dd>'Estrone glucuronide'</dd>
	<dt>C11134</dt>
		<dd>'Testosterone glucuronide'</dd>
	<dt>C11135</dt>
		<dd>'Androsterone glucuronide'</dd>
	<dt>C11136</dt>
		<dd>'Etiocholan-3alpha-ol-17-one 3-glucuronide'</dd>
	<dt>C13712</dt>
		<dd>'Allopregnanolone'</dd>
	<dt>C13713</dt>
		<dd>'Allotetrahydrodeoxycorticosterone'</dd>
	<dt>C18038</dt>
		<dd>'7alpha-Hydroxypregnenolone'</dd>
	<dt>C18039</dt>
		<dd>'Aldosterone hemiacetal'</dd>
	<dt>C18040</dt>
		<dd>'5alpha-Dihydrodeoxycorticosterone'</dd>
	<dt>C18041</dt>
		<dd>'5alpha-Pregnan-20alpha-ol-3-one'</dd>
	<dt>C18042</dt>
		<dd>'5alpha-Pregnane-3alpha,20alpha-diol'</dd>
	<dt>C18043</dt>
		<dd>'Cholesterol sulfate'</dd>
	<dt>C18044</dt>
		<dd>'3beta-Hydroxypregn-5-en-20-one sulfate'</dd>
	<dt>C18045</dt>
		<dd>'7alpha-Hydroxydehydroepiandrosterone'</dd>
	<dt>C18075</dt>
		<dd>'11beta,17beta-Dihydroxy-4-androsten-3-one'</dd>
</dl>
</dd>
	<dt>$KO_PATHWAY</dt>
		<dd>'ko00140'</dd>
	<dt>$REFERENCE</dt>
		<dd><ol>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:10998348'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Penning TM, Burczynski ME, Jez JM, Hung CF, Lin HK, Ma H, Moore M, Palackal N, Ratnam K'</dd>
	<dt>$TITLE</dt>
		<dd>'Human 3alpha-hydroxysteroid dehydrogenase isoforms (AKR1C1-AKR1C4) of the aldo-keto reductase superfamily: functional plasticity and tissue distribution reveals roles in the inactivation and formation of male and female sex hormones.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'Biochem J 351:67-77 (2000)'</li>
	<li>'DOI:10.1042/bj3510067'</li>
</ol>
</dd>
</dl>
</li>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:11469811'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Tomlinson JW, Stewart PM'</dd>
	<dt>$TITLE</dt>
		<dd>'Cortisol metabolism and the role of 11beta-hydroxysteroid dehydrogenase.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'Best Pract Res Clin Endocrinol Metab 15:61-78 (2001)'</li>
	<li>'DOI:10.1053/beem.2000.0119'</li>
</ol>
</dd>
</dl>
</li>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:12604236'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Higaki Y, Usami N, Shintani S, Ishikura S, El-Kabbani O, Hara A'</dd>
	<dt>$TITLE</dt>
		<dd>'Selective and potent inhibitors of human 20alpha-hydroxysteroid dehydrogenase (AKR1C1) that metabolizes neurosteroids derived from progesterone.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'Chem Biol Interact 143-144:503-13 (2003)'</li>
	<li>'DOI:10.1016/S0009-2797(02)00206-5'</li>
</ol>
</dd>
</dl>
</li>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:12832414'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Thomas JL, Duax WL, Addlagatta A, Brandt S, Fuller RR, Norris W'</dd>
	<dt>$TITLE</dt>
		<dd>'Structure/function relationships responsible for coenzyme specificity and the isomerase activity of human type 1 3 beta-hydroxysteroid dehydrogenase/isomerase.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'J Biol Chem 278:35483-90 (2003)'</li>
	<li>'DOI:10.1074/jbc.M304752200'</li>
</ol>
</dd>
</dl>
</li>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:15642792'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Cavigelli SA, Monfort SL, Whitney TK, Mechref YS, Novotny M, McClintock MK'</dd>
	<dt>$TITLE</dt>
		<dd>'Frequent serial fecal corticoid measures from rats reflect circadian and ovarian corticosterone rhythms.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'J Endocrinol 184:153-63 (2005)'</li>
	<li>'DOI:10.1677/joe.1.05935'</li>
</ol>
</dd>
</dl>
</li>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:16413106'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Holmes MC, Seckl JR'</dd>
	<dt>$TITLE</dt>
		<dd>'The role of 11beta-hydroxysteroid dehydrogenases in the brain.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'Mol Cell Endocrinol 248:9-14 (2006)'</li>
	<li>'DOI:10.1016/j.mce.2005.12.002'</li>
</ol>
</dd>
</dl>
</li>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:16547389'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Matsunaga T, Shintani S, Hara A'</dd>
	<dt>$TITLE</dt>
		<dd>'Multiplicity of mammalian reductases for xenobiotic carbonyl compounds.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'Drug Metab Pharmacokinet 21:1-18 (2006)'</li>
	<li>'DOI:10.2133/dmpk.21.1'</li>
</ol>
</dd>
</dl>
</li>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:17459698'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Morris DJ, Latif SA, Hardy MP, Brem AS'</dd>
	<dt>$TITLE</dt>
		<dd>'Endogenous inhibitors (GALFs) of 11beta-hydroxysteroid dehydrogenase isoforms 1 and 2: derivatives of adrenally produced corticosterone and cortisol.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'J Steroid Biochem Mol Biol 104:161-8 (2007)'</li>
	<li>'DOI:10.1016/j.jsbmb.2007.03.020'</li>
</ol>
</dd>
</dl>
</li>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:17475203'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Sanai M, Endo S, Matsunaga T, Ishikura S, Tajima K, El-Kabbani O, Hara A'</dd>
	<dt>$TITLE</dt>
		<dd>'Rat NAD+-dependent 3alpha-hydroxysteroid dehydrogenase (AKR1C17): a member of the aldo-keto reductase family highly expressed in kidney cytosol.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'Arch Biochem Biophys 464:122-9 (2007)'</li>
	<li>'DOI:10.1016/j.abb.2007.04.003'</li>
</ol>
</dd>
</dl>
</li>
	<li><dl>
	<dt>$REFERENCE</dt>
		<dd>'PMID:17926129'</dd>
	<dt>$AUTHORS</dt>
		<dd>'Ghayee HK, Auchus RJ'</dd>
	<dt>$TITLE</dt>
		<dd>'Basic concepts and recent developments in human steroid hormone biosynthesis.'</dd>
	<dt>$JOURNAL</dt>
		<dd><ol class=list-inline>
	<li>'Rev Endocr Metab Disord 8:289-300 (2007)'</li>
	<li>'DOI:10.1007/s11154-007-9052-2'</li>
</ol>
</dd>
</dl>
</li>
</ol>
</dd>
</dl>
</li>
</ol>



# Enrichment


```R
library(clusterProfiler)
data(gcSample)
genes <-gcSample[[4]] # a list of five sets
genes
kegg_enrichment <- enrichKEGG(genes, pvalueCutoff = 0.01)
kegg_enrichment
as.data.frame(kegg_enrichment)
```


<ol class=list-inline>
	<li>'5573'</li>
	<li>'7453'</li>
	<li>'5245'</li>
	<li>'23450'</li>
	<li>'6500'</li>
	<li>'4926'</li>
	<li>'6427'</li>
	<li>'813'</li>
	<li>'10960'</li>
	<li>'5048'</li>
	<li>'3920'</li>
	<li>'2950'</li>
	<li>'5351'</li>
	<li>'9611'</li>
	<li>'55905'</li>
	<li>'473'</li>
	<li>'6749'</li>
	<li>'10606'</li>
	<li>'5214'</li>
	<li>'2804'</li>
	<li>'5701'</li>
	<li>'7057'</li>
	<li>'1434'</li>
	<li>'5425'</li>
	<li>'1363'</li>
	<li>'5226'</li>
	<li>'9867'</li>
	<li>'7078'</li>
	<li>'8553'</li>
	<li>'10054'</li>
	<li>'262'</li>
	<li>'5641'</li>
	<li>'6513'</li>
	<li>'5704'</li>
	<li>'7296'</li>
	<li>'231'</li>
	<li>'11047'</li>
	<li>'327'</li>
	<li>'7153'</li>
	<li>'55233'</li>
	<li>'10969'</li>
	<li>'2012'</li>
	<li>'908'</li>
	<li>'6631'</li>
	<li>'10181'</li>
	<li>'1809'</li>
	<li>'7072'</li>
	<li>'9520'</li>
	<li>'9184'</li>
	<li>'1736'</li>
	<li>'5955'</li>
	<li>'1075'</li>
	<li>'10105'</li>
	<li>'7874'</li>
	<li>'3912'</li>
	<li>'7334'</li>
	<li>'5412'</li>
	<li>'56681'</li>
	<li>'4172'</li>
	<li>'56255'</li>
	<li>'4681'</li>
	<li>'51780'</li>
	<li>'3339'</li>
	<li>'2181'</li>
	<li>'10051'</li>
	<li>'4860'</li>
	<li>'1786'</li>
	<li>'5713'</li>
	<li>'7283'</li>
	<li>'2037'</li>
	<li>'9929'</li>
	<li>'120'</li>
	<li>'9813'</li>
	<li>'1717'</li>
	<li>'3930'</li>
	<li>'5571'</li>
	<li>'9779'</li>
	<li>'79888'</li>
	<li>'9604'</li>
	<li>'51097'</li>
	<li>'23429'</li>
	<li>'994'</li>
	<li>'23300'</li>
	<li>'5552'</li>
	<li>'6059'</li>
	<li>'6241'</li>
	<li>'7319'</li>
	<li>'10217'</li>
	<li>'8502'</li>
	<li>'4175'</li>
	<li>'10769'</li>
	<li>'214'</li>
	<li>'23077'</li>
	<li>'4678'</li>
	<li>'4651'</li>
	<li>'6251'</li>
	<li>'1389'</li>
	<li>'2131'</li>
	<li>'23013'</li>
	<li>'6768'</li>
	<li>'6422'</li>
	<li>'6611'</li>
	<li>'9202'</li>
	<li>'224'</li>
	<li>'3836'</li>
	<li>'10133'</li>
	<li>'22906'</li>
	<li>'5899'</li>
	<li>'3476'</li>
	<li>'4171'</li>
	<li>'66008'</li>
	<li>'25937'</li>
	<li>'4061'</li>
	<li>'10659'</li>
	<li>'1387'</li>
	<li>'9337'</li>
	<li>'7716'</li>
	<li>'5590'</li>
	<li>'9961'</li>
	<li>'8985'</li>
	<li>'8522'</li>
	<li>'6732'</li>
	<li>'10123'</li>
	<li>'10902'</li>
	<li>'143'</li>
	<li>'23564'</li>
	<li>'2633'</li>
	<li>'5159'</li>
	<li>'10274'</li>
	<li>'22862'</li>
	<li>'6890'</li>
	<li>'6720'</li>
	<li>'613'</li>
	<li>'1445'</li>
	<li>'7374'</li>
	<li>'7320'</li>
	<li>'5718'</li>
	<li>'6695'</li>
	<li>'84747'</li>
	<li>'4820'</li>
	<li>'9665'</li>
	<li>'10204'</li>
	<li>'5118'</li>
	<li>'11044'</li>
	<li>'11052'</li>
	<li>'10982'</li>
	<li>'7127'</li>
	<li>'9474'</li>
	<li>'3659'</li>
	<li>'8874'</li>
	<li>'25949'</li>
	<li>'6840'</li>
	<li>'7298'</li>
	<li>'8204'</li>
	<li>'9874'</li>
	<li>'55837'</li>
	<li>'10513'</li>
	<li>'5440'</li>
	<li>'7128'</li>
	<li>'5699'</li>
	<li>'5604'</li>
	<li>'8566'</li>
	<li>'7375'</li>
	<li>'558'</li>
	<li>'6632'</li>
	<li>'11051'</li>
	<li>'7372'</li>
	<li>'5430'</li>
	<li>'95'</li>
	<li>'2634'</li>
	<li>'23428'</li>
	<li>'23518'</li>
	<li>'8625'</li>
	<li>'901'</li>
	<li>'3155'</li>
	<li>'27338'</li>
	<li>'22908'</li>
	<li>'8192'</li>
	<li>'3689'</li>
	<li>'10617'</li>
	<li>'10614'</li>
	<li>'9648'</li>
	<li>'5106'</li>
	<li>'2184'</li>
	<li>'22918'</li>
	<li>'8697'</li>
	<li>'1520'</li>
	<li>'7466'</li>
	<li>'6279'</li>
	<li>'7525'</li>
	<li>'713'</li>
	<li>'11065'</li>
	<li>'23303'</li>
	<li>'2192'</li>
	<li>'3632'</li>
	<li>'10535'</li>
	<li>'8260'</li>
	<li>'10401'</li>
	<li>'857'</li>
	<li>'22796'</li>
	<li>'9693'</li>
	<li>'3800'</li>
	<li>'5156'</li>
	<li>'9674'</li>
	<li>'5819'</li>
	<li>'11215'</li>
	<li>'11103'</li>
	<li>'9650'</li>
	<li>'353'</li>
	<li>'7088'</li>
	<li>'6302'</li>
	<li>'5050'</li>
	<li>'6310'</li>
	<li>'10641'</li>
	<li>'79096'</li>
	<li>'6416'</li>
	<li>'2162'</li>
	<li>'9823'</li>
	<li>'8624'</li>
	<li>'5493'</li>
	<li>'9743'</li>
	<li>'8614'</li>
	<li>'64601'</li>
	<li>'10788'</li>
	<li>'10651'</li>
	<li>'4016'</li>
	<li>'8634'</li>
	<li>'8603'</li>
	<li>'9724'</li>
	<li>'6817'</li>
	<li>'56902'</li>
	<li>'51704'</li>
	<li>'9332'</li>
	<li>'2230'</li>
	<li>'8722'</li>
	<li>'10206'</li>
	<li>'6387'</li>
	<li>'5577'</li>
	<li>'596'</li>
	<li>'4350'</li>
	<li>'6239'</li>
	<li>'8324'</li>
	<li>'10127'</li>
	<li>'3708'</li>
	<li>'9933'</li>
	<li>'2067'</li>
	<li>'3149'</li>
	<li>'6925'</li>
	<li>'1519'</li>
	<li>'9936'</li>
	<li>'51747'</li>
	<li>'22864'</li>
	<li>'6594'</li>
	<li>'4137'</li>
	<li>'6182'</li>
	<li>'4318'</li>
	<li>'9371'</li>
	<li>'9111'</li>
	<li>'2530'</li>
	<li>'5074'</li>
	<li>'9836'</li>
	<li>'5984'</li>
	<li>'11130'</li>
	<li>'4129'</li>
	<li>'9338'</li>
	<li>'9583'</li>
	<li>'8884'</li>
	<li>'5025'</li>
	<li>'902'</li>
	<li>'22795'</li>
	<li>'2791'</li>
	<li>'962'</li>
	<li>'11259'</li>
	<li>'55556'</li>
	<li>'1164'</li>
	<li>'1120'</li>
	<li>'4335'</li>
	<li>'2643'</li>
	<li>'4005'</li>
	<li>'4316'</li>
	<li>'5698'</li>
	<li>'5507'</li>
	<li>'977'</li>
	<li>'6183'</li>
	<li>'6286'</li>
	<li>'7188'</li>
	<li>'25913'</li>
	<li>'341'</li>
	<li>'2581'</li>
	<li>'8321'</li>
	<li>'1909'</li>
	<li>'10040'</li>
	<li>'771'</li>
	<li>'5480'</li>
	<li>'11138'</li>
	<li>'10610'</li>
	<li>'55856'</li>
	<li>'54677'</li>
	<li>'9891'</li>
	<li>'10573'</li>
	<li>'1462'</li>
	<li>'3169'</li>
	<li>'9482'</li>
	<li>'7049'</li>
	<li>'2237'</li>
	<li>'7060'</li>
	<li>'66036'</li>
	<li>'4832'</li>
	<li>'10468'</li>
	<li>'5002'</li>
	<li>'4600'</li>
	<li>'4013'</li>
	<li>'549'</li>
	<li>'4124'</li>
	<li>'7444'</li>
	<li>'6038'</li>
	<li>'5930'</li>
	<li>'2737'</li>
	<li>'3730'</li>
	<li>'9997'</li>
	<li>'8864'</li>
	<li>'3625'</li>
	<li>'5203'</li>
	<li>'8309'</li>
	<li>'1675'</li>
	<li>'26137'</li>
	<li>'3014'</li>
	<li>'10950'</li>
	<li>'10557'</li>
	<li>'80736'</li>
	<li>'7177'</li>
	<li>'1311'</li>
	<li>'638'</li>
	<li>'7286'</li>
	<li>'3790'</li>
	<li>'7884'</li>
	<li>'6627'</li>
	<li>'7644'</li>
	<li>'3159'</li>
	<li>'10625'</li>
	<li>'23061'</li>
	<li>'22977'</li>
	<li>'51230'</li>
	<li>'30008'</li>
	<li>'5143'</li>
	<li>'3223'</li>
	<li>'4782'</li>
	<li>'1466'</li>
	<li>'3189'</li>
	<li>'25875'</li>
	<li>'22878'</li>
	<li>'23499'</li>
	<li>'9736'</li>
	<li>'10195'</li>
	<li>'23524'</li>
	<li>'4616'</li>
	<li>'23047'</li>
	<li>'10370'</li>
	<li>'118'</li>
	<li>'8527'</li>
	<li>'114876'</li>
	<li>'51599'</li>
	<li>'10111'</li>
	<li>'1500'</li>
	<li>'8073'</li>
	<li>'7430'</li>
	<li>'9601'</li>
	<li>'7267'</li>
	<li>'23741'</li>
	<li>'5036'</li>
	<li>'5591'</li>
	<li>'54499'</li>
	<li>'10284'</li>
	<li>'1453'</li>
	<li>'284119'</li>
	<li>'4176'</li>
	<li>'6731'</li>
	<li>'1107'</li>
	<li>'2805'</li>
	<li>'1312'</li>
	<li>'5127'</li>
	<li>'54107'</li>
	<li>'8428'</li>
	<li>'1452'</li>
	<li>'831'</li>
	<li>'5162'</li>
	<li>'4758'</li>
	<li>'8405'</li>
	<li>'5531'</li>
	<li>'3964'</li>
	<li>'22929'</li>
	<li>'8678'</li>
	<li>'23185'</li>
	<li>'3428'</li>
	<li>'57019'</li>
	<li>'57658'</li>
	<li>'65018'</li>
	<li>'10492'</li>
	<li>'10957'</li>
	<li>'10938'</li>
	<li>'9987'</li>
	<li>'8490'</li>
	<li>'80781'</li>
	<li>'7336'</li>
	<li>'8648'</li>
	<li>'10103'</li>
	<li>'123'</li>
	<li>'9100'</li>
	<li>'8611'</li>
	<li>'54870'</li>
	<li>'79143'</li>
	<li>'8660'</li>
	<li>'2353'</li>
	<li>'7852'</li>
	<li>'8543'</li>
	<li>'10979'</li>
	<li>'10436'</li>
	<li>'10061'</li>
	<li>'23438'</li>
	<li>'2186'</li>
	<li>'1410'</li>
	<li>'3400'</li>
	<li>'11321'</li>
	<li>'10767'</li>
	<li>'7355'</li>
	<li>'11168'</li>
	<li>'3551'</li>
	<li>'2874'</li>
	<li>'7508'</li>
	<li>'9470'</li>
	<li>'4200'</li>
	<li>'10418'</li>
	<li>'10553'</li>
	<li>'57819'</li>
	<li>'23291'</li>
	<li>'18'</li>
	<li>'953'</li>
	<li>'25901'</li>
	<li>'54861'</li>
	<li>'10248'</li>
	<li>'11030'</li>
	<li>'23598'</li>
	<li>'5435'</li>
	<li>'1384'</li>
	<li>'6873'</li>
	<li>'8767'</li>
	<li>'26190'</li>
	<li>'9655'</li>
	<li>'7041'</li>
	<li>'8881'</li>
	<li>'1070'</li>
	<li>'11068'</li>
	<li>'10560'</li>
	<li>'7043'</li>
	<li>'1632'</li>
	<li>'55526'</li>
	<li>'2191'</li>
	<li>'6576'</li>
	<li>'4116'</li>
	<li>'2534'</li>
	<li>'1635'</li>
	<li>'158'</li>
	<li>'4646'</li>
	<li>'983'</li>
	<li>'5921'</li>
	<li>'8289'</li>
	<li>'11257'</li>
	<li>'28639'</li>
	<li>'10142'</li>
	<li>'5191'</li>
	<li>'4690'</li>
	<li>'367'</li>
	<li>'27332'</li>
	<li>'1022'</li>
	<li>'8662'</li>
	<li>'2260'</li>
	<li>'4179'</li>
	<li>'7337'</li>
	<li>'64963'</li>
	<li>'3597'</li>
	<li>'4680'</li>
	<li>'58495'</li>
	<li>'1639'</li>
	<li>'9112'</li>
	<li>'220988'</li>
	<li>'53635'</li>
	<li>'10618'</li>
	<li>'5595'</li>
	<li>'9289'</li>
	<li>'90861'</li>
	<li>'26123'</li>
	<li>'25777'</li>
	<li>'23107'</li>
	<li>'23126'</li>
	<li>'11260'</li>
	<li>'57498'</li>
	<li>'7422'</li>
	<li>'25957'</li>
	<li>'23389'</li>
	<li>'9581'</li>
	<li>'23198'</li>
	<li>'57017'</li>
	<li>'4131'</li>
	<li>'5295'</li>
	<li>'7277'</li>
	<li>'667'</li>
	<li>'55568'</li>
	<li>'26058'</li>
	<li>'9444'</li>
	<li>'23063'</li>
	<li>'27346'</li>
	<li>'23048'</li>
	<li>'91754'</li>
	<li>'22998'</li>
	<li>'25940'</li>
	<li>'23028'</li>
	<li>'2618'</li>
	<li>'23091'</li>
	<li>'89910'</li>
	<li>'1997'</li>
	<li>'219654'</li>
	<li>'9778'</li>
	<li>'253782'</li>
	<li>'26128'</li>
	<li>'200734'</li>
	<li>'23522'</li>
	<li>'1289'</li>
	<li>'23030'</li>
	<li>'22982'</li>
	<li>'80308'</li>
	<li>'285527'</li>
	<li>'10486'</li>
	<li>'23258'</li>
	<li>'23295'</li>
	<li>'114882'</li>
	<li>'10000'</li>
	<li>'80205'</li>
	<li>'440026'</li>
	<li>'3842'</li>
	<li>'54505'</li>
	<li>'23338'</li>
	<li>'25976'</li>
	<li>'23177'</li>
	<li>'4763'</li>
	<li>'26472'</li>
	<li>'140890'</li>
	<li>'23210'</li>
	<li>'23336'</li>
	<li>'6935'</li>
	<li>'7090'</li>
	<li>'84162'</li>
	<li>'57037'</li>
	<li>'83700'</li>
	<li>'1955'</li>
	<li>'8760'</li>
	<li>'6655'</li>
	<li>'23405'</li>
	<li>'124565'</li>
	<li>'29'</li>
	<li>'23517'</li>
	<li>'23270'</li>
	<li>'9847'</li>
	<li>'23078'</li>
	<li>'23158'</li>
	<li>'323'</li>
	<li>'22934'</li>
	<li>'23452'</li>
	<li>'23189'</li>
	<li>'9527'</li>
	<li>'4781'</li>
	<li>'513'</li>
	<li>'253959'</li>
	<li>'5576'</li>
	<li>'90993'</li>
	<li>'79882'</li>
	<li>'57493'</li>
	<li>'23234'</li>
	<li>'1287'</li>
	<li>'29015'</li>
	<li>'4238'</li>
	<li>'25903'</li>
	<li>'23067'</li>
	<li>'23247'</li>
	<li>'23361'</li>
	<li>'10424'</li>
	<li>'57613'</li>
	<li>'157680'</li>
	<li>'10963'</li>
	<li>'1028'</li>
	<li>'112479'</li>
	<li>'124152'</li>
	<li>'222161'</li>
	<li>'1992'</li>
	<li>'4494'</li>
	<li>'3831'</li>
	<li>'4141'</li>
	<li>'65258'</li>
	<li>'728498'</li>
	<li>'26156'</li>
	<li>'3643'</li>
	<li>'5912'</li>
	<li>'51491'</li>
	<li>'6421'</li>
	<li>'3693'</li>
	<li>'2066'</li>
	<li>'987'</li>
	<li>'5108'</li>
	<li>'9236'</li>
	<li>'8725'</li>
	<li>'57326'</li>
	<li>'7703'</li>
	<li>'27122'</li>
	<li>'30'</li>
	<li>'4839'</li>
	<li>'8833'</li>
	<li>'5898'</li>
	<li>'5387'</li>
	<li>'7586'</li>
	<li>'10781'</li>
	<li>'891'</li>
	<li>'388677'</li>
	<li>'54432'</li>
	<li>'8801'</li>
	<li>'55660'</li>
	<li>'5228'</li>
	<li>'30968'</li>
	<li>'56271'</li>
	<li>'2747'</li>
	<li>'5688'</li>
	<li>'3913'</li>
	<li>'100423062'</li>
	<li>'1277'</li>
	<li>'166'</li>
	<li>'11180'</li>
	<li>'51155'</li>
	<li>'79811'</li>
	<li>'64755'</li>
	<li>'51133'</li>
	<li>'53373'</li>
	<li>'28977'</li>
	<li>'55245'</li>
	<li>'55914'</li>
	<li>'10055'</li>
	<li>'23469'</li>
	<li>'58478'</li>
	<li>'29105'</li>
	<li>'54927'</li>
	<li>'51292'</li>
	<li>'55718'</li>
	<li>'64834'</li>
	<li>'79137'</li>
	<li>'51203'</li>
	<li>'28998'</li>
	<li>'10328'</li>
	<li>'55082'</li>
	<li>'51647'</li>
	<li>'11124'</li>
	<li>'56654'</li>
	<li>'3267'</li>
	<li>'55608'</li>
	<li>'55858'</li>
	<li>'51071'</li>
	<li>'55148'</li>
	<li>'65993'</li>
	<li>'54884'</li>
	<li>'55177'</li>
	<li>'54815'</li>
	<li>'55696'</li>
	<li>'51185'</li>
	<li>'55893'</li>
	<li>'79581'</li>
	<li>'64847'</li>
	<li>'79139'</li>
	<li>'56995'</li>
	<li>'51447'</li>
	<li>'79624'</li>
	<li>'60314'</li>
	<li>'712'</li>
	<li>'51118'</li>
	<li>'23560'</li>
	<li>'51111'</li>
	<li>'63901'</li>
	<li>'64777'</li>
	<li>'79665'</li>
	<li>'51642'</li>
	<li>'55316'</li>
	<li>'27342'</li>
	<li>'65982'</li>
	<li>'11108'</li>
	<li>'27244'</li>
	<li>'51053'</li>
	<li>'55204'</li>
	<li>'79641'</li>
	<li>'25959'</li>
	<li>'54585'</li>
	<li>'51397'</li>
	<li>'56942'</li>
	<li>'55325'</li>
	<li>'79887'</li>
	<li>'582'</li>
	<li>'79709'</li>
	<li>'56912'</li>
	<li>'50650'</li>
	<li>'55759'</li>
	<li>'51306'</li>
	<li>'152006'</li>
	<li>'54463'</li>
	<li>'55109'</li>
	<li>'55781'</li>
	<li>'64761'</li>
	<li>'51115'</li>
	<li>'55268'</li>
	<li>'29995'</li>
	<li>'79600'</li>
	<li>'10161'</li>
	<li>'27440'</li>
	<li>'55847'</li>
	<li>'54872'</li>
	<li>'10186'</li>
	<li>'8322'</li>
	<li>'23753'</li>
	<li>'55144'</li>
	<li>'64285'</li>
	<li>'55638'</li>
	<li>'51309'</li>
	<li>'29107'</li>
	<li>'54778'</li>
	<li>'53407'</li>
	<li>'29911'</li>
	<li>'51205'</li>
	<li>'10451'</li>
	<li>'55227'</li>
	<li>'10848'</li>
	<li>'27242'</li>
	<li>'64798'</li>
	<li>'51728'</li>
	<li>'79682'</li>
	<li>'79763'</li>
	<li>'64418'</li>
	<li>'64924'</li>
	<li>'79269'</li>
	<li>'57124'</li>
	<li>'55164'</li>
	<li>'54433'</li>
	<li>'54845'</li>
	<li>'79170'</li>
	<li>'56947'</li>
	<li>'51192'</li>
	<li>'65003'</li>
	<li>'64236'</li>
	<li>'55833'</li>
	<li>'51182'</li>
	<li>'54961'</li>
	<li>'80310'</li>
	<li>'79652'</li>
	<li>'51001'</li>
	<li>'10365'</li>
	<li>'51816'</li>
	<li>'55841'</li>
	<li>'29922'</li>
	<li>'93664'</li>
	<li>'757'</li>
	<li>'79818'</li>
	<li>'56935'</li>
	<li>'26502'</li>
	<li>'51340'</li>
	<li>'4054'</li>
	<li>'51022'</li>
	<li>'84914'</li>
	<li>'51478'</li>
	<li>'63933'</li>
	<li>'80127'</li>
	<li>'51078'</li>
	<li>'51287'</li>
	<li>'55764'</li>
	<li>'79939'</li>
	<li>'78996'</li>
	<li>'112398'</li>
	<li>'81855'</li>
	<li>'81603'</li>
	<li>'81035'</li>
	<li>'81034'</li>
	<li>'79961'</li>
	<li>'1154'</li>
	<li>'81562'</li>
	<li>'83468'</li>
	<li>'6648'</li>
	<li>'9334'</li>
	<li>'81611'</li>
	<li>'84065'</li>
	<li>'6468'</li>
	<li>'83483'</li>
	<li>'79365'</li>
	<li>'1978'</li>
	<li>'3275'</li>
	<li>'3778'</li>
	<li>'4329'</li>
	<li>'29100'</li>
	<li>'28971'</li>
	<li>'79135'</li>
	<li>'23603'</li>
	<li>'57148'</li>
	<li>'8266'</li>
	<li>'4121'</li>
	<li>'221037'</li>
	<li>'55603'</li>
	<li>'2004'</li>
	<li>'387263'</li>
	<li>'57698'</li>
	<li>'9367'</li>
	<li>'84148'</li>
	<li>'90806'</li>
	<li>'5911'</li>
	<li>'100129361'</li>
	<li>'55793'</li>
	<li>'93129'</li>
	<li>'57535'</li>
	<li>'55152'</li>
	<li>'4173'</li>
	<li>'29127'</li>
	<li>'64743'</li>
	<li>'65084'</li>
	<li>'3669'</li>
	<li>'203069'</li>
	<li>'6448'</li>
	<li>'6405'</li>
	<li>'6015'</li>
	<li>'8859'</li>
	<li>'8623'</li>
	<li>'375449'</li>
	<li>'4059'</li>
	<li>'23038'</li>
	<li>'9620'</li>
	<li>'219699'</li>
	<li>'54985'</li>
	<li>'54512'</li>
	<li>'79874'</li>
	<li>'6772'</li>
</ol>




    #
    # over-representation test
    #
    #...@organism 	 hsa 
    #...@ontology 	 KEGG 
    #...@keytype 	 kegg 
    #...@gene 	 chr [1:838] "5573" "7453" "5245" "23450" "6500" "4926" "6427" "813" ...
    #...pvalues adjusted by 'BH' with cutoff <0.01 
    #...6 enriched terms found
    'data.frame':	6 obs. of  9 variables:
     $ ID         : chr  "hsa04110" "hsa05215" "hsa03030" "hsa01521" ...
     $ Description: chr  "Cell cycle" "Prostate cancer" "DNA replication" "EGFR tyrosine kinase inhibitor resistance" ...
     $ GeneRatio  : chr  "20/388" "17/388" "10/388" "14/388" ...
     $ BgRatio    : chr  "124/7387" "97/7387" "36/7387" "79/7387" ...
     $ pvalue     : num  6.03e-06 9.73e-06 1.06e-05 5.23e-05 5.35e-05 ...
     $ p.adjust   : num  0.000987 0.000987 0.000987 0.002996 0.002996 ...
     $ qvalue     : num  0.000812 0.000812 0.000812 0.002466 0.002466 ...
     $ geneID     : chr  "6500/9184/4172/994/4175/4171/1387/10274/8697/902/4616/5591/4176/8881/7043/983/1022/1028/891/4173" "2950/1387/5159/5604/5156/596/4318/3551/367/2260/5595/5295/10000/6935/6655/90993/80310" "5425/4172/4175/4171/10535/5984/2237/4176/54107/4173" "5159/5604/558/5156/596/9470/5595/7422/5295/10000/4763/6655/80310/1978" ...
     $ Count      : int  20 17 10 14 19 16
    #...Citation
      Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
      clusterProfiler: an R package for comparing biological themes among
      gene clusters. OMICS: A Journal of Integrative Biology
      2012, 16(5):284-287 
    



<table>
<thead><tr><th></th><th scope=col>ID</th><th scope=col>Description</th><th scope=col>GeneRatio</th><th scope=col>BgRatio</th><th scope=col>pvalue</th><th scope=col>p.adjust</th><th scope=col>qvalue</th><th scope=col>geneID</th><th scope=col>Count</th></tr></thead>
<tbody>
	<tr><th scope=row>hsa04110</th><td>hsa04110                                                                                        </td><td>Cell cycle                                                                                      </td><td>20/388                                                                                          </td><td>124/7387                                                                                        </td><td>6.027647e-06                                                                                    </td><td>0.0009866915                                                                                    </td><td>0.0008123513                                                                                    </td><td>6500/9184/4172/994/4175/4171/1387/10274/8697/902/4616/5591/4176/8881/7043/983/1022/1028/891/4173</td><td>20                                                                                              </td></tr>
	<tr><th scope=row>hsa05215</th><td>hsa05215                                                                                        </td><td>Prostate cancer                                                                                 </td><td>17/388                                                                                          </td><td>97/7387                                                                                         </td><td>9.732023e-06                                                                                    </td><td>0.0009866915                                                                                    </td><td>0.0008123513                                                                                    </td><td>2950/1387/5159/5604/5156/596/4318/3551/367/2260/5595/5295/10000/6935/6655/90993/80310           </td><td>17                                                                                              </td></tr>
	<tr><th scope=row>hsa03030</th><td>hsa03030                                                                                        </td><td>DNA replication                                                                                 </td><td>10/388                                                                                          </td><td>36/7387                                                                                         </td><td>1.057169e-05                                                                                    </td><td>0.0009866915                                                                                    </td><td>0.0008123513                                                                                    </td><td>5425/4172/4175/4171/10535/5984/2237/4176/54107/4173                                             </td><td>10                                                                                              </td></tr>
	<tr><th scope=row>hsa01521</th><td>hsa01521                                                                                        </td><td>EGFR tyrosine kinase inhibitor resistance                                                       </td><td>14/388                                                                                          </td><td>79/7387                                                                                         </td><td>5.233863e-05                                                                                    </td><td>0.0029956171                                                                                    </td><td>0.0024663164                                                                                    </td><td>5159/5604/558/5156/596/9470/5595/7422/5295/10000/4763/6655/80310/1978                           </td><td>14                                                                                              </td></tr>
	<tr><th scope=row>hsa04068</th><td>hsa04068                                                                                        </td><td>FoxO signaling pathway                                                                          </td><td>19/388                                                                                          </td><td>132/7387                                                                                        </td><td>5.349316e-05                                                                                    </td><td>0.0029956171                                                                                    </td><td>0.0024663164                                                                                    </td><td>7874/5571/10769/1387/5604/901/5106/4616/8660/3551/7043/5595/5295/10000/6655/3643/891/10365/6648 </td><td>19                                                                                              </td></tr>
	<tr><th scope=row>hsa00240</th><td>hsa00240                                                                                        </td><td>Pyrimidine metabolism                                                                           </td><td>16/388                                                                                          </td><td>102/7387                                                                                        </td><td>7.297947e-05                                                                                    </td><td>0.0034057087                                                                                    </td><td>0.0028039482                                                                                    </td><td>5425/4860/6241/7298/5440/7372/5430/9583/4832/54107/953/5435/1635/55718/51728/29922              </td><td>16                                                                                              </td></tr>
</tbody>
</table>




```R
plot(kegg_enrichment)  # ? no plot
```


```R

```
