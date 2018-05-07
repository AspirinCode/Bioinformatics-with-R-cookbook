
# 8、Enriching RNAseq data with GO terms

The RNAseq data coming out of NGS provides great detail on the cellular transcriptional process. Besides measuring expression levels of the transcripts, they provide information on alternative splicing, allele-specific expressions, and so on. Thus, the RNAseq data gives a more comprehensive picture of differential expression in cells. However, adding the functional aspects can refine this analysis further with statistical robustness. This section aims to
demonstrate the highlighting of the RNAseq data in terms of GO terms.<br>
Take the data from the `goseq` package for analysis.<br>

来自NGS的RNAseq数据提供了细胞转录过程的细节。除了测量转录本的表达水平之外，它们还提供有关可变剪接，等位基因特异性表达等的信息。因此，RNAseq数据提供了细胞中差异表达的更全面的图片。然而，增加它的功能方面可以进一步细化这种分析，并具有统计稳健性。这部分旨在展示GO术语中RNAseq数据的表示。<br>


```R
# 1、Loading the required packages, `goseq` and `edgeR`
library(goseq)
library(edgeR)
supportedGenomes() # Take a look at the genomes that are supported by the package
```


<table>
<thead><tr><th></th><th scope=col>db</th><th scope=col>species</th><th scope=col>date</th><th scope=col>name</th><th scope=col>AvailableGeneIDs</th></tr></thead>
<tbody>
	<tr><th scope=row>2</th><td>hg19                                                                                                                             </td><td>Human                                                                                                                            </td><td>Feb. 2009                                                                                                                        </td><td>Genome Reference Consortium GRCh37                                                                                               </td><td>ccdsGene,ensGene,exoniphy,geneSymbol,knownGene,nscanGene,refGene,xenoRefGene                                                     </td></tr>
	<tr><th scope=row>3</th><td>hg18                                                                                                                             </td><td>Human                                                                                                                            </td><td>Mar. 2006                                                                                                                        </td><td>NCBI Build 36.1                                                                                                                  </td><td>acembly,acescan,ccdsGene,ensGene,exoniphy,geneSymbol,geneid,genscan,knownGene,knownGeneOld3,refGene,sgpGene,sibGene,xenoRefGene  </td></tr>
	<tr><th scope=row>4</th><td>hg17                                                                                                                             </td><td>Human                                                                                                                            </td><td>May 2004                                                                                                                         </td><td>NCBI Build 35                                                                                                                    </td><td>acembly,acescan,ccdsGene,ensGene,exoniphy,geneSymbol,geneid,genscan,knownGene,refGene,sgpGene,vegaGene,vegaPseudoGene,xenoRefGene</td></tr>
	<tr><th scope=row>5</th><td>hg16                                                                                                                             </td><td>Human                                                                                                                            </td><td>Jul. 2003                                                                                                                        </td><td>NCBI Build 34                                                                                                                    </td><td>acembly,ensGene,exoniphy,geneSymbol,geneid,genscan,knownGene,refGene,sgpGene                                                     </td></tr>
	<tr><th scope=row>19</th><td>felCat3                                                                                                                          </td><td>Cat                                                                                                                              </td><td>Mar. 2006                                                                                                                        </td><td>Broad Institute Release 3                                                                                                        </td><td>ensGene,geneSymbol,geneid,genscan,nscanGene,refGene,sgpGene,xenoRefGene                                                          </td></tr>
	<tr><th scope=row>23</th><td>panTro2                                                                                                                          </td><td>Chimp                                                                                                                            </td><td>Mar. 2006                                                                                                                        </td><td>CGSC Build 2.1                                                                                                                   </td><td>ensGene,geneSymbol,genscan,nscanGene,refGene,xenoRefGene                                                                         </td></tr>
	<tr><th scope=row>24</th><td>panTro1                                                                                                                          </td><td>Chimp                                                                                                                            </td><td>Nov. 2003                                                                                                                        </td><td>CGSC Build 1.1                                                                                                                   </td><td>ensGene,geneid,genscan,xenoRefGene                                                                                               </td></tr>
	<tr><th scope=row>31</th><td>bosTau4                                                                                                                          </td><td>Cow                                                                                                                              </td><td>Oct. 2007                                                                                                                        </td><td>Baylor College of Medicine HGSC Btau_4.0                                                                                         </td><td>ensGene,geneSymbol,genscan,nscanGene,refGene                                                                                     </td></tr>
	<tr><th scope=row>32</th><td>bosTau3                                                                                                                          </td><td>Cow                                                                                                                              </td><td>Aug. 2006                                                                                                                        </td><td>Baylor College of Medicine HGSC Btau_3.1                                                                                         </td><td>ensGene,geneSymbol,geneid,genscan,refGene,sgpGene                                                                                </td></tr>
	<tr><th scope=row>33</th><td>bosTau2                                                                                                                          </td><td>Cow                                                                                                                              </td><td>Mar. 2005                                                                                                                        </td><td>Baylor College of Medicine HGSC Btau_2.0                                                                                         </td><td>geneSymbol,geneid,genscan,refGene,sgpGene                                                                                        </td></tr>
	<tr><th scope=row>36</th><td>canFam2                                                                                                                          </td><td>Dog                                                                                                                              </td><td>May 2005                                                                                                                         </td><td>Broad Institute v2.0                                                                                                             </td><td>ensGene,geneSymbol,genscan,nscanGene,refGene,xenoRefGene                                                                         </td></tr>
	<tr><th scope=row>37</th><td>canFam1                                                                                                                          </td><td>Dog                                                                                                                              </td><td>Jul. 2004                                                                                                                        </td><td>Broad Institute v1.0                                                                                                             </td><td>ensGene,geneSymbol,genscan,nscanGene,refGene,xenoRefGene                                                                         </td></tr>
	<tr><th scope=row>39</th><td>loxAfr3                                                                                                                          </td><td>Elephant                                                                                                                         </td><td>Jul. 2009                                                                                                                        </td><td>Broad Institute LoxAfr3                                                                                                          </td><td>xenoRefGene                                                                                                                      </td></tr>
	<tr><th scope=row>50</th><td>cavPor3                                                                                                                          </td><td>Guinea pig                                                                                                                       </td><td>Feb. 2008                                                                                                                        </td><td>Broad Institute cavPor3                                                                                                          </td><td>ensGene,genscan,nscanGene,xenoRefGene                                                                                            </td></tr>
	<tr><th scope=row>53</th><td>equCab2                                                                                                                          </td><td>Horse                                                                                                                            </td><td>Sep. 2007                                                                                                                        </td><td>Broad Institute EquCab2                                                                                                          </td><td>ensGene,geneSymbol,nscanGene,refGene,xenoRefGene                                                                                 </td></tr>
	<tr><th scope=row>54</th><td>equCab1                                                                                                                          </td><td>Horse                                                                                                                            </td><td>Jan. 2007                                                                                                                        </td><td>Broad Institute EquCab1                                                                                                          </td><td>geneSymbol,geneid,nscanGene,refGene,sgpGene                                                                                      </td></tr>
	<tr><th scope=row>59</th><td>calJac1                                                                                                                          </td><td>Marmoset                                                                                                                         </td><td>Jun. 2007                                                                                                                        </td><td>WUSTL Callithrix_jacchus-v2.0.2                                                                                                  </td><td>genscan,nscanGene,xenoRefGene                                                                                                    </td></tr>
	<tr><th scope=row>64</th><td>mm9                                                                                                                              </td><td>Mouse                                                                                                                            </td><td>Jul. 2007                                                                                                                        </td><td>NCBI Build 37                                                                                                                    </td><td>acembly,ccdsGene,ensGene,exoniphy,geneSymbol,geneid,genscan,knownGene,nscanGene,refGene,sgpGene,xenoRefGene                      </td></tr>
	<tr><th scope=row>65</th><td>mm8                                                                                                                              </td><td>Mouse                                                                                                                            </td><td>Feb. 2006                                                                                                                        </td><td>NCBI Build 36                                                                                                                    </td><td>ccdsGene,ensGene,geneSymbol,geneid,genscan,knownGene,nscanGene,refGene,sgpGene,sibGene,xenoRefGene                               </td></tr>
	<tr><th scope=row>66</th><td>mm7                                                                                                                              </td><td>Mouse                                                                                                                            </td><td>Aug. 2005                                                                                                                        </td><td>NCBI Build 35                                                                                                                    </td><td>ensGene,geneSymbol,geneid,genscan,knownGene,refGene,sgpGene,xenoRefGene                                                          </td></tr>
	<tr><th scope=row>71</th><td>monDom5                                                                                                                          </td><td>Opossum                                                                                                                          </td><td>Oct. 2006                                                                                                                        </td><td>Broad Institute release MonDom5                                                                                                  </td><td>ensGene,geneSymbol,genscan,nscanGene,refGene,xenoRefGene                                                                         </td></tr>
	<tr><th scope=row>72</th><td>monDom4                                                                                                                          </td><td>Opossum                                                                                                                          </td><td>Jan. 2006                                                                                                                        </td><td>Broad Institute release MonDom4                                                                                                  </td><td>ensGene,geneSymbol,genscan,nscanGene,refGene,xenoRefGene                                                                         </td></tr>
	<tr><th scope=row>73</th><td>monDom1                                                                                                                          </td><td>Opossum                                                                                                                          </td><td>Oct. 2004                                                                                                                        </td><td>Broad Institute release MonDom1                                                                                                  </td><td>genscan                                                                                                                          </td></tr>
	<tr><th scope=row>74</th><td>ponAbe2                                                                                                                          </td><td>Orangutan                                                                                                                        </td><td>Jul. 2007                                                                                                                        </td><td>WUSTL Pongo_albelii-2.0.2                                                                                                        </td><td>ensGene,geneSymbol,genscan,nscanGene,refGene,xenoRefGene                                                                         </td></tr>
	<tr><th scope=row>82</th><td>ornAna1                                                                                                                          </td><td>Platypus                                                                                                                         </td><td>Mar. 2007                                                                                                                        </td><td>WUSTL v5.0.1                                                                                                                     </td><td>ensGene,geneSymbol,refGene,xenoRefGene                                                                                           </td></tr>
	<tr><th scope=row>87</th><td>rn4                                                                                                                              </td><td>Rat                                                                                                                              </td><td>Nov. 2004                                                                                                                        </td><td>Baylor College of Medicine HGSC v3.4                                                                                             </td><td>ensGene,geneSymbol,geneid,genscan,knownGene,nscanGene,refGene,sgpGene,xenoRefGene                                                </td></tr>
	<tr><th scope=row>88</th><td>rn3                                                                                                                              </td><td>Rat                                                                                                                              </td><td>Jun. 2003                                                                                                                        </td><td>Baylor College of Medicine HGSC v3.1                                                                                             </td><td>ensGene,geneSymbol,geneid,genscan,knownGene,nscanGene,refGene,sgpGene,xenoRefGene                                                </td></tr>
	<tr><th scope=row>91</th><td>rheMac2                                                                                                                          </td><td>Rhesus                                                                                                                           </td><td>Jan. 2006                                                                                                                        </td><td>Baylor College of Medicine HGSC v1.0 Mmul_051212                                                                                 </td><td>ensGene,geneSymbol,geneid,nscanGene,refGene,sgpGene,xenoRefGene                                                                  </td></tr>
	<tr><th scope=row>114</th><td>galGal3                                                                                                                          </td><td>Chicken                                                                                                                          </td><td>May 2006                                                                                                                         </td><td>WUSTL Gallus-gallus-2.1                                                                                                          </td><td>ensGene,geneSymbol,genscan,nscanGene,refGene,xenoRefGene                                                                         </td></tr>
	<tr><th scope=row>115</th><td>galGal2                                                                                                                          </td><td>Chicken                                                                                                                          </td><td>Feb. 2004                                                                                                                        </td><td>WUSTL Gallus-gallus-1.0                                                                                                          </td><td>ensGene,geneSymbol,geneid,genscan,refGene,sgpGene                                                                                </td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>161</th><td>droEre1                                    </td><td>D. erecta                                  </td><td>Aug. 2005                                  </td><td>Agencourt Arachne release                  </td><td>genscan,xenoRefGene                        </td></tr>
	<tr><th scope=row>162</th><td>droGri1                                    </td><td>D. grimshawi                               </td><td>Aug. 2005                                  </td><td>Agencourt Arachne release                  </td><td>genscan,xenoRefGene                        </td></tr>
	<tr><th scope=row>164</th><td>dm3                                        </td><td>D. melanogaster                            </td><td>Apr. 2006                                  </td><td>BDGP Release 5                             </td><td>geneSymbol,nscanPasaGene,refGene           </td></tr>
	<tr><th scope=row>165</th><td>dm2                                        </td><td>D. melanogaster                            </td><td>Apr. 2004                                  </td><td>BDGP Release 4                             </td><td>geneSymbol,geneid,genscan,nscanGene,refGene</td></tr>
	<tr><th scope=row>166</th><td>dm1                                        </td><td>D. melanogaster                            </td><td>Jan. 2003                                  </td><td>BDGP Release 3                             </td><td>geneSymbol,genscan,refGene                 </td></tr>
	<tr><th scope=row>167</th><td>droMoj2                                    </td><td>D. mojavensis                              </td><td>Aug. 2005                                  </td><td>Agencourt Arachne release                  </td><td>genscan,xenoRefGene                        </td></tr>
	<tr><th scope=row>168</th><td>droMoj1                                    </td><td>D. mojavensis                              </td><td>Aug. 2004                                  </td><td>Agencourt Arachne release                  </td><td>geneid,genscan,xenoRefGene                 </td></tr>
	<tr><th scope=row>169</th><td>droPer1                                    </td><td>D. persimilis                              </td><td>Oct. 2005                                  </td><td>Broad Institute release                    </td><td>genscan,xenoRefGene                        </td></tr>
	<tr><th scope=row>170</th><td>dp3                                        </td><td>D. pseudoobscura                           </td><td>Nov. 2004                                  </td><td>FlyBase Release 1.0                        </td><td>geneid,genscan,xenoRefGene                 </td></tr>
	<tr><th scope=row>171</th><td>dp2                                        </td><td>D. pseudoobscura                           </td><td>Aug. 2003                                  </td><td>Baylor College of Medicine HGSC Freeze 1   </td><td>genscan,xenoRefGene                        </td></tr>
	<tr><th scope=row>172</th><td>droSec1                                    </td><td>D. sechellia                               </td><td>Oct. 2005                                  </td><td>Broad Institute Release 1.0                </td><td>genscan,xenoRefGene                        </td></tr>
	<tr><th scope=row>173</th><td>droSim1                                    </td><td>D. simulans                                </td><td>Apr. 2005                                  </td><td>WUSTL Release 1.0                          </td><td>geneid,genscan,xenoRefGene                 </td></tr>
	<tr><th scope=row>174</th><td>droVir2                                    </td><td>D. virilis                                 </td><td>Aug. 2005                                  </td><td>Agencourt Arachne release                  </td><td>genscan,xenoRefGene                        </td></tr>
	<tr><th scope=row>175</th><td>droVir1                                    </td><td>D. virilis                                 </td><td>Jul. 2004                                  </td><td>Agencourt Arachne release                  </td><td>geneid,genscan,xenoRefGene                 </td></tr>
	<tr><th scope=row>176</th><td>droYak2                                    </td><td>D. yakuba                                  </td><td>Nov. 2005                                  </td><td>WUSTL Release 2.0                          </td><td>genscan,xenoRefGene                        </td></tr>
	<tr><th scope=row>177</th><td>droYak1                                    </td><td>D. yakuba                                  </td><td>Apr. 2004                                  </td><td>WUSTL Release 1.0                          </td><td>geneid,genscan,xenoRefGene                 </td></tr>
	<tr><th scope=row>178</th><td>caePb2                                     </td><td>C. brenneri                                </td><td>Feb. 2008                                  </td><td>WUSTL 6.0.1                                </td><td>xenoRefGene                                </td></tr>
	<tr><th scope=row>179</th><td>caePb1                                     </td><td>C. brenneri                                </td><td>Jan. 2007                                  </td><td>WUSTL 4.0                                  </td><td>xenoRefGene                                </td></tr>
	<tr><th scope=row>180</th><td>cb3                                        </td><td>C. briggsae                                </td><td>Jan. 2007                                  </td><td>WUSTL Cb3                                  </td><td>xenoRefGene                                </td></tr>
	<tr><th scope=row>181</th><td>cb1                                        </td><td>C. briggsae                                </td><td>Jul. 2002                                  </td><td>WormBase v. cb25.agp8                      </td><td>xenoRefGene                                </td></tr>
	<tr><th scope=row>184</th><td>ce6                                        </td><td>C. elegans                                 </td><td>May 2008                                   </td><td>WormBase v. WS190                          </td><td>ensGene,geneSymbol,refGene,xenoRefGene     </td></tr>
	<tr><th scope=row>185</th><td>ce4                                        </td><td>C. elegans                                 </td><td>Jan. 2007                                  </td><td>WormBase v. WS170                          </td><td>geneSymbol,refGene,xenoRefGene             </td></tr>
	<tr><th scope=row>186</th><td>ce2                                        </td><td>C. elegans                                 </td><td>Mar. 2004                                  </td><td>WormBase v. WS120                          </td><td>geneSymbol,geneid,refGene                  </td></tr>
	<tr><th scope=row>187</th><td>caeJap1                                    </td><td>C. japonica                                </td><td>Mar. 2008                                  </td><td>WUSTL 3.0.2                                </td><td>xenoRefGene                                </td></tr>
	<tr><th scope=row>188</th><td>caeRem3                                    </td><td>C. remanei                                 </td><td>May 2007                                   </td><td>WUSTL 15.0.1                               </td><td>xenoRefGene                                </td></tr>
	<tr><th scope=row>189</th><td>caeRem2                                    </td><td>C. remanei                                 </td><td>Mar. 2006                                  </td><td>WUSTL 1.0                                  </td><td>xenoRefGene                                </td></tr>
	<tr><th scope=row>190</th><td>priPac1                                    </td><td>P. pacificus                               </td><td>Feb. 2007                                  </td><td>WUSTL 5.0                                  </td><td>xenoRefGene                                </td></tr>
	<tr><th scope=row>191</th><td>aplCal1                                    </td><td>Sea Hare                                   </td><td>Sep. 2008                                  </td><td>Broad Release Aplcal2.0                    </td><td>xenoRefGene                                </td></tr>
	<tr><th scope=row>193</th><td>sacCer2                                    </td><td>Yeast                                      </td><td>June 2008                                  </td><td>SGD June 2008 sequence                     </td><td>ensGene                                    </td></tr>
	<tr><th scope=row>194</th><td>sacCer1                                    </td><td>Yeast                                      </td><td>Oct. 2003                                  </td><td>SGD 1 Oct 2003 sequence                    </td><td>ensGene                                    </td></tr>
</tbody>
</table>



2、Load the data used for analysis from the `goseq` package, and the first four columns in the data are controls and the last three are the treatment samples. Assign these attributes to the data through `factor()`.<br>
载入`goseq`中用来分析的数据，前四列是对照组，后三列是实验组。通过`factor()`为这些数据设置属性值。


```R
myData <- read.table(system.file("extdata", "Li_sum.txt", package="goseq"), 
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE,row.names=1)
head(myData)
# Alternatively, you can read the data directly from the file provided with the code files
# myData <- read.table("path/to/code/files/Li_sum.txt", sep ='\t', header = TRUE, stringsAsFactors = FALSE,row.names=1)
myTreat <- factor(rep(c("Control","Treatment"),times = c(4,3)))
```


<table>
<thead><tr><th></th><th scope=col>lane1</th><th scope=col>lane2</th><th scope=col>lane3</th><th scope=col>lane4</th><th scope=col>lane5</th><th scope=col>lane6</th><th scope=col>lane8</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000215688</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000215689</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000220823</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000242499</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000224938</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000239242</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>



3、Create a `DGElist` object via the `edgeR` library using all the count data and treatment information, and then  estimate the dispersion in the data followed by a Fisher's exact test.<br>
使用所有的计数数据和处理信息，通过`edgeR`库创建`DGElist`对象，然后估计数据中的离差，然后进行Fisher精确检验。


```R
myDG <- DGEList(myData,lib.size = colSums(myData),group = myTreat)    #create a `DGElist` object
myDG
```


<dl>
	<dt>$counts</dt>
		<dd><table>
<thead><tr><th></th><th scope=col>lane1</th><th scope=col>lane2</th><th scope=col>lane3</th><th scope=col>lane4</th><th scope=col>lane5</th><th scope=col>lane6</th><th scope=col>lane8</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000215688</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000215689</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000220823</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000242499</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000224938</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000239242</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000243140</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000240187</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000241444</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000242468</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000241826</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000239150</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000220023</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000238303</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000239060</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000238664</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000238600</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000234317</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000225412</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000233023</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000231129</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000233149</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000225459</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000226732</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000223844</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000229493</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000234100</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000234372</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000224230</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000225859</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>ENSG00000229284</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206437</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206308</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000232902</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235223</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000236934</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206439</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235757</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000237559</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000224044</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000228104</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000183574</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206510</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000243594</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235814</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206458</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000215425</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206366</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000229353</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000242685</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000221974</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235630</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000235952</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000229199</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000232345</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000231030</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000206296</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000212866</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000228904</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>ENSG00000226795</th><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td></tr>
</tbody>
</table>
</dd>
	<dt>$samples</dt>
		<dd><table>
<thead><tr><th></th><th scope=col>group</th><th scope=col>lib.size</th><th scope=col>norm.factors</th></tr></thead>
<tbody>
	<tr><th scope=row>lane1</th><td>Control  </td><td>1178832  </td><td>1        </td></tr>
	<tr><th scope=row>lane2</th><td>Control  </td><td>1384945  </td><td>1        </td></tr>
	<tr><th scope=row>lane3</th><td>Control  </td><td>1716355  </td><td>1        </td></tr>
	<tr><th scope=row>lane4</th><td>Control  </td><td>1767927  </td><td>1        </td></tr>
	<tr><th scope=row>lane5</th><td>Treatment</td><td>2127868  </td><td>1        </td></tr>
	<tr><th scope=row>lane6</th><td>Treatment</td><td>2142158  </td><td>1        </td></tr>
	<tr><th scope=row>lane8</th><td>Treatment</td><td> 816171  </td><td>1        </td></tr>
</tbody>
</table>
</dd>
</dl>




```R
myDisp <- estimateCommonDisp(myDG)                                    #estimate the disperision
mytest <- exactTest(myDisp)                                           #Fisher's exact test 
```

4、Use the genes from this analysis for enrichment. So, extract the genes with the desired p-values and log fold change condition with their corresponding gene names, creating a named vector.<br>
The vector created in the previous step bears the status of every gene in terms of DE or non-DE tags. Check the corresponding numbers by `table()`.<br>
使用此分析中的基因进行富集。因此，用相应的基因名称提取具有所需p值和对数倍数变化条件的基因，创建一个命名好了的变量。
上一步创建的变量在DE或非DE标记方面表示每个基因的状态。可以通过table（）检查相应的数字。


```R
myTags <- as.integer(p.adjust(mytest$table$PValue[mytest$table$logFC!=0], method="BH") < 0.05)
names(myTags) <- row.names(mytest$table[mytest$table$logFC!=0,])
table(myTags)
```


    myTags
        0     1 
    19535  3208 


5、Compute the probability weighting function for a set of genes based on their status, and use the results obtained to compute the enrichment. Take a look at the enrichment.<br>
根据它们的状态计算一组基因的概率加权函数，并使用获得的结果计算富集程度。查看富集结果。


```R
wtFunc <- nullp(myTags,"hg19","ensGene")      #Compute the probability weighting function
head(wtFunc)
```

    Loading hg19 length data...
    Warning message in pcls(G):
    "initial point very close to some inequality constraints"


<table>
<thead><tr><th></th><th scope=col>DEgenes</th><th scope=col>bias.data</th><th scope=col>pwf</th></tr></thead>
<tbody>
	<tr><th scope=row>ENSG00000230758</th><td>0         </td><td> 247      </td><td>0.03757470</td></tr>
	<tr><th scope=row>ENSG00000182463</th><td>0         </td><td>3133      </td><td>0.20436865</td></tr>
	<tr><th scope=row>ENSG00000124208</th><td>0         </td><td>1978      </td><td>0.16881769</td></tr>
	<tr><th scope=row>ENSG00000230753</th><td>0         </td><td> 466      </td><td>0.06927243</td></tr>
	<tr><th scope=row>ENSG00000224628</th><td>0         </td><td>1510      </td><td>0.15903532</td></tr>
	<tr><th scope=row>ENSG00000125835</th><td>0         </td><td> 954      </td><td>0.12711992</td></tr>
</tbody>
</table>




[output](https://github.com/Chengshu21/Chapter-8-Analyzing-NGS-Data/blob/master/MD/pic/output_8、Enriching%20RNAseq%20data%20with%20GO%20terms.png)



```R
myEnrich_wall <- goseq(wtFunc,"hg19","ensGene", test.cats=c("GO:BP"))    # compute the enrichment
head(myEnrich_wall)                                                      # take a look 
```

    Fetching GO annotations...
    For 9616 genes, we could not find any categories. These genes will be excluded.
    To force their use, please run with use_genes_without_cat=TRUE (see documentation).
    This was the default behavior for version 1.15.1 and earlier.
    Calculating the p-values...
    'select()' returned 1:1 mapping between keys and columns
    


<table>
<thead><tr><th></th><th scope=col>category</th><th scope=col>over_represented_pvalue</th><th scope=col>under_represented_pvalue</th><th scope=col>numDEInCat</th><th scope=col>numInCat</th><th scope=col>term</th><th scope=col>ontology</th></tr></thead>
<tbody>
	<tr><th scope=row>62</th><td>GO:0000278                                                 </td><td>1.824185e-09                                               </td><td>1                                                          </td><td>254                                                        </td><td> 888                                                       </td><td>mitotic cell cycle                                         </td><td>BP                                                         </td></tr>
	<tr><th scope=row>1516</th><td>GO:0006613                                                 </td><td>2.420576e-08                                               </td><td>1                                                          </td><td> 35                                                        </td><td>  96                                                       </td><td>cotranslational protein targeting to membrane              </td><td>BP                                                         </td></tr>
	<tr><th scope=row>1642</th><td>GO:0006793                                                 </td><td>2.976380e-08                                               </td><td>1                                                          </td><td>614                                                        </td><td>2508                                                       </td><td>phosphorus metabolic process                               </td><td>BP                                                         </td></tr>
	<tr><th scope=row>1517</th><td>GO:0006614                                                 </td><td>3.084971e-08                                               </td><td>1                                                          </td><td> 33                                                        </td><td>  91                                                       </td><td>SRP-dependent cotranslational protein targeting to membrane</td><td>BP                                                         </td></tr>
	<tr><th scope=row>1801</th><td>GO:0007049                                                 </td><td>4.390085e-08                                               </td><td>1                                                          </td><td>406                                                        </td><td>1570                                                       </td><td>cell cycle                                                 </td><td>BP                                                         </td></tr>
	<tr><th scope=row>6943</th><td>GO:0045047                                                 </td><td>5.952232e-08                                               </td><td>1                                                          </td><td> 35                                                        </td><td>  99                                                       </td><td>protein targeting to ER                                    </td><td>BP                                                         </td></tr>
</tbody>
</table>



The above table shows the GO enrichment of the data.

6、Use the `GO.db` package to look at the meaning of a GO category.<br>
使用`GO.db` 包老查看GO分类的意思。


```R
library(GO.db)
GOTERM[[myEnrich_wall$category[1]]]
```


    GOID: GO:0000278
    Term: mitotic cell cycle
    Ontology: BP
    Definition: Progression through the phases of the mitotic cell cycle,
        the most common eukaryotic cell cycle, which canonically comprises
        four successive phases called G1, S, G2, and M and includes
        replication of the genome and the subsequent segregation of
        chromosomes into daughter cells. In some variant cell cycles
        nuclear replication or nuclear division may not be followed by cell
        division, or G1 and G2 phases may be absent.
    Synonym: mitosis
    Synonym: GO:0007067
    Secondary: GO:0007067


This part is similar to the previous part, `dealing with the edgeR package`. Till step 5, here computed the genes or tags of interest and we have the Ensembl IDs for the tags. The `goseq` function fetched the GO category for the Ensembl genes obtained in step 5. In this part, we used only the biological process category, but for other categories, we can use the values `GO:CC` or `GO:MF`. Besides fetching the GO annotation, the `goseq` function also computes the enrichment in terms of over and under representation of p-values. Finally, the use of the `GO.db` package gives the actual GO term details for the categories retrieved.<br>
这部分与前一部分“处理edgeR软件包”相似。直到步骤5，在这里计算感兴趣的基因或标签，并且我们会有标签的Ensembl ID。`goseq`函数为步骤5中获得的Ensembl基因提取了GO类别。在本部分中，我们仅使用生物过程类别，但对于其他类别，我们可以使用值`GO：CC`或`GO：MF`。除了提取GO注释之外，`goseq`函数还可以计算p值的过度表示和不足表示的富集。最后，使用`GO.db`包给出了检索类别的实际GO术语详细信息。


```R

```
