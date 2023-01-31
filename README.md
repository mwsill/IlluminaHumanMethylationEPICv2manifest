# IlluminaHumanMethylationEPICv2manifest

```r
 library(devtools)
 install_github("mwsill/IlluminaHumanMethylationEPICv2manifest") 
 install_github("mwsill/minfi")
 
 library(minfi)
 RGset <- read.metharray("./206891110005/206891110005_R02C01")
 RGset
 Mset <- preprocessIllumina(RGset)
 Mset
 getBeta(Mset)
 
```


