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
 # I have removed the 6397 duplicated probes, trying to select those under manifest_probe_match==TRUE when possible
 any(duplicated(rownames(Mset)))
 getBeta(Mset)
```


