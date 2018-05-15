# simBasis
misc code to simulate basis scenarios and test cupcake R library and method 15/04/2018. Olly Burren

This code is under active development and is not for public use.

`createAnnoFile.R` - This creates a table of SNPs from which to simulate including MAF's positions etc. It is used by software so the getting of haplotypes using `bcftools` is much more efficient. We also can use it to select SNPs that are in a current basis so that we can simulate resolution. I added some addition code so that we don't necessarily need to have the CV in the basis but can still simulate it as long as within r2 of 0.8.

`createScenarioYML.R` - This creates a YAML file that describes a set of simulations. It is envisaged that these are run using `simFullGWAS.R`. Currently it implements a scenario where there is no sharing between basis diseases. Diseases to be projected are as follows:

* Same as GWAS10 but with different case control sizes
* 50% like GWAS10 and 50% like GWAS1 with different case control sizes
* Random randomly selected CV's that may or may not be CV's in the basis already.



`simFullGWAS.R` - This uses Chris Wallace and Mary Fortunes `simGWAS` R package to simulate a GWAS based on a YML file (`createScenarioYML.R`). It then use `cupcake` to compute various PCA weightings etc. Currently it works only for one scenario see above.

`ldFilter.R` - work with `createAnnoFile.R` so that we include variants not in basis but in LD as possible CV.

`malhalanobis.R` - code for me to understand malhalanobis distance. Not used in sims
