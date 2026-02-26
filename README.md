# Additional functions for gamlss2 in growth charting (gamlss2charts)

**Author**: Andrew A. Chen, chenandr@musc.edu 

**Maintainer**: Andrew A. Chen, chenandr@musc.edu

**License**: MIT License

`gamlss2charts` includes additional functions based on gamlss2, with backwards compatibility for gamlss. Currently, our package contains `predict_score` for out-of-sample scoring (i.e. centile scoring) based on previous works (see references). 

This package is still a work-in-progress and will be continually supported and updated to include additional functionalities.

## 1. Installation
The latest version can be installed via `remotes` by running the following code

```
# install.packages("remotes")
remotes::install_github("andy1764/gamlss2charts")
```

Then, you can load this package via

```
library(gamlss2charts)
```

## 2. Out-of-sample centile scoring
A toy example for getting out-of-sample centile scores in the iris dataset is shown below

```
library(gamlss2)
ft <- gamlss2(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species | Species,
              family = BCCG(), data = iris[1:100,])
predict_score(ft, iris[-(1:100),], adjust = FALSE, rm.term = "Species") |> hist()
predict_score(ft, iris[-(1:100),], rm.term = "Species") |> hist()
```

## 3. Citations
Please cite the following papers for out-of-sample centile scoring:

> Dinga, R., Fraza, C. J., Bayer, J. M. M., Kia, S. M., Beckmann, C. F., & Marquand, A. F. (2021). Normative modeling of neuroimaging data using generalized additive models of location scale and shape (p. 2021.06.14.448106). bioRxiv. https://doi.org/10.1101/2021.06.14.448106
>
> Bethlehem, R. A. I., Seidlitz, J., White, S. R., Vogel, J. W., Anderson, K. M., Adamson, C., Adler, S., Alexopoulos, G. S., Anagnostou, E., Areces-Gonzalez, A., Astle, D. E., Auyeung, B., Ayub, M., Bae, J., Ball, G., Baron-Cohen, S., Beare, R., Bedford, S. A., Benegal, V., … Alexander-Bloch, A. F. (2022). Brain charts for the human lifespan. Nature, 604(7906), Article 7906. https://doi.org/10.1038/s41586-022-04554-y

