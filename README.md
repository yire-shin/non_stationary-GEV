An R package for L-moment-based estimation of nonstationary generalized extreme value model

## Installation

Install the latest development version (on GitHub) via `{remotes}`:

``` r
remotes::install_github("yire-shin/non_stationary-GEV")
```

## Getting started

### Estimates

``` r
  library(nsgev)
  data("Trehafod")
  result1 <- nsgev(Trehafod$r1)
  print(result1)

  result2 <- gado.prop_11(Trehafod$r1)
  print(result2)
```

## Learn more

...

## Authors

- Yire Shin; Ph.D, Department of Mathematics and Statistics, Chonnam National University, Republic of Korea, [shinyire87@gmail.com]
- Yonggwan Shin; Ph.D, Research & Development Center, XRAI Inc., Gwangju 61186, Korea, [syg.stat@gmail.com]
- Jeong-Soo Park; Professor, Department of Mathematics and Statistics, Chonnam National University, Republic of Korea, [jspark@jnu.ac.kr]

## Citation

Shin, Y., Shin, Y., Park, J. S. (2024). L-moment-based algorithm for nonstationary generalized extreme value model.
-----
