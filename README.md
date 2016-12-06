
isin
----

Support the NASA / GlobColour / CCI ISIN grid used for MODIS L3BIN satellite products. This package can be installed as follow:

``` r
devtools::install_github("PMassicotte/isin")
```

Examples
--------

``` r
library(isin)

## Create test dataset
bin_num = c(245535, 245536, 247290, 249046, 249047, 250809)
```

### Convert bin index to latlon

The `bin2latlon()` function converts a MODIS bin index to longitude and latitude.

``` r
df <- bin2latlon(bin_num)
df
#>   longitude  latitude
#> 1 -162.2057 -78.31250
#> 2 -162.0000 -78.31250
#> 3 -161.2415 -78.27084
#> 4 -161.3159 -78.22916
#> 5 -161.1117 -78.22916
#> 6 -161.3793 -78.18750
```

### Convert latlon to bin index

The `latlon2bin()` function converts a pair of longitude/latitude to a MODIS bin index.

``` r
latlon2bin(lat = df$latitude, lon = df$longitude)
#>   bin_index
#> 1    245535
#> 2    245536
#> 3    247290
#> 4    249046
#> 5    249047
#> 6    250809
```

### Credit

The C/C++ code is almost entierly based on [oceancolor algorithm](http://oceancolor.gsfc.nasa.gov/staff/norman/swreadl3b/swreadl3b.c). Credit should be given to them.
