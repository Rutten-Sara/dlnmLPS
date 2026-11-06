# Spatial boundaries of the 10 Barcelona districts

Spatial boundaries of the 10 administrative districts of Barcelona,
provided as an `sf` object. The dataset is used alongside the simulated
`barcelona_temp` dataset to facilitate spatial visualisation and
examples in the `dlnmLPS` package.

## Usage

``` r
shapefile_bcn
```

## Format

An `sf` object with 10 features and 2 columns:

- District:

  District identifier, ranging from 01 to 10.

- geometry:

  Spatial geometry of each Barcelona district, represented as `POLYGON`
  or `MULTIPOLYGON` features.

The geometries use the ETRS89 / UTM zone 31N coordinate reference system
(EPSG:25831).

## Source

Barcelona City Council, Open Data BCN. *Districtes i barris*.
<https://opendata-ajuntament.barcelona.cat/data/dataset/20170706-districtes-barris/resource/11851135-6919-4dcb-91ed-821e5e87a428>

## Details

The district boundaries correspond to the 10 administrative districts of
Barcelona and can be linked to the `District` variable in
`barcelona_temp`.

Unlike `barcelona_temp`, which contains simulated data, the spatial
boundaries in this dataset originate from the Barcelona City Council
Open Data BCN portal.

## See also

[`barcelona_temp`](https://rutten-sara.github.io/dlnmLPS/reference/barcelona_temp.md)

## Examples

``` r
data(shapefile_bcn)

library(sf)
#> Linking to GEOS 3.14.1, GDAL 3.12.2, PROJ 9.7.1; sf_use_s2() is TRUE
plot(shapefile_bcn["District"])
```
