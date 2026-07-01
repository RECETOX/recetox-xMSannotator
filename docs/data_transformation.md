## recetox-aplcms to xMSannotator
The recetox-aplcms tool outputs the peak table in a format which is not directly compatible with recetox-xMSannotator or xMSannotator. A conversion is required. Example file layouts are given below.

### recetox-aplcms Example Output

| feature      | mz   | rt | sample | sample_rt | sample_intensity |
|--------------|------|----|--------|-----------|------------------|
|1	|110.034534364169 |	246.052915606977 |	Tribrid_201005_001-MeOH_NEG_MU |NA|	0 |
|1|	110.034534364169|	246.052915606977|	Tribrid_201005_002-Norman_NEG_MU|	729.439375046056|	256906.674516574|
|2|	113.007573530083|	370.352087874597|	Tribrid_201005_001-MeOH_NEG_MU|	317.947801280858|	2406770.19370756|
| ... | ... | ... | ... | ...| ... |

### recetox-xMSannotator Example Input

| peak      | mz   | rt | Tribrid_201005_001-MeOH_NEG_MU | Tribrid_201005_002-Norman_NEG_MU |
|-----------|------|----|------------------------------------------|--------------------------------------------|
| 1 |110.034534364169 | 246.052915606977 | 0 | 256906.674516574 |
| ... | ... | ... | ... | ... |

### xMSannotator Example Input
| mz   | time | intensity_Tribrid_201005_001-MeOH_NEG_MU | intensity_Tribrid_201005_002-Norman_NEG_MU |
|------|----|------------------------------------------|--------------------------------------------|
|110.034534364169 | 246.052915606977 | 0 | 256906.674516574 |
| ... | ... | ... | ... |
### R Sample Code

``` R
#' @param file Path to parquet file to be loaded into memory
#'
#' @export
load_parquet <- function(file) {
    message <- paste("The file", toString(file), "seems not to be a valid Parquet file.")
    rlang::with_handlers(
        arrow::read_parquet(file),
        error = ~ rlang::abort(message, parent = .)
    )
}


#' @param peak_table Peak table from recetox.aplcms package.
#'
#' @export
#' @import dplyr
rcx_aplcms_to_xmsannotator <- function(peak_table) {
    col_base <- c("feature", "mz", "rt")
    output_table <- peak_table %>% distinct(across(any_of(col_base)))

    for (level in levels(peak_table$sample)) {
        subdata <- peak_table %>%
            filter(sample == level) %>%
            select(any_of(c(col_base, "sample_intensity"))) %>%
            rename(!!paste0("intensity_", level) := "sample_intensity")
        output_table <- inner_join(output_table, subdata, by = col_base)
    }

    output_table <- output_table  %>% select(-feature) %>% rename(time = rt)
    return(output_table)
}

```

## recetox-aplcms to recetox-xMSannotator
This is the same as to the recetox-xMSannotator, but the `feature` column has to be renamed to peak, and `rt` can stay the same.

### R Sample Code

```R
#' @param peak_table Peak table from recetox.aplcms package.
#'
#' @export
#' @import dplyr
rcx_aplcms_to_rcx_xmsannotator <- function(peak_table) {
    col_base <- c("feature", "mz", "rt")
    output_table <- peak_table %>% distinct(across(any_of(col_base)))

    for (level in levels(peak_table$sample)) {
        subdata <- peak_table %>%
            filter(sample == level) %>%
            select(any_of(c(col_base, "sample_intensity"))) %>%
            rename(!!level := "sample_intensity")
        output_table <- inner_join(output_table, subdata, by = col_base)
    }

    output_table <- output_table %>% rename(peak = feature)
    return(output_table)
}

```

## HMDB: Rda to Parquet
These are simple instructions how to convert `RDA` file containing the HMDB annotation database (used in xMSannotator) to `parquet` format.

The important steps are defining `filter_keys` and `target_keys`, which are vectors containing keywords to be selected and subsequently renamed, in the respective order.

### R Sample Code

```R
save_parquet <- function(data, filename) {
  invisible(arrow::write_parquet(data, filename))
}

load_rda <- function(filename) {
  get(load(filename))
}

convert_rda_to_parquet <- function(rda, filter_keys, target_keys) {
  target_data <- rda[filter_keys]
  names(target_data) <- target_keys
  return(target_data)
}

rda <- load_rda("filename.rda")
filter_keys <- c("Name", "InChiKey", "Formula", "MonoisotopicMass")
target_keys <- c("iupac_name", "iupac_inchikey", "molecular_formula", "monoisotopic_mass")
result <- convert_rda_to_parquet(rda, filter_keys, target_keys)
save_parquet(result, "filename.parquet")

```
