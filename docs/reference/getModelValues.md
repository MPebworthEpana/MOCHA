# Get a data.frame of model values from the output of linear modeling

\`r lifecycle::badge("deprecated")\` This function is deprecated -
improved modeling functions can be found in the package "ChAI" at
https://github.com/aifimmunology/ChAI `getModelValues` extracts a
data.frame of model values (slope, significance, and std.error) for a
given factor from the SummarizedExperiment output of runZIGLMM.

## Usage

``` r
getModelValues(object, specificVariable)
```

## Arguments

- object:

  A SummarizedExperiment object generated from runZIGLMM.

- specificVariable:

  A string, describing the factor of influence.

## Value

A data.frame of slopes, significance, and standard error for one factor.

## Examples

``` r
if (FALSE) { # \dontrun{
age_df <- getModelValues(runZIGLMM_output, "Age")
} # }
```
