# We can call peaks independent of ArchR

    Code
      SummarizedExperiment::assays(S4Vectors::metadata(tiles)$summarizedData)[[
        "FragmentCounts"]]
    Output
         PBMCSmall
      C2    117146
      C5    171033

