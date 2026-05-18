# We can call peaks independent of ArchR

    Code
      SummarizedExperiment::assays(S4Vectors::metadata(tiles)$summarizedData)[[
        "CellCounts"]]
    Output
         PBMCSmall
      C2       152
      C5       201

