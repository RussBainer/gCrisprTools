Changes in version 2.0.0 (2021-05-18)
### *MAJOR UPDATE Adding the Following Functionality:*

+ *Format, Library, and Testing Improvements*
  o Enable processing of libraries with 1:many reagent:target assignments
  o Standardization and clarification of Annotation objects and symbol/identifier relationships
  o Implementation of factored quantile normalization for timecourse screens
  o Introduction of the `simpleResult` format and integration with associated functions
  o Conditional testing framework for quantifying and visualizing signal agreement between contrasts

+ *Transition to gene set enrichment testing via `Sparrow`*
  o Implement wrappers and provide recommendations fopr geneset enrichment testing in pooled screens
  o Implementation of GREAT-style pathway mapping for libraries with heterogenous target:gene mappings
  o Summarization tools for comparing enrichment signals across screens

+ *New Visualization and Interpretation Tools*
  o Signal Summary Barchart (Single or Multiple Contrasts)
  o Waterfall reagent/target/pathway visualization (Single Contrast)
  o Contrast comparison plots: 
    - Concordance at the Top (CAT)
    - Probability Space scatter plots
    - UpSet plots with conditional overlap framework
