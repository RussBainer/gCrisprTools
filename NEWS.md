### Version 0.99.0: 
Initial bioconductor submission. 

### Version 1.0.0: 
First official release to Bioconductor. 

### Version 1.1.0: 
Miscellaneous minor report edits and performance upgrades. Bugfix to P-value based scoring algorithm. Peter Haverty designated as official maintainer. 

### Version 1.15.0: 
Various code updates, new functions added for visualization, compatibility with modern R objects. Russell Bainer redesignated as official maintainer. 

### Version 1.20.0: 
Inclusion of factored quantile normalization and additional convenience functions, typo revisions in the vignettes. Suppression of results dataframe reordering by ct.generateResults().  Implementation of alternative reagent:target assignment, useful for situations where reagents are known to target multiple valid targets. 

### Version 2.0.0:
_*MAJOR UPDATE Adding the Following Functionality:*_
* *Format, Library, and Testing Improvements*
  o Enable processing of libraries with 1:many reagent:target assignments
  o Standardization and clarification of Annotation objects and symbol/identifier relationships
  o Implementation of factored quantile normalization for timecourse screens
  o Introduction of the `simpleResult` format and integration with associated functions
  o Conditional testing framework for quantifying and visualizing signal agreement between contrasts

* *Transition to gene set enrichment testing via `Sparrow`*

  o Implement wrappers and provide recommendations fopr geneset enrichment testing in pooled screens
  o Implementation of GREAT-style pathway mapping for libraries with heterogenous target:gene mappings
  o Summarization tools for comparing enrichment signals across screens

* *New Visualization and Interpretation Tools*

  o Signal Summary Barchart (Single or Multiple Contrasts)
  o Waterfall reagent/target/pathway visualization (Single Contrast)
  o Contrast comparison plots: 
    - Concordance at the Top (CAT)
    - Probability Space scatter plots
    - UpSet plots with conditional overlap framework
    
### Version 2.3.3 
Minor updates: 

- Modified ct.viewGuides to suppress reordering of provided gRNA signals. 
- Enabled use of SummarizedExperiment objects across all functions. 