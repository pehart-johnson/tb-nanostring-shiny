# TB NanoString Shiny App

This app was developed for Homework 3 in Machine Learning for Biomedical Data. It compares PCA, UMAP, and rank-2 NMF representations of gene-expression measurements from 179 samples and 107 genes.

Users can select the dimension-reduction method and the number of k-means clusters. The app displays samples colored by cluster and shaped by TB status, cluster centers, a cluster-versus-status table, and the adjusted Rand index.

## Files

- `app.R`: Complete Shiny application.
- `TBnanostring.rds`: Dataset required by the app.
- `README.md`: Project description and instructions.

## How to run

1. Download the repository using **Code → Download ZIP**, then extract the ZIP file.
2. Keep `app.R` and `TBnanostring.rds` together in the extracted folder.
3. Open RStudio. If the required packages are not installed, run:

```r
install.packages(c("shiny", "ggplot2", "uwot", "mclust", "BiocManager"))  # Install the main packages and Bioconductor installer.
BiocManager::install(c("Biobase", "NMF"), ask = FALSE, update = FALSE)  # Install NMF and its Bioconductor dependency.
```

4. Open `app.R` in RStudio and click **Run App**.
5. Wait for the initial dimension-reduction calculations to finish.

## Using the app

Choose PCA, UMAP, or NMF from the method menu and select a cluster count from 1 to 9. K-means operates on the two displayed dimensions. Changing the cluster count updates the clustering, table, and adjusted Rand index without recalculating the dimension reductions.

Higher adjusted Rand index values indicate closer agreement with TB status. This comparison describes the current dataset and does not measure prediction accuracy on new samples.
