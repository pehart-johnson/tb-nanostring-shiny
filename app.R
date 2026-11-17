# Load the packages used to build the interface and plot.
library(shiny)  # Provide the functions used to build and run the app.
library(ggplot2)  # Provide the plotting functions.

# Load the dataset directly so app.R can run independently.
tb_data <- readRDS("TBnanostring.rds")  # Read the dataset from the app's folder.
tb_status <- factor(tb_data$TB_Status, levels = c("TB", "LTBI"))  # Store the status labels in sample order.
gene_matrix <- as.matrix(tb_data[, names(tb_data) != "TB_Status"])  # Keep only the numeric gene measurements.

# Calculate PCA once, before starting the server.
pca_model <- prcomp(  # Fit PCA to the gene measurements.
  gene_matrix,  # Use samples as rows and genes as columns.
  center = TRUE,  # Subtract each gene's mean.
  scale. = FALSE  # Keep the existing scaling; these genes already have unit variance.
)  # Finish fitting PCA.

# Calculate UMAP once using the settings from Task 2.1.
set.seed(123)  # Make the UMAP calculation reproducible.
umap_coordinates <- uwot::umap(  # Calculate two UMAP coordinates for each sample.
  X = gene_matrix,  # Supply the gene measurements.
  n_neighbors = 15,  # Use a neighborhood size of 15.
  n_components = 2,  # Produce two coordinates per sample.
  min_dist = 0.1,  # Set how closely nearby points can be packed.
  metric = "euclidean",  # Measure distances using the gene measurements.
  scale = FALSE,  # Keep the existing gene scaling.
  n_threads = 1,  # Use one thread for supporting calculations.
  n_sgd_threads = 1,  # Use one thread when arranging points for reproducibility.
  verbose = FALSE  # Hide routine progress messages.
)  # Finish calculating UMAP.

# Prepare non-negative measurements for NMF.
gene_minimums <- apply(gene_matrix, 2, min)  # Find the minimum of each gene column.
nmf_input <- sweep(gene_matrix, 2, gene_minimums, FUN = "-")  # Subtract each column's minimum.

# Calculate NMF once using the settings from Question 3.
nmf_model <- NMF::nmf(  # Fit a non-negative matrix factorization.
  x = nmf_input,  # Supply the shifted gene measurements.
  rank = 2,  # Represent each sample using two components.
  method = "brunet",  # Use the Brunet NMF algorithm.
  nrun = 10,  # Try ten initializations and retain the best fit.
  seed = 123,  # Set the seed for reproducibility.
  .options = list(parallel = FALSE)  # Run the fits sequentially.
)  # Finish fitting NMF.

# Store the three sets of sample coordinates for reuse.
embeddings <- list(  # Create a named collection of reduced datasets.
  PCA = pca_model$x[, 1:2, drop = FALSE],  # Keep PC1 and PC2.
  UMAP = umap_coordinates,  # Keep the two UMAP coordinates.
  NMF = NMF::basis(nmf_model)  # Keep the two columns of W, with one row per sample.
)  # Finish creating the collection.

# Build the interface.
ui <- fluidPage(  # Create a page that adjusts to the window size.
  titlePanel("TB gene expression: dimension reduction and k-means"),  # Add the app title.
  sidebarLayout(  # Place controls and results beside the plot.
    sidebarPanel(  # Create the sidebar.
      width = 3,  # Give the sidebar three of the twelve available layout columns.
      selectInput(  # Add a menu for choosing the dimension-reduction method.
        inputId = "method",  # Name the value that the server will read.
        label = "Dimension-reduction method",  # Label the menu.
        choices = c("PCA", "UMAP", "NMF"),  # Provide the three required methods.
        selected = "PCA"  # Display PCA when the app starts.
      ),  # Finish creating the method menu.
      numericInput(  # Add a control for the number of clusters.
        inputId = "clusters",  # Name the value that the server will read.
        label = "Number of clusters (k)",  # Label the control.
        value = 2,  # Start with two clusters.
        min = 1,  # Allow at least one cluster.
        max = 9,  # Match the tutorial's upper limit.
        step = 1  # Change the value in whole-number increments.
      ),  # Finish creating the cluster control.
      h4("Clusters vs. TB status"),  # Add a heading above the comparison table.
      tableOutput("cluster_table"),  # Reserve space for the cluster counts.
      h4("Adjusted Rand index"),  # Add a heading above the agreement score.
      textOutput("ari")  # Reserve space for the agreement score.
    ),  # Finish creating the sidebar.
    mainPanel(  # Create the main plotting area.
      width = 9,  # Give the plot the remaining nine layout columns.
      plotOutput("plot1", height = "550px")  # Reserve space for the scatterplot.
    )  # Finish creating the main panel.
  )  # Finish arranging the sidebar and main panel.
)  # Finish defining the interface.

# Define how the app responds to the controls.
server <- function(input, output, session) {  # Begin the server function.
  
  selectedData <- reactive({  # Update the selected coordinates when the method changes.
    coordinates <- embeddings[[input$method]]  # Retrieve the previously calculated coordinates.
    data.frame(  # Give all three methods consistent column names.
      Dimension1 = coordinates[, 1],  # Store the first coordinate for each sample.
      Dimension2 = coordinates[, 2]  # Store the second coordinate for each sample.
    )  # Finish creating the selected dataset.
  })  # Finish defining the reactive dataset.
  
  clusters <- reactive({  # Recalculate clustering when the method or k changes.
    req(input$clusters)  # Wait until a cluster count is available.
    k <- input$clusters  # Store the requested number of clusters.
    validate(  # Check that the requested number is valid.
      need(k >= 1 && k <= 9 && k == floor(k), "Enter a whole number from 1 to 9.")  # Reject invalid values.
    )  # Finish checking the input.
    set.seed(123)  # Make clustering reproducible for the selected method and k.
    kmeans(  # Divide samples into clusters using the displayed coordinates.
      x = selectedData(),  # Use the two reduced dimensions.
      centers = k,  # Request the selected number of clusters.
      nstart = 25,  # Try 25 starting configurations.
      iter.max = 100  # Allow up to 100 iterations per configuration.
    )  # Return the fitted clustering result.
  })  # Finish defining the reactive clustering.
  
  output$plot1 <- renderPlot({  # Create and update the scatterplot.
    fit <- clusters()  # Retrieve the shared clustering result.
    plot_data <- data.frame(  # Combine coordinates, cluster assignments, and status labels.
      selectedData(),  # Include the two selected coordinates.
      Cluster = factor(fit$cluster),  # Treat cluster numbers as categories for coloring.
      TB_Status = tb_status  # Include the known status for choosing point shapes.
    )  # Finish creating the plotting dataset.
    centers <- as.data.frame(fit$centers)  # Extract the cluster-center coordinates.
    axis_labels <- switch(  # Choose axis labels appropriate for the selected method.
      input$method,  # Read the selected method.
      PCA = c("PC1 score", "PC2 score"),  # Label PCA scores.
      UMAP = c("UMAP 1", "UMAP 2"),  # Label UMAP coordinates.
      NMF = c("NMF component 1 (W1)", "NMF component 2 (W2)")  # Label NMF weights.
    )  # Finish choosing the axis labels.
    
    ggplot(  # Begin constructing the scatterplot.
      plot_data,  # Supply the sample coordinates and labels.
      aes(x = Dimension1, y = Dimension2, color = Cluster, shape = TB_Status)  # Map position, color, and shape.
    ) +  # Finish the initial setup and add plot layers.
      geom_point(size = 2.5, alpha = 0.85) +  # Draw one point per sample.
      geom_point(  # Add a separate layer for cluster centers.
        data = centers,  # Supply the center coordinates.
        aes(x = Dimension1, y = Dimension2),  # Position each center.
        inherit.aes = FALSE,  # Avoid applying the sample-specific color and shape mappings.
        shape = 4, color = "black", size = 6, stroke = 1.5  # Draw centers as large black X marks.
      ) +  # Finish the center layer.
      scale_shape_manual(values = c(TB = 16, LTBI = 17)) +  # Show TB as circles and LTBI as triangles.
      coord_equal() +  # Give equal coordinate distances equal physical lengths.
      labs(  # Add the plot title, subtitle, and labels.
        title = paste(input$method, "with k-means clustering"),  # Identify the selected method.
        subtitle = paste("k =", input$clusters, "| Black X marks cluster centers"),  # Explain the centers.
        x = axis_labels[1],  # Label the horizontal axis.
        y = axis_labels[2],  # Label the vertical axis.
        color = "Cluster",  # Label the color legend.
        shape = "TB status"  # Label the shape legend.
      ) +  # Finish adding labels.
      theme_minimal(base_size = 13) +  # Use a simple background and readable text.
      theme(legend.position = "bottom")  # Place both legends below the plot.
  })  # Finish defining the plot output.
  
  output$cluster_table <- renderTable({  # Create and update the comparison table.
    counts <- table(  # Count samples in each cluster and status combination.
      Cluster = clusters()$cluster,  # Use the shared cluster assignments.
      TB_Status = tb_status  # Compare them with the known TB status.
    )  # Finish counting samples.
    data.frame(  # Prepare the counts for display.
      Cluster = rownames(counts),  # Include the cluster numbers as a visible column.
      as.data.frame.matrix(counts),  # Add separate count columns for TB and LTBI.
      row.names = NULL  # Remove automatic row labels.
    )  # Return the display table.
  }, striped = TRUE, bordered = TRUE, digits = 0)  # Format the table with integer counts.
  
  output$ari <- renderText({  # Create and update the agreement score.
    score <- mclust::adjustedRandIndex(  # Calculate agreement adjusted for chance.
      clusters()$cluster,  # Supply the shared cluster assignments.
      as.character(tb_status)  # Supply the known status labels.
    )  # Finish calculating the adjusted Rand index.
    sprintf("%.4f", score)  # Display the score with four decimal places.
  })  # Finish defining the agreement-score output.
  
}  # Finish defining the server.

shinyApp(ui = ui, server = server)  # Connect the interface and server to launch the app.