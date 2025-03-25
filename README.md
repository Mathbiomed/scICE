# Single-Cell Analysis with scLENS and scICE

This project provides a pipeline for single-cell data analysis using the Julia packages `scLENS` and `scICE`. 
The `scICE` package, which stands for *Single Cell Inconsistency Clustering Estimator*, is specifically designed to perform multiple clustering runs and extract only the reliable labels that consistently appear across the runs.

It includes data preprocessing, embedding, clustering, and visualization of results. The code is designed to leverage CUDA for GPU acceleration when available.

## Requirements

To run this project, you will need the following:

- **Julia** (version 1.6 or higher recommended)
- **Python** (required for certain dependencies)

### GPU Requirements for CUDA
To use CUDA, you must have an NVIDIA GPU with CUDA capability, and the appropriate NVIDIA drivers must be installed on your system.

## Installation

1. **Clone the repository**:
   ```bash
   git clone TBA
   cd scICE

2. **Activate the project environment**:
   In the Julia REPL, navigate to the project directory and run:
   ```julia
   import Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```
   
## Usage (example.jl)

The following steps describe how to use this project to analyze single-cell data.

1. **Configure Processing Cores**:
   Set the number of CPU cores for processing:
   ```julia
   ENV["NUM_CORES"] = "12"
   ```

2. **Set up the environment**:
   Load the necessary packages and include the local scICE file:
   ```julia
   using CUDA, CSV, scLENS
   include("src/scICE.jl")
   using CairoMakie
   CairoMakie.activate!(type="png")
   ```

3. **Device Selection**:
   The device (CPU or GPU) is automatically selected based on CUDA availability:
   ```julia
   cur_dev = if CUDA.has_cuda()
       "gpu"
   else
       "cpu"
   end
   ```

4. **Data Preprocessing**:
   Load your single-cell data in compressed CSV format and preprocess it:
   ```julia
   ndf = scLENS.read_file(raw"data/Z8eq.csv.gz")
   pre_df = scLENS.preprocess(ndf)
   ```

5. **Embedding Creation**:
   Create an embedding for the preprocessed data using `scLENS`:
   ```julia
   sclens_embedding = scLENS.sclens(pre_df, device_=cur_dev)
   CSV.write("out/pca.csv", sclens_embedding[:pca_n1])
   ```

6. **UMAP Transformation**:
   Apply UMAP to the embedding and save the results:
   ```julia
   scLENS.apply_umap!(sclens_embedding)
   CSV.write("out/umap.csv", DataFrame(sclens_embedding[:umap], :auto))
   ```

7. **Visualization**:
   Plot the UMAP distribution and save the output:
   ```julia
   panel_0 = scLENS.plot_embedding(sclens_embedding)
   save("out/umap_dist.png",panel_0)
   ```

8. **Applying scICE**:
   Apply `scICE` clustering to the embedding:
   ```julia
   clustering!(sclens_embedding)
   ```

9. **Inconsistency Coefficient Visualization**:
   Visualize the Inconsistency coefficient and save it:
   ```julia
   panel_1 = plot_ic(sclens_embedding)
   save("out/ic_plot.png",panel_1)
   ```

10. **Consistent Cluster Label Extraction**:
    Extract consistent cluster labels and save them to a CSV file:
    ```julia
    label_out = get_rlabel!(sclens_embedding)
    CSV.write("out/consistent_labels.csv", label_out)
    ```

11. **Cluster Visualization**:
    Set the number of clusters and visualize them with labels:
    ```julia
    n_clusters = 9
    panel_2 = scLENS.plot_embedding(sclens_embedding, label_out[!, "l_$n_clusters"])
    save("out/umap_dist_with_label$n_clusters.png",panel_2)
    ```

12. **Save result as AnnData**:
   ```julia
   scLENS.save_anndata("out/test.h5ad",sclens_embedding)
   ```

## Output Files

The following output files will be generated:
- `out/pca.csv` - PCA results
- `out/umap.csv` - UMAP results
- `out/umap_dist.png` - UMAP distribution visualization
- `out/ic_plot.png` - Inconsistency coefficient plot
- `out/consistent_labels.csv` - Consistent cluster labels
- `out/umap_dist_with_label<n_clusters>.png` - UMAP distribution with cluster labels
