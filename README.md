# Single-Cell Analysis with scLENS and scICE

This project provides a pipeline for single-cell data analysis using the Julia packages `scLENS` and `scICE`.
The `scICE` package, which stands for *Single Cell Inconsistency Clustering Estimator*, is specifically designed to perform multiple clustering runs and extract only the reliable labels that consistently appear across the runs.

It includes data preprocessing, embedding, clustering, and visualization of results. The code is designed to leverage CUDA for GPU acceleration when available.

> [!CAUTION] **Important: Julia Version Compatibility** > There is a known **crash issue when using `PyCall` with Julia version 1.12.3**. To ensure a stable analysis environment, please **use a Julia version lower than 1.12.3** (e.g., 1.11.x).

## Requirements

To run this project, you will need the following:

-   **Julia**: Version 1.6 ~ 1.11.x is highly recommended. - Note: Avoid version 1.12.3 due to PyCall compatibility issues.
-   **Python**: Required for certain dependencies. Tested with Python 3.12.3.
-   **Operating Systems**: Tested on Windows 11 and Ubuntu 22.04.
-   **Key Julia Packages**: `CUDA`, `CSV`, `scLENS`, `CairoMakie` (These will be automatically installed via `Pkg.instantiate()` from the project environment).

### GPU Requirements for CUDA

-   **Hardware**: NVIDIA GPU with CUDA capability.
-   **Software**: Appropriate NVIDIA drivers must be installed.
-   **CUDA Toolkit**: Tested with CUDA Toolkit 12.2.

## Installation

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/Mathbiomed/scICE
    cd scICE
    ```

2.  **Activate the project environment**:
    In the Julia REPL, navigate to the project directory (`scICE`) and run:
    ```julia
    import Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    ```
    *Note: Installation, including package downloads and precompilation via `Pkg.instantiate()`, typically takes **5 to 15 minutes** depending on your system and network connection.*

## Usage (example.jl)

The following steps describe how to use this project to analyze single-cell data using the example script `example.jl`.

1.  **Configure Processing Cores**:
    Set the number of CPU cores for processing:
    ```julia
    ENV["NUM_CORES"] = "12"
    ```

2.  **Set up the environment**:
    Load the necessary packages and include the local scICE file:
    ```julia
    using CUDA, CSV, DataFrames, scLENS
    include("src/scICE.jl")
    using CairoMakie
    CairoMakie.activate!(type="png")
    ```

3.  **Device Selection**:
    The device (CPU or GPU) is automatically selected based on CUDA availability:
    ```julia
    cur_dev = if CUDA.has_cuda()
        "gpu"
    else
        "cpu"
    end
    ```

4.  **Data Preprocessing**:
    Load your single-cell data (example uses compressed CSV) and preprocess it:
    ```julia
    ndf = scLENS.read_file(raw"data/Z8eq.csv.gz")
    pre_df = scLENS.preprocess(ndf)
    ```

5.  **Embedding Creation**:
    Create an embedding for the preprocessed data using `scLENS`:
    ```julia
    sclens_embedding = scLENS.sclens(pre_df, device_=cur_dev)
    CSV.write("out/pca.csv", sclens_embedding[:pca_n1])
    ```

6.  **UMAP Transformation**:
    Apply UMAP to the embedding and save the results:
    ```julia
    scLENS.apply_umap!(sclens_embedding)
    CSV.write("out/umap.csv", DataFrame(sclens_embedding[:umap], :auto))
    ```

7.  **Visualization**:
    Plot the UMAP distribution and save the output:
    ```julia
    panel_0 = scLENS.plot_embedding(sclens_embedding)
    save("out/umap_dist.png",panel_0)
    ```

8. **Applying scICE**: Apply `scICE` clustering to the embedding:
    ```julia
    clustering!(sclens_embedding)
    ```
    
    By default, `scICE` explores cluster numbers ranging from 1 to 20 (this is the default value for the optional second argument `r`, as seen in the function signature `clustering!(a_dict, r=[1,20]; ...)`). If you wish to focus the analysis on a specific range of cluster numbers, for instance, from 5 to 10 clusters, you provide this range as the second argument:
    ```julia
    clustering!(sclens_embedding, [5,10])
    ```
    
    This enables you to find consistent cluster labels more efficiently within an anticipated range.

9.  **Inconsistency Coefficient Visualization**:
    Visualize the Inconsistency coefficient and save it:
    ```julia
    panel_1 = plot_ic(sclens_embedding)
    save("out/ic_plot.png",panel_1)
    ```


10. **Consistent Cluster Label Extraction**: Extract consistent cluster labels using `get_rlabel!` and save them to a CSV file. This function filters labels based on an Inconsistency Coefficient (IC) threshold.
    ```julia
    label_out = get_rlabel!(sclens_embedding)
    CSV.write("out/consistent_labels.csv", label_out)
    ```
    
    The IC threshold parameter (`th`) defaults to `1.005`. This value is passed as the optional second argument to `get_rlabel!` and can be adjusted if needed. For example, to change the threshold to `1.01`, you would call the function like this:
    ```julia
    label_out = get_rlabel!(sclens_embedding, 1.01)
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

### Expected Runtime (Example Data)

Running the `example.jl` script with the provided sample data with ~10,000 cells:

-   **scLENS embedding (`sclens` function):** Approximately **2-5 minutes**.
-   **scICE clustering (`clustering!` function):** Approximately **10-15 minutes**.

*Note: Runtimes can vary significantly based on your hardware (CPU/GPU specifics, RAM), the number of cores configured, and the size/complexity of the input data.*

## Output Files

Running the `example.jl` script will generate the following files in the `out/` directory:

-   `pca.csv`: PCA results.
-   `umap.csv`: UMAP coordinates.
-   `umap_dist.png`: Visualization of the UMAP embedding.
-   `ic_plot.png`: Plot of the inconsistency coefficient.
-   `consistent_labels.csv`: Consistent cluster labels generated by scICE.
-   `umap_dist_with_label<n_clusters>.png`: UMAP embedding colored by consistent cluster labels (e.g., `umap_dist_with_label9.png`).
-   `test.h5ad`: Output saved in AnnData format for compatibility with Python tools.

