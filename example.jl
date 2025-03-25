import Pkg
Pkg.activate(".") # Activate the project environment in the current directory
Pkg.instantiate() # Install dep pacakges

# Load packages
using CUDA
using CSV
using scLENS
using CairoMakie
using CodecLz4
CairoMakie.activate!(type="png")

# Environment configuration (Set the number of cores to use (for multi-core processing))
ENV["NUM_CORES"] = "15"
include("src/scICE.jl")

# Device selection for GPU or CPU processing
cur_dev = if CUDA.has_cuda()
    "gpu"
else
    "cpu"
end

# Load the compressed CSV file into a dataframe
# ndf = scLENS.read_file("data/Z8eq.csv.gz")
ndf = scLENS.read_file("data/ZhengMix_2410.jld2")
# Perform data preprocessing
pre_df = scLENS.preprocess(ndf)
# Create an embedding using scLENS
sclens_embedding = scLENS.sclens(pre_df,device_=cur_dev)

# Apply UMAP transformation
scLENS.apply_umap!(sclens_embedding)

# Visualize the embedding and save the UMAP distribution as an image
l_true = pre_df.cell
panel_0 = scLENS.plot_embedding(sclens_embedding,l_true)
# scLENS.plot_mpdist(sclens_embedding)
# scLENS.plot_stability(sclens_embedding)
# CSV.write("out/pca.csv",sclens_embedding[:pca_n1]) # Save the PCA results to a CSV file
# CSV.write("out/umap.csv",DataFrame(sclens_embedding[:umap],:auto)) # Save the UMAP results to a CSV file
save("out/umap_dist.png",panel_0)

# Apply scICE to the umap graph
clustering!(sclens_embedding)

# Visualize Inconsistency Coefficient (IC) and save it as an image
panel_1 = plot_ic(sclens_embedding)
save("out/ic_plot.png",panel_1)

# Retrieve the consistent cluster labels
label_out = get_rlabel!(sclens_embedding)
# Save the cluster labels to a CSV file
CSV.write("out/consistent_labels.csv",label_out)

# Set the number of clusters and visualize with cluster labels
n_clusters = 9
panel_2 = scLENS.plot_embedding(sclens_embedding,label_out[!,"l_$n_clusters"])
save("out/umap_dist_with_label$n_clusters.png",panel_2)

# Save output as an AnnData
scLENS.save_anndata("out/test.h5ad",sclens_embedding)
