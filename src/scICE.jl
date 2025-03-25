using PyCall
try
    pyimport("sys")
    println("Python is properly installed and accessible.")
    return true
catch e
    println("""
    Python is not properly installed or not accessible.
    Error: $e

    To install Python, follow the instructions below for your operating system:

    Windows:
      1. Download the Python installer from https://www.python.org/downloads/
      2. Run the installer and make sure to select 'Add Python to PATH'.
      3. After installation, restart your terminal or IDE.

    Linux (Debian/Ubuntu):
      1. Open a terminal.
      2. Run the following commands:
         sudo apt update
         sudo apt install python3

    MacOS:
      1. Open a terminal.
      2. Run the following command:
         brew install python
         (Note: You need Homebrew installed. Visit https://brew.sh/ if not installed)
    """)
    return false
end

function pip_installation(pkg_name::String)
    try
        pyimport(pkg_name)
        println("Python package '$pkg_name' is already installed.")
    catch
        println("Python package '$pkg_name' is not installed. Attempting to install.")
        try
            run(`$(PyCall.python) -m pip install $pkg_name`)
            println("Python package '$pkg_name' has been successfully installed.")
            pyimport(pkg_name)
        catch e
            println("Failed to install Python package '$pkg_name'.")
            println("Error: ", e)
            return nothing
        end
    end
end
pip_installation("igraph")
pip_installation("scipy")

using DataFrames
using ProgressMeter
using StatsBase
using Distributed
using SparseArrays
using LinearAlgebra: dot, I as L_I
using CairoMakie

n_cores = parse(Int, get(ENV, "NUM_CORES", "10")) 
addprocs(n_cores - nprocs() + 1)
@everywhere using SimpleWeightedGraphs: SimpleWeightedGraph
@everywhere using PyCall
@everywhere using StatsBase:mean
@everywhere const ig = pyimport("igraph")

function graph2ig(adj_m,weighted=true)
    sp_sparse = pyimport("scipy.sparse")
    a,b,c = findnz(adj_m)
    coo_mat = sp_sparse.coo_matrix((c,(a .-1,b .-1)),shape=[size(adj_m,1),size(adj_m,1)])
    if weighted
        return ig.Graph.Weighted_Adjacency(coo_mat,mode="undirected")
    else
        return ig.Graph.Adjacency(coo_mat,mode="undirected")
    end
end

@everywhere function clust_graph(igg;gamma=0.8,objective_function="CPM",n_iter=5,beta=0.1,init_mem=nothing)
    if igg.is_weighted()
        igl = igg.community_leiden(resolution_parameter=gamma,weights="weight",objective_function=objective_function,n_iterations=n_iter,beta=beta,initial_membership=init_mem)
        Vector{Int16}(igl.membership)
    else
        igl = igg.community_leiden(resolution_parameter=gamma,objective_function=objective_function,n_iterations=n_iter,beta=beta,initial_membership=init_mem)
        Vector{Int16}(igl.membership)
    end
end

function extract_arr(inp_arr)
    a_ = [Vector(s) for s in eachrow(inp_arr)]
    c_dict_ = countmap(a_)
    n_c = collect(values(c_dict_))
    arr_t = collect(keys(c_dict_))

    prob_arr = n_c ./ sum(n_c)
    sort_i = sortperm(prob_arr,rev=true)
    arr_ = arr_t[sort_i]
    n_c = n_c[sort_i]
    prob_arr2 = prob_arr[sort_i]
    
    prob_arr2 = prob_arr2 ./ sum(prob_arr2)
    arr_ = arr_

    return Dict(:arr => arr_, :parr => prob_arr2)
end

@everywhere function simmat_v2(inp_a::Vector{Int16}, inp_b::Vector{Int16}; d::Float64=0.9, flag::Bool=true)::Union{Float64, Vector{Float64}}
    n = length(inp_a)
    g_idx_a = unique(inp_a)
    g_idx_b = unique(inp_b)

    gg_idx_a = [Int[] for _ in g_idx_a]
    gg_idx_b = [Int[] for _ in g_idx_b]
    for i in 1:n
        id1_ = inp_a[i] + 1
        id2_ = inp_b[i] + 1
        push!(gg_idx_a[id1_], i)
        push!(gg_idx_b[id2_], i)
    end
    c_size1 = d./ length.(gg_idx_a)
    c_size2 = d./ length.(gg_idx_b)
   
    unique_ecs_vals = fill(NaN, length(g_idx_a), length(g_idx_b))
    ecs = zeros(n)
    ppr1 = zeros(n)
    ppr2 = zeros(n)

    for i in 1:n
        i1 = inp_a[i] + 1
        i2 = inp_b[i] + 1
        
        if isnan(unique_ecs_vals[i1, i2])
            # d_div_c_size1 = d / c_size1[i1]
            nei1 = gg_idx_a[i1]
            nei2 = gg_idx_b[i2]
            all_ = BitSet(vcat(nei1,nei2))
            
            for idx in nei1 
                @inbounds ppr1[idx] = c_size1[i1]
            end
            ppr1[i] = 1.0 - d + c_size1[i1]

            for idx in nei2
                @inbounds ppr2[idx] = c_size2[i2]
            end
            ppr2[i] = 1.0 - d + c_size2[i2]
            
            escore = 0
            for j in all_
                escore += abs(ppr2[j] - ppr1[j])
            end
            ecs[i] = escore

            for idx in nei1
                @inbounds ppr1[idx] = 0.0
            end
            for idx in nei2
                @inbounds ppr2[idx] = 0.0
            end
            unique_ecs_vals[i1, i2] = ecs[i]
        else
            ecs[i] = unique_ecs_vals[i1, i2]
        end
    end

    if flag
        return mean(1 .- 1 / (2 * d) .* ecs)
    else
        return 1 .- 1 / (2 * d) .* ecs
    end
    
end

function get_mei_from_array(inp,n_workers=nworkers())
    if length(inp[:arr]) == 1
        return ones(length(inp[:arr][1]))
    else
        wp = CachingPool(workers()[1:n_workers])
        tmp_S = let nu_mem=inp[:arr], p_mem=inp[:parr]
            hcat(vcat(pmap(i -> [simmat_v2(nu_mem[i],nu_mem[j],flag=false).*(p_mem[i]+p_mem[j]) for j=i+1:lastindex(nu_mem)], wp, 1:lastindex(nu_mem))...)...)
        end
        clear!(wp)
        return sum(tmp_S,dims=2)[:]./(length(inp[:parr])-1)    
    end
end

@everywhere function get_ic2(inp_c; n_workers=nworkers())
    nu_mem = hcat(inp_c[:arr]...)'
    wp = CachingPool(workers()[1:n_workers])
    tmp_S = let nu_mem=nu_mem
        pmap(i -> [simmat_v2(nu_mem[i,:],nu_mem[j,:]) for j=i+1:size(nu_mem,1)], wp, 1:size(nu_mem,1))
    end
    clear!(wp)
    prob_arr = if sum(inp_c[:parr]) != 1
        inp_c[:parr]./sum(inp_c[:parr])
    else
        inp_c[:parr]
    end

    S_ab = zeros(size(nu_mem,1),size(nu_mem,1))
    for i = 1:size(nu_mem,1)
        S_ab[i,i+1:end] = tmp_S[i]
    end
    S_ab = S_ab + S_ab' + L_I

    return nu_mem, S_ab, prob_arr, dot(S_ab*prob_arr,prob_arr)
end

function get_best_l(inp)
    nu_mem = hcat(inp[:arr]...)'
    wp = CachingPool(workers()[1:nworkers()])
    tmp_S = let nu_mem=nu_mem
            pmap(i -> [simmat_v2(nu_mem[i,:],nu_mem[j,:]) for j=i+1:size(nu_mem,1)], wp, 1:size(nu_mem,1))
    end
    clear!(wp)

    S_ab = zeros(size(nu_mem,1),size(nu_mem,1))
    for i = 1:size(nu_mem,1)
        S_ab[i,i+1:end] = tmp_S[i]
    end
    S_ab = S_ab + S_ab' + L_I

    nu_mem[argmax(sum(S_ab,dims=2)[:]),:]
end

function clustering!(a_dict,r=[1,20];n_workers=nworkers(),l_ground=nothing,n_steps=11,n_trials=15,n_boot=100,β=0.1,Nc=10,dN=2,max_iter=150,remove_threshold=1.15,g_type="umap",obj_fun="CPM",val_tol=1e-8)
    t_range=collect(r[1]:r[2])
    rmt_g,rigg,_ = if typeof(a_dict) <: Dict
        if g_type == "umap"
            SimpleWeightedGraph(a_dict[:umap_obj].graph), graph2ig(a_dict[:umap_obj].graph), Array(a_dict[:umap_obj].embedding'), a_dict[:umap_obj], a_dict[:umap_obj].graph
        elseif g_type == "snn"
            ig_ = graph2ig(adjacency_matrix(a_dict[:snn_obj][1]))
            if haskey(a_dict,:umap)
                a_dict[:snn_obj][1], ig_, a_dict[:umap]
            else
                a_dict[:snn_obj][1], ig_, nothing
            end
        elseif g_type == "knn"
            ig_ = graph2ig(adjacency_matrix(a_dict[:knn_obj][1]))
            if haskey(a_dict,:umap)
                a_dict[:knn_obj][1], ig_, a_dict[:umap]
            else
                a_dict[:knn_obj][1], ig_, nothing
            end
        elseif g_type == "harmony"
            if haskey(a_dict,:umap)
                a_dict[:harmony_graph][1], a_dict[:harmony_graph][2], a_dict[:harmony_umap]
            else
                a_dict[:harmony_graph][1], a_dict[:harmony_graph][2], nothing
            end
        elseif g_type == "testg"
            ig_ = graph2ig(adjacency_matrix(a_dict[:testg_obj][1]))
            if haskey(a_dict,:umap)
                a_dict[:testg_obj][1], ig_, a_dict[:umap]
            else
                a_dict[:testg_obj][1], ig_, nothing
            end
        end
    elseif typeof(a_dict[2]) <: PyObject && ((typeof(a_dict[1]) <: SimpleWeightedGraph) | (typeof(a_dict[1]) <: SimpleGraph))
        a_dict[1],a_dict[2], nothing
    end

    wp = CachingPool(workers()[1:n_workers])
    @everywhere rg_all = $rigg

    start_g, end_g = if obj_fun == "modularity"
        0, 10
    elseif obj_fun == "CPM"
        min(log(val_tol),-13), 0
    end

    g_dict = Dict{Int,Vector{Float64}}()
    left, right = start_g, end_g
    left_b, right_b = start_g, end_g
    a,b = start_g, end_g
    gam_xarr = Matrix{Float64}(undef,3,0)
    Nc_pre=3;N_cls=10;β_d=0.01;
    @showprogress "range_searching..." for i = 1:Int(ceil(length(t_range)/2))
        y = t_range[i]
        if haskey(g_dict,y)
            break
        end

        left, right = a, b
        flag_ = if obj_fun == "modularity"
            !isapprox(left,right,atol=val_tol)
        else obj_fun == "CPM"
            !isapprox(exp(left),exp(right),atol=val_tol)
        end
        while flag_
            mid = (left + right) / 2
            g_ = if obj_fun == "modularity"
                mid
            elseif obj_fun == "CPM"
                exp(mid)
            end

            X = let gam=g_
                Matrix(hcat(pmap(i -> clust_graph(rg_all,gamma=gam,objective_function=obj_fun,n_iter=Nc_pre,beta=β_d),wp,1:N_cls)...)')
            end
            ff_i = maximum(X,dims=2)[:] .+ 1
            nc_ = median(ff_i)
            if nc_ < y
                left = mid 
            else
                right = mid
            end
            gam_xarr = hcat(gam_xarr,[nc_,mid,get_ic2(extract_arr(X),n_workers=n_workers)[end]^-1])
            
            flag_ = if obj_fun == "modularity"
                !isapprox(left,right,atol=val_tol)
            else obj_fun == "CPM"
                !isapprox(exp(left),exp(right),atol=val_tol)
            end
        end
        
        left_b = copy(right)
        flag_ = if obj_fun == "modularity"
            !isapprox(left_b,right_b,atol=val_tol)
        else obj_fun == "CPM"
            !isapprox(exp(left_b),exp(right_b),atol=val_tol)
        end
        while flag_
            mid = (left_b + right_b) / 2
            g_ = if obj_fun == "modularity"
                mid
            elseif obj_fun == "CPM"
                exp(mid)
            end

            X = let gam=g_
                Matrix(hcat(pmap(i -> clust_graph(rg_all,gamma=gam,objective_function=obj_fun,n_iter=Nc_pre,beta=β_d),wp,1:N_cls)...)')
            end
            ff_i = maximum(X,dims=2)[:] .+ 1
            nc_ = median(ff_i)
            if nc_ > y
                right_b = mid
            else
                left_b = mid 
            end
            gam_xarr = hcat(gam_xarr,[nc_,mid,get_ic2(extract_arr(X),n_workers=n_workers)[end]^-1])

            flag_ = if obj_fun == "modularity"
                !isapprox(left,right,atol=val_tol)
            else obj_fun == "CPM"
                !isapprox(exp(left_b),exp(right_b),atol=val_tol)
            end
        end
        g_dict[y] = [left, right_b]

        left = sum(g_dict[y])/2
        right = b
        y = t_range[end-(i-1)]
        if haskey(g_dict,y)
            break
        end

        flag_ = if obj_fun == "modularity"
            !isapprox(left,right,atol=val_tol)
        else obj_fun == "CPM"
            !isapprox(exp(left),exp(right),atol=val_tol)
        end
        while flag_
            mid = (left + right) / 2
            g_ = if obj_fun == "modularity"
                mid
            elseif obj_fun == "CPM"
                exp(mid)
            end

            X = let gam=g_
                Matrix(hcat(pmap(i -> clust_graph(rg_all,gamma=gam,objective_function=obj_fun,n_iter=Nc_pre,beta=β_d),wp,1:N_cls)...)')
            end
            ff_i = maximum(X,dims=2)[:] .+ 1
            nc_ = median(ff_i)
            if nc_ < y
                left = mid 
            else
                right = mid
            end
            gam_xarr = hcat(gam_xarr,[nc_,mid,get_ic2(extract_arr(X),n_workers=n_workers)[end]^-1])
            flag_ = if obj_fun == "modularity"
                !isapprox(left,right,atol=val_tol)
            else obj_fun == "CPM"
                !isapprox(exp(left),exp(right),atol=val_tol)
            end
        end

        left_b = right
        right_b = b
        flag_ = if obj_fun == "modularity"
            !isapprox(left_b,right_b,atol=val_tol)
        else obj_fun == "CPM"
            !isapprox(exp(left_b),exp(right_b),atol=val_tol)
        end
        while flag_
            mid = (left_b + right_b) / 2
            g_ = if obj_fun == "modularity"
                mid
            elseif obj_fun == "CPM"
                exp(mid)
            end

            X = let gam=g_
                Matrix(hcat(pmap(i -> clust_graph(rg_all,gamma=gam,objective_function=obj_fun,n_iter=Nc_pre,beta=β_d),wp,1:N_cls)...)')
            end
            ff_i = maximum(X,dims=2)[:] .+ 1
            nc_ = median(ff_i)
            if nc_ > y
                right_b = mid
            else
                left_b = mid 
            end
            gam_xarr = hcat(gam_xarr,[nc_,mid,get_ic2(extract_arr(X),n_workers=n_workers)[end]^-1])
            flag_ = if obj_fun == "modularity"
                !isapprox(left,right,atol=val_tol)
            else obj_fun == "CPM"
                !isapprox(exp(left_b),exp(right_b),atol=val_tol)
            end
        end
        g_dict[y] = [left, right_b]
        
        left, = sum(g_dict[t_range[i]])/2
        right = sum(g_dict[t_range[end-(i-1)]])/2
        left_b = sum(g_dict[t_range[i]])/2
        right_b = sum(g_dict[t_range[end-(i-1)]])/2
        a, b = left, right_b
    end
    gam_xarr = gam_xarr[:,sortperm(gam_xarr[2,:])]
    gam_range = if obj_fun == "CPM" 
        Dict(s.first => (s.second[2] - s.second[1]) > 0 ? exp.(s.second) : exp.([g_dict[max(1,s.first-1)][2],
        g_dict[min(maximum(t_range),s.first+1)][1]]) for s in g_dict)
    elseif obj_fun == "modularity" 
        gam_range = Dict(s.first => (s.second[2] - s.second[1]) > 0 ? s.second : [g_dict[max(1,s.first-1)][2],
        g_dict[min(maximum(t_range),s.first+1)][1]] for s in g_dict)
    else
        println("Warning!")
        return nothing
    end
    
    tmp_dict = Dict(i=>Matrix{Float64}(undef,0,3) for i in t_range)
    for i = 1:lastindex(gam_xarr,2)
        g_ = if obj_fun == "CPM"
            exp(gam_xarr[2,i])
        elseif obj_fun =="modularity"
            gam_xarr[2,i]
        end
        n_ = [x.first for x in gam_range if x.second[1] <= g_ <= x.second[2]]
        if isempty(n_)
            continue
        end
        for n in n_
            if n in t_range
                tmp_dict[n] = vcat(tmp_dict[n],gam_xarr[:,i]')
            end
        end
    end
    filter!(x -> size(x.second,1) > 1 ? (minimum(x.second[:,3]) >= remove_threshold) : true, tmp_dict)

    ex_n = sort(collect(keys(tmp_dict)))
    println("$ex_n are removed from the candidates")

    list_l = Vector{Dict{Symbol, Vector}}([])
    list_l_best = Vector{Vector{Int16}}([])
    list_incons = Vector{Float64}([])
    list_allic = Vector{Vector{Float64}}([])
    list_gam = Vector{Float64}([])
    list_mn = Vector{Float64}([])
    list_k = Vector{Int}([])
    all_n = sort(intersect(collect(keys(gam_range)),setdiff(t_range,vcat(1,ex_n))))
    pbar = Progress(length(all_n),desc="zoomin_steps: ",showspeed=true)    
    tmp_get_ic = x -> x[end]^-1
    for i in all_n
        st_g, en_g = gam_range[i]
        d_g = if obj_fun == "modularity"
            (en_g-st_g)./n_steps
        elseif obj_fun == "CPM"
            (log(en_g)-log(st_g))./n_steps
        end
        
        g_x = try
            if obj_fun == "modularity"
                LinRange(round(st_g,digits=2),round(en_g,digits=2),n_steps)
                if st_g != en_g
                    st_g:d_g:en_g
                else
                    st_g-d_g : d_g /n_steps*2 : st_g + d_g
                end

            elseif obj_fun == "CPM"
                if st_g != en_g
                    exp.(log(st_g):d_g:log(en_g))
                else
                    exp.(log(st_g)-d_g : d_g /n_steps*2 : log(st_g) + d_g)
                end
            end
        catch
            []
        end
        mn_arr = zeros(size(g_x,1))

        X_list = Vector{Matrix}(undef,size(g_x,1))
        for (ii,g) in enumerate(g_x)
            X_list[ii] = let gam=g
                Matrix(hcat(pmap(i -> clust_graph(rg_all,gamma=gam,objective_function=obj_fun,n_iter=Nc,beta=β),wp,1:n_trials)...)')
            end
            mn_arr[ii] = median(maximum(X_list[ii],dims=2))+1
        end
        clear!(wp)
        
        b_idx = mn_arr .== i
        g_x = g_x[b_idx]
        X_list = X_list[b_idx]
        mn_arr = mn_arr[b_idx]
        
        e_arr = extract_arr.(X_list)
        b_ic = tmp_get_ic.(get_ic2.(e_arr,n_workers=n_workers))
        one_i = findfirst(b_ic .== 1)
        k = Nc

        if isempty(g_x)
            ProgressMeter.next!(pbar, showvalues=[(:searching,i)])
            continue
        elseif !isnothing(one_i)
            g_x = g_x[one_i]
            X_list = X_list[one_i]
            mn_arr = mn_arr[one_i]

            ic_mat = zeros(n_boot)
            for j = 1:n_boot
                smp_i = sample(1:size(X_list,1),size(X_list,1))
                X_arr = extract_arr(X_list[smp_i,:])
                ic_mat[j] = get_ic2(X_arr,n_workers=n_workers)[end]^-1
            end
            ic_median = median(ic_mat)
        
            X_earr = extract_arr(X_list)
            best_l = get_best_l(X_earr)
        
            push!(list_gam,g_x)
            push!(list_mn,mn_arr)
            push!(list_incons,ic_median)
            push!(list_allic,ic_mat)
            push!(list_l_best,best_l)
            push!(list_l,X_earr)
            push!(list_k,k)
            ProgressMeter.next!(pbar, showvalues=[(:searching,i),(:γ, g_x),(:nClusterings, size(X_list,1)),
            (:IC, median(ic_mat)), (:N_Iterations, k)])
            continue
        end

        nk = 1
        tank_bic_all = repeat(ones(length(b_ic)).*2,1,10)
        tank_bic = repeat(b_ic,10)
        ustable_id = ones(Bool,length(g_x))
        while true
            k += dN
            for i = 1:lastindex(g_x)
                X_list[i] = if ustable_id[i]
                    let gam=g_x[i], init_mem = X_list[i], dN = dN
                        Matrix(hcat(pmap(i -> clust_graph(rg_all,gamma=gam,objective_function=obj_fun,n_iter=dN,beta=β,init_mem=init_mem[i,:]),wp,1:n_trials)...)')
                    end
                else
                    X_list[i]
                end
            end
            clear!(wp)
            
            mn_arr = median.(maximum.(X_list,dims=2))[:].+1
            e_arr = extract_arr.(X_list)
            b_ic_t = tmp_get_ic.(get_ic2.(e_arr,n_workers=n_workers))
            
            tank_bic_all[:,1:end-1] = tank_bic_all[:,2:end]
            tank_bic_all[:,end] = b_ic_t
            last_bic = b_ic_t
            diff_bic = diff(tank_bic_all,dims=2)
            stable_idx = iszero.(diff_bic)
            ustable_id = (sum(stable_idx,dims=2)[:,1] .!= size(stable_idx,2)) .& (tank_bic_all[:,1] .- tank_bic_all[:,end] .>= 0)
            one_i = tank_bic_all[:,end] .== 1
            
            b_idx = (last_bic .<= quantile(last_bic,0.5)) .| ustable_id
            b_idx[argmin(last_bic)] = true
            if length(b_idx) == 1
                g_x = g_x[1]
                X_list = X_list[1]
                mn_arr = mn_arr[1]
                b_ic_t = b_ic_t[1]
                tank_bic = tank_bic_all[1,:]
                break
            elseif any(one_i)
                sel_i = findfirst(one_i)
                g_x = g_x[sel_i]
                X_list = X_list[sel_i]
                mn_arr = mn_arr[sel_i]
                b_ic_t = b_ic_t[sel_i]
                tank_bic = tank_bic_all[sel_i,:]
                break
            elseif all(stable_idx)
                sel_i = argmin(b_ic_t)
                g_x = g_x[sel_i]
                X_list = X_list[sel_i]
                mn_arr = mn_arr[sel_i]
                b_ic_t = b_ic_t[sel_i]
                tank_bic = tank_bic_all[sel_i,:]
                break
            elseif k >= max_iter
                sel_i = argmin(b_ic_t)
                g_x = g_x[sel_i]
                X_list = X_list[sel_i]
                mn_arr = mn_arr[sel_i]
                b_ic_t = b_ic_t[sel_i]
                tank_bic = tank_bic_all[sel_i,:]
                break
            elseif k > 100
                if median(b_ic_t) > 1.1
                    sel_i = argmin(b_ic_t)
                    g_x = g_x[sel_i]
                    X_list = X_list[sel_i]
                    mn_arr = mn_arr[sel_i]
                    b_ic_t = b_ic_t[sel_i]
                    tank_bic = tank_bic_all[sel_i,:]
                    break
                end
            else
                g_x = g_x[b_idx]
                X_list = X_list[b_idx]
                mn_arr = mn_arr[b_idx]
                b_ic_t = b_ic_t[b_idx]
                tank_bic_all = tank_bic_all[b_idx,:]
            end
            nk += 1
        end

        if typeof(g_x) <: Vector
            g_x = g_x[1]
            X_list = X_list[1]
            mn_arr = mn_arr[1]
        end
        
        ic_mat = zeros(n_boot)
        for j = 1:n_boot
            smp_i = sample(1:size(X_list,1),size(X_list,1))
            X_arr = extract_arr(X_list[smp_i,:])
            ic_mat[j] = get_ic2(X_arr,n_workers=n_workers)[end]^-1
        end
        ic_median = median(ic_mat)
       
        X_earr = extract_arr(X_list)
        best_l = get_best_l(X_earr)

        push!(list_gam,g_x)
        push!(list_mn,mn_arr)
        push!(list_incons,ic_median)
        push!(list_allic,ic_mat)
        push!(list_l_best,best_l)
        push!(list_l,X_earr)
        push!(list_k,k)
        ProgressMeter.next!(pbar, showvalues=[(:searching,i),(:γ, g_x),(:nClusterings, size(X_list,1)),
        (:IC, median(ic_mat)), (:N_Iterations, k),(:N_loops, nk)])
    end

    if isempty(findall(mean.(list_incons) .< 2))
        @everywhere begin
            rigg=nothing
            GC.gc()
        end
        return nothing
    else
        a_dict[:gamma] = list_gam
        a_dict[:labels] = list_l
        a_dict[:ic] = list_incons
        a_dict[:ic_vec] = list_allic
        a_dict[:n_cluster] = list_mn
        a_dict[:l_ground] = l_ground
        a_dict[:graph] = rmt_g
        a_dict[:best_l] = list_l_best
        a_dict[:n_iter] = list_k
        a_dict[:mei] = get_mei_from_array.(list_l,Ref(n_workers))

        @everywhere begin
            rigg=nothing
            GC.gc()
        end

    end
end


function get_rlabel!(input,th=1.005)
    b_idx = input[:ic] .< th
    labels_ = Dict(Int.(input[:n_cluster][b_idx]) .=> input[:best_l][b_idx])
    c_names = input[:pca].cell
    out_df = DataFrame(Dict(["l_"*string(s.first) => s.second for s in labels_]))
    insertcols!(out_df,1,:cell_id => c_names)
    input[:l_df] = out_df
    out_df
end

function plot_ic(out_,th=1.005;fig_size=(550,500))
    x = repeat(out_[:n_cluster],inner=length(out_[:ic_vec][1]))
    y = vcat(out_[:ic_vec]...)
    fig = Figure(size=fig_size)
    ax = Axis(fig[1, 1], xlabel="Number of clusters", ylabel="IC",
        limits=((0.5, maximum(out_[:n_cluster]) + 0.5), nothing))
    hlines!(ax, [th]; xmin = 0.02, xmax = 0.98, color=:red)
    boxplot!(ax, x, y)
    fig
end