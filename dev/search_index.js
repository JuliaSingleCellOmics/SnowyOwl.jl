var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SnowyOwl","category":"page"},{"location":"#SnowyOwl","page":"Home","title":"SnowyOwl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SnowyOwl]","category":"page"},{"location":"","page":"Home","title":"Home","text":"SnowyOwl.Preprocess.quality_control_metrics\nSnowyOwl.Preprocess.quality_control_metrics!","category":"page"},{"location":"#SnowyOwl.Preprocess.quality_control_metrics","page":"Home","title":"SnowyOwl.Preprocess.quality_control_metrics","text":"quality_control_metrics(p; kwargs...)\nquality_control_metrics(X, obs, var; obsindex, varindex, kwargs...)\n\nCalculate quality control metrics.\n\nArguments\n\np::Profile: The profile object to calculate on.\nX::AbstractMatrix: The count matrix to calculate on.\nobs::DataFrame: The feature information matrix with cell information.\nvar::DataFrame: The feature information matrix with gene information.\n\nCommon Keyword arguments\n\nobsindex::Symbol=:barcode: The index of dataframe obs.\nvarindex::Symbol=:gene_symbols: The index of dataframe var.\nuse_log1p::Bool: Computing log1p-transformed quality control metric.\nqc_vars::AbstractVector{String}=[\"mt\"]: The boolean mask for identifying variables you could calculate on.\npercent_top=nothing: The proportions of top genes to calculate on.\n\nExample\n\nSee also quality_control_metrics! for non-inplace operation.\n\n\n\n\n\n","category":"function"},{"location":"#SnowyOwl.Preprocess.quality_control_metrics!","page":"Home","title":"SnowyOwl.Preprocess.quality_control_metrics!","text":"quality_control_metrics!(p; kwargs...)\nquality_control_metrics!(X, obs, var; kwargs...)\n\nCalculate quality control metrics.\n\nArguments\n\np::AnnotatedProfile: The profile object to calculate on.\nX::AbstractMatrix: The count matrix to calculate on.\nobs::DataFrame: The feature information matrix with cell information.\nvar::DataFrame: The feature information matrix with gene information.\n\nCommon Keyword arguments\n\nuse_log1p::Bool: Computing log1p-transformed quality control metric.\nqc_vars::AbstractVector{String}=[\"mt\"]: The boolean mask for identifying variables you could calculate on.\npercent_top=nothing: The proportions of top genes to calculate on.\n\nExample\n\nSee also quality_control_metrics for non-inplace operation.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"SnowyOwl.Preprocess.filter_genes!\nSnowyOwl.Preprocess.filter_cells!","category":"page"},{"location":"#SnowyOwl.Preprocess.filter_genes!","page":"Home","title":"SnowyOwl.Preprocess.filter_genes!","text":"Filter out genes which are not satisfy specific criteria.\n\n\n\n\n\n","category":"function"},{"location":"#SnowyOwl.Preprocess.filter_cells!","page":"Home","title":"SnowyOwl.Preprocess.filter_cells!","text":"Filter out cells which are not satisfy specific criteria.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"SnowyOwl.Preprocess.logarithmize\nSnowyOwl.Preprocess.logarithmize!","category":"page"},{"location":"#SnowyOwl.Preprocess.logarithmize","page":"Home","title":"SnowyOwl.Preprocess.logarithmize","text":"logarithmize(prof; omicsname=:RNA, layer=:count, kwargs...)\nlogarithmize(X; kwargs...)\n\nLogarithmize the data matrix.\n\nArguments\n\nprof::AnnotatedProfile: The profile object to calculate on.\nX::AbstractMatrix: The count matrix to calculate on.\n\nSpecific keyword arguments\n\nomicsname::Symbol: The OmicsProfile specified to calculate on.\nlayer::Symbol: The layer specified to calculate on.\n\nCommon keyword arguments\n\nbase::Real=ℯ: Base of the logarithm.\n\nExamples\n\nSee also logarithmize! for inplace operation.\n\n\n\n\n\n","category":"function"},{"location":"#SnowyOwl.Preprocess.logarithmize!","page":"Home","title":"SnowyOwl.Preprocess.logarithmize!","text":"logarithmize!(prof; omicsname=:RNA, layer=:count, kwargs...)\nlogarithmize!(X; kwargs...)\n\nLogarithmize the data matrix.\n\nArguments\n\nprof::AnnotatedProfile: The profile object to calculate on.\nX::AbstractMatrix: The count matrix to calculate on.\n\nSpecific keyword arguments\n\nomicsname::Symbol: The OmicsProfile specified to calculate on.\nlayer::Symbol: The layer specified to calculate on.\n\nCommon keyword arguments\n\nbase::Real=ℯ: Base of the logarithm.\n\nExamples\n\nSee also logarithmize for non-inplace operation.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"SnowyOwl.Preprocess.highly_variable_genes\nSnowyOwl.Preprocess.highly_variable_genes!","category":"page"},{"location":"#SnowyOwl.Preprocess.highly_variable_genes","page":"Home","title":"SnowyOwl.Preprocess.highly_variable_genes","text":"highly_variable_genes(prof, method; omicsname=:RNA, layer=:count, kwargs...)\nhighly_variable_genes(X, var, method; varname=:gene_symbols, kwargs...)\n\nCalculate highly variable genes and return a new DataFrame with column :highlyvariable, which selects highly variable genes from given gene set. Additional information including :means, :dispersions and `:dispersionsnormwill be added to returnedDataFrame`.\n\nArguments\n\nprof::AnnotatedProfile: The profile object to calculate on.\nX::AbstractMatrix: The count matrix to calculate on.\nvar::DataFrame: The feature information matrix with gene information.\nmethod::Symbol: Method to calculate highly variable genes, available for :cellranger,   :seurat and :seuratv3.\n\nSpecific keyword arguments\n\nomicsname::Symbol: The OmicsProfile specified to calculate on.\nlayer::Symbol: The layer specified to calculate on.\nvarname::Symbol: The variable name to be specified as identifier for genes.\n\nCommon keyword arguments\n\nntop_genes::Int=-1: Number of top variable genes to be selected. Specify -1 to switch   to selection by mean and dispersion. Available for :cellranger and :seurat methods.\nmin_disp::Real=0.5: Minimum dispersion for selecting highly variable genes. Available   for :cellranger and :seurat methods.\nmax_disp::Real=Inf: Maximum dispersion for selecting highly variable genes. Available   for :cellranger and :seurat methods.\nmin_mean::Real=0.0125: Minimum mean for selecting highly variable genes. Available for   :cellranger and :seurat methods.\nmax_mean::Real=3.: Maximum mean for selecting highly variable genes. Available for   :cellranger and :seurat methods.\n\nExamples\n\nSee also highly_variable_genes! for inplace operation.\n\n\n\n\n\n","category":"function"},{"location":"#SnowyOwl.Preprocess.highly_variable_genes!","page":"Home","title":"SnowyOwl.Preprocess.highly_variable_genes!","text":"highly_variable_genes!(prof, method; omicsname=:RNA, layer=:count, kwargs...)\nhighly_variable_genes!(X, var, method; kwargs...)\n\nCalculate highly variable genes and modify var directly by adding column :highlyvariable, which selects highly variable genes from given gene set. Additional information including :means, :dispersions and `:dispersionsnormwill be added tovar`.\n\nArguments\n\nprof::AnnotatedProfile: The profile object to calculate on.\nX::AbstractMatrix: The count matrix to calculate on.\nvar::DataFrame: The feature information matrix with gene information.\nmethod: Method to calculate highly variable genes, available for :cellranger,   :seurat and :seuratv3.\n\nSpecific keyword arguments\n\nomicsname::Symbol: The OmicsProfile specified to calculate on.\nlayer::Symbol: The layer specified to calculate on.\n\nCommon keyword arguments\n\nntop_genes::Int=-1: Number of top variable genes to be selected. Specify -1 to switch   to selection by mean and dispersion. Available for :cellranger and :seurat methods.\nmin_disp::Real=0.5: Minimum dispersion for selecting highly variable genes. Available   for :cellranger and :seurat methods.\nmax_disp::Real=Inf: Maximum dispersion for selecting highly variable genes. Available   for :cellranger and :seurat methods.\nmin_mean::Real=0.0125: Minimum mean for selecting highly variable genes. Available for   :cellranger and :seurat methods.\nmax_mean::Real=3.: Maximum mean for selecting highly variable genes. Available for   :cellranger and :seurat methods.\n\nExamples\n\nSee also highly_variable_genes for non-inplace operation.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"SnowyOwl.Analysis.project\nSnowyOwl.Analysis.project!","category":"page"},{"location":"#SnowyOwl.Analysis.project","page":"Home","title":"SnowyOwl.Analysis.project","text":"project(prof, method; omicsname=:RNA, layer=:count, kwargs...)\nproject(X, method; kwargs...)\n\nProject the count matrix from original space to low-dimensional space.\n\nArguments\n\nprof::AnnotatedProfile: The profile object to calculate on.\nX::AbstractMatrix: The count matrix to calculate on. The dimension of X should be   (nfeature, nsample) where nfeature is number of features and nsample is number   of samples.\nmethod::AnalysisMethod: Method to project count matrix, available for PCAMethod() and   UMAPMethod().\n\nSpecific keyword arguments\n\nomicsname::Symbol: The OmicsProfile specified to calculate on.\nlayer::Symbol: The layer specified to calculate on.\n\nCommon keyword arguments\n\ndims::Int=50: The dimension of target space. Should lower or equal to the dimension of   feature from count matrix which is dims <= nfeature.\nn_neighbors::Int=20: The number of local neighborhood points used in approximations of   manifold structure. It is only effective when method = UMAPMethod().\nmin_dist::Real=0.5: This controls how tightly the embedding is allowed compress points   together. Larger values prevent points from packing together and will prone to the   preservation of global structure instead. Smaller values result in a more clustered/clumped   embedding and prone to the preservation of local structure. It is only effective when   method = UMAPMethod().\n\nExamples\n\nSee also project! for inplace operation.\n\n\n\n\n\n","category":"function"},{"location":"#SnowyOwl.Analysis.project!","page":"Home","title":"SnowyOwl.Analysis.project!","text":"project!(prof, method; omicsname=:RNA, layer=:count, kwargs...)\nproject!(X, method; kwargs...)\n\nProject the data matrix from original space to (usually) low-dimensional space.\n\nArguments\n\nprof::AnnotatedProfile: The profile object to calculate on.\nX::AbstractMatrix: The count matrix to calculate on. The dimension of X should be   (nfeature, nsample) where nfeature is number of features and nsample is number   of samples.\nmethod::AnalysisMethod: Method to project count matrix, available for PCAMethod() and   UMAPMethod().\n\nSpecific keyword arguments\n\nomicsname::Symbol: The OmicsProfile specified to calculate on.\nlayer::Symbol: The layer specified to calculate on.\n\nCommon keyword arguments\n\ndims::Int=50: The dimension of target space. Should lower or equal to the dimension of   feature from count matrix which is dims <= nfeature.\nn_neighbors::Int=20: The number of local neighborhood points used in approximations of   manifold structure. It is only effective when method = UMAPMethod().\nmin_dist::Real=0.5: This controls how tightly the embedding is allowed compress points   together. Larger values prevent points from packing together and will prone to the   preservation of global structure instead. Smaller values result in a more clustered/clumped   embedding and prone to the preservation of local structure. It is only effective when   method = UMAPMethod().\n\nExamples\n\nSee also project for non-inplace operation.\n\n\n\n\n\n","category":"function"}]
}
