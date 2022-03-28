%% Q0: path configurations

% Read files to work
expr = sc_readmtxfile("../data/matrix.mtx");
meta = readtable("../data/metadata.csv");

barcode = tdfread("../data/barcodes.tsv").barcode;
barcode = cellstr(barcode);

gene_list = tdfread("../data/features.tsv").genes;
gene_list = cellstr(gene_list);


%% Q1: Dataset QC

% 1.1
dim = size(expr);
disp(['# genes: ', num2str(dim(1)), '; # cells: ', num2str(dim(2))]);

% 1.2
disp(expr(1:10, 1:10)); 
expr_norm = sc_norm(expr);

% 1.3
[libsize, sortIdx] = sort(sum(expr, 1), 'descend');
libsize = log10(libsize);  % take log on libsize
barcode = barcode(sortIdx);


% --plot--
fig1 = figure;
idxs = linspace(0, length(libsize), length(libsize));
grid on;
plot(log10(idxs), libsize);
xlabel("Sorted barcode index (descending)");
ylabel("log10 library size");
title("log10 total # transcripts vs. Barcode");


fig2 = figure;
bar(libsize);
grid on;
xticklabels(barcode); % todo: convert barcode to single cell
fig2.WindowState = 'maximized';
xlabel("Sorted barcode (descending)");
ylabel("log10 library size");
title("log10 total # transcripts vs. Barcode")


% 1.4
sample_lbls = unique(meta.sample);
n_samples = length(sample_lbls);
disp(['# unique samples: ', num2str(n_samples)]);
libsize = sum(expr, 1);  % order of original barcode (before sorting)

% median # transcript & nonzero genes per sample
median_transcript = zeros(n_samples, 1);
nonzero_gexp = zeros(n_samples, 1);

for i=1:n_samples
    idx = meta.sample == string(sample_lbls(i));
    median_transcript(i) = median(libsize(idx));
    nonzero_gexp(i) = median(sum(expr(:,idx) > 0, 1));
end
    

%% Q2: UMAP
% 2.1
mean_gexp = mean(expr_norm, 2);
std_gexp = std(expr_norm, 0, 2);
dispersion = (std_gexp ./ mean_gexp).^2;

plot(log10(mean_gexp), log10(dispersion), '.')
grid on;
xlabel("Mean gexp (logscale)");
ylabel("Dispersion gexp (logscale)");
title("Dispersion vs. Mean gexp");

% 2.2

% choose top 2000 hvgs
[T] = sc_hvg(expr_norm, gene_list, true, true);
hvgs = T.genes(1:2000);
idxs = zeros(length(gene_list), 1);

for i=1:length(gene_list)
    if ( find( string(hvgs) == string(gene_list(i))) > 0 )
        idxs(i) = 1;
    end
end

expr_hvgs = expr_norm(find(idxs), :);


% umap: we color the labels by sample source
%[s] = sc_umap(expr_hvgs, 2, false, true);

fig3 = figure;
gscatter(s(:,1), s(:,2), meta.sample, '', '', 8);
xlabel("Dim 1");
ylabel("Dim 2");
title("umap - 2.2");


% 2.3
s_list = {};

dists = linspace(0.1, 0.625, 4);
n_nbrs = linspace(5, 45, 5);

for i=1:length(dists)
    for j=1:length(n_nbrs) 
        dist = dists(i);
        n_nbr = n_nbrs(j);
        [curr_s] = sc_umap(expr_hvgs, 3, true, true, dist, n_nbr);
        s_list{end+1} = curr_s;
    end
end

fig4 = figure;
idx = 1;
for i=1:length(dists)
    for j=1:length(n_nbrs)
        dist = dists(i);
        n_nbr = n_nbrs(j);
        
        subplot(4, 5, idx);
        % 2D scatter (first 2 axes)
        %gscatter(s_list{idx}(:,1), s_list{idx}(:,2), meta.sample, '', '.', 4, 'off');
        % 3D scatter 
        gscatter3b(s_list{idx}(:,1), s_list{idx}(:,2), s_list{idx}(:,3), meta.sample, jet(5), '.', 4, 'off');
        title(sprintf('dist: %.2f; nbrs: %d', dist, n_nbr));
        
        idx = idx+1;
    end
end

clear idx;

%2.4
% (1).
marker_genes = ["ITGAM", "CTLA4", "TNFRSF18", "MOG"];
marker_idxs = zeros(4, 1);

for i=1:length(marker_genes)
    marker_idxs(i) = find( string(gene_list) == marker_genes(i) );
end

fig5 = figure;
for i=1:length(idxs)
    subplot(1,4,i);
    gscatter(s(:,1), s(:,2), expr_norm(idxs(i), :), viridis_white, '.', 4, 'off');
    title(["Gene expression - " marker_genes(i)]);
end

fig5.Position = [390, 317, 2500, 420];
saveas(fig5, "results/2_4_1.png");


% (3).
cell_type = meta.cell_assignment;
unique_cts = unique(cell_type);

mp_idxs = find(cellstr(cell_type) == string(unique_cts(1)));  % maacrophage
tumor_idxs = find(cellstr(cell_type) == string(unique_cts(2))); % malignant
ol_idxs = find(cellstr(cell_type) == string(unique_cts(3))); % oligodendrocyte
t_idxs = find(cellstr(cell_type) == string(unique_cts(4))); % t-cell


% plot
for i=1:length(marker_genes)
    marker_idx = marker_idxs(i);
    
    mp_counts = sum(expr_norm(marker_idx, mp_idxs) > 0) / size(expr_norm, 2);
    tumor_counts = sum(expr_norm(marker_idx, tumor_idxs) > 0) / size(expr_norm, 2);
    ol_counts = sum(expr_norm(marker_idx, ol_idxs) > 0) / size(expr_norm, 2);
    t_counts = sum(expr_norm(marker_idx, t_idxs) > 0) / size(expr_norm, 2);
    
    cell_counts = [mp_counts, tumor_counts, ol_counts, t_counts] .* 100;
    
    fig = figure;
    bar(categorical(unique_cts), cell_counts);
    ylabel("% nonzero cells");
    title(marker_genes(i));
    
end


%% Q3 DEGs

% 3.1
% subset of malignant cell expressions
% adult cells
pedia_idxs = find(string(meta.GBM_type) == 'Pediatric' & string(cell_type) == 'Malignant');
adult_idxs = find(string(meta.GBM_type) == 'Adult' & string(cell_type) == 'Malignant');

expr_pedia = expr_norm(:, pedia_idxs);
expr_adult = expr_norm(:, adult_idxs);
expr_maglinant = expr_norm(:, cat(1, pedia_idxs, adult_idxs));


% identify DEGs

[T_pedia,Tup_pedia,Tdn_pedia] = sc_deg(expr_pedia, expr_adult, gene_list);
[logfc_pedia, logfc_pedia_idxs] = sort(Tup_pedia.avg_log2FC, 'descend');

[T_adult,Tup_adult,Tdn_adult] = sc_deg(expr_adult, expr_pedia, gene_list);
[logfc_adult, logfc_adult_idxs] = sort(Tup_adult.avg_log2FC, 'descend');

% discard Infs
idx = find(logfc_pedia ~=Inf & logfc_pedia ~= -Inf);
top_degs_pedia = Tup_pedia.gene(logfc_pedia_idxs(idx(1):idx(1)+4));

idx = find(logfc_adult ~=Inf & logfc_adult ~= -Inf);
top_degs_adult = Tup_adult.gene(logfc_adult_idxs(idx(1):idx(1)+4));


% 3.2
tumor_subtypes = meta.GBM_type(malignant_idxs);

% First plot UMAP colored with GBM assignments
gscatter(s_malignant(:, 1), s_malignant(:, 2), tumor_subtypes, '', '.', 4);
title("GBM assignments: pediatric / adult");

% Color with top 2 pediatric genes
pedia_g1_idx = find(string(gene_list) == top_degs_pedia(1));
pedia_g2_idx = find(string(gene_list) == top_degs_pedia(2));
adult_g1_idx = find(string(gene_list) == top_degs_adult(1));
adult_g2_idx = find(string(gene_list) == top_degs_adult(2));

fig = figure;
subplot(1,4,1)
gscatter(s_malignant(:, 1), s_malignant(:, 2), expr_norm(pedia_g1_idx, malignant_idxs), viridis_white, '.', 4, 'off');
title(top_degs_pedia(1));

subplot(1,4,2)
gscatter(s_malignant(:, 1), s_malignant(:, 2), expr_norm(pedia_g2_idx, malignant_idxs), viridis_white, '.', 4, 'off');
title(top_degs_pedia(2));

subplot(1,4,3)
gscatter(s_malignant(:, 1), s_malignant(:, 2), expr_norm(adult_g1_idx, malignant_idxs), viridis_white, '.', 4, 'off');
title(top_degs_adult(1));

subplot(1,4,4)
gscatter(s_malignant(:, 1), s_malignant(:, 2), expr_norm(adult_g2_idx, malignant_idxs), viridis_white, '.', 4, 'off');
title(top_degs_adult(2));
