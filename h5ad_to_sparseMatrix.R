library(rhdf5)
library(Matrix)

root <- rhdf5::H5Fopen('/home/jupyter/notebooks/data/pdac_19-587/nucSeq/inferCNV_files/inputs/all/infercnv_anndata.h5ad')

target <- "/X"


i_path <- paste0(target, "/indices")
p_path <- paste0(target, "/indptr")
x_path <- paste0(target, "/data")

print('reading indices')
i <- as.vector(unlist(rhdf5::h5read(root, i_path)))
print('reading pointers')
p <- as.vector(unlist(rhdf5::h5read(root, p_path)))
print('reading values')
x <- as.vector(unlist(rhdf5::h5read(root, x_path)))
print('reading observations')
o <- as.vector(rhdf5::h5read(root, "/obs/_index"))
print('reading variables')
v <- as.vector(rhdf5::h5read(root, "/var/_index"))

print("Reading dimensions")
dims <- c(length(v), length(o))

print("Assembling dgCMatrix")
m <- Matrix::sparseMatrix(i = i,
                            p = p,
                            x = x,
                            index1 = FALSE,
                            dims = dims)

rownames(m) <- v
colnames(m) <- o

saveRDS(m, '/home/jupyter/notebooks/data/pdac_19-587/nucSeq/inferCNV_files/inputs/all/sparse_counts.Rds')
