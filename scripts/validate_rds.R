library("Seurat")

# If running this in Docker, mount the folder with your RDS to /data
# docker run -v /path/to/folder:/data docker_image
f <- readRDS(file = "/data/dataset.rds")

if (f@misc$schema_version != "3.0.0") {
    print("Wrong schema version")
} else {
    print("Correct schema version")
}


if("donor_id" %in% colnames(f@meta.data)) {
    print("donor_id is in")
} else {
    print("donor_id missing")
}