library("Seurat")

# If running this in Docker, mount the folder with your RDS to /data
# docker run -v /path/to/folder:/data docker_image
f <- readRDS(file = "/data/dataset2.rds")

if (f@misc$schema_version != "3.0.0") {
    print("ERROR: wrong schema version")
} else {
    print("Correct schema version")
}


if("donor_id" %in% colnames(f@meta.data)) {
    print("donor_id is in")
} else {
    print("ERROR: donor_id missing")
}

if("suspension_type" %in% colnames(f@meta.data)) {
    print("suspension_type is in")
} else {
    print("ERROR: suspension_type missing")
}

if(!"self_reported_ethnicity" %in% colnames(f@meta.data)) {
    print("ERROR: self_reported_ethnicity is missing")
} 

if("ethnicity" %in% colnames(f@meta.data)) {
    print("ERROR: ethnicity is in. this is a deprecated field")
} 