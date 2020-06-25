export interface Contributor {
  institution: string;
  name: string;
}

export interface Project {
  assays: string[];
  biosample_categories: string[];
  cell_count: number;
  contributors: Contributor[];
  cxg_enabled: boolean;
  description: string;
  diseases: string[];
  id: string;
  label: string;
  organs: string[];
  paired_end: string[];
  publication_title: string;
  species: string[];
  title: string;
}

// id: "HCA-AnalysisFile-001a3fdd-c9e6-4871-9e48-a840d52ecdf9"
// filename: "b29eeb85-53ff-4786-9101-3e241e6dc250_rsem.bam"
// file_format: "bam"
// file_size: 0
// species: "homo sapiens"
// library_construction_method_ontology: ""
// tissue_ontology: ""
export interface File {
  file_format: string;
  file_size: number;
  filename: string;
  id: string;
  library_construction_method_ontology: string;
  species: string;
  tissue_ontology: string;
}
