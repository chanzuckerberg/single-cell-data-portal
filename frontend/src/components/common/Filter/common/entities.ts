import { Collection, IS_PRIMARY_DATA, Ontology } from "src/common/entities";

/* Join of collection and aggregated dataset information optimized for filtering */
export interface FilterableCollection {
  assayAggregated: Ontology[];
  cell_typeAggregated: Ontology[];
  diseaseAggregated: Ontology[];
  is_primary_dataAggregated: IS_PRIMARY_DATA[];
  organismAggregated: Ontology[];
  published_at: number;
  revised_at?: number;
  sexAggregated: Ontology[];
  tissueAggregated: Ontology[];
}

/* Join of dataset and collection information optimized for filtering */
export interface FilterableDataset {
  assay: Ontology[];
  cell_count: number | null;
  cell_type: Ontology[];
  collection_id: Collection["id"];
  collection_name: Collection["name"];
  disease: Ontology[];
  filterableCollection?: FilterableCollection; // Only specified when filtering over collections TODO(cc) revisit
  id: string;
  is_primary_data: IS_PRIMARY_DATA;
  name: string;
  organism: Ontology[];
  published_at: number;
  revised_at?: number;
  sex: Ontology[];
  tissue: Ontology[];
}
