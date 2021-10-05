export const INCOMPLETE_GENE_MIGRATIONS_COLLECTIONS = [
  // Actual Collections:
  // "5d445965-6f1a-4b68-ba3a-b8f765155d3a",
  // "a3c37598-b8e8-47db-92a9-eb91a77c9529",
  // "558385a4-b7b7-4eca-af0c-9e54d010e8dc",
  // "ae1420fe-6630-46ed-8b3d-cc6056a66467",
  // "45f0f67d-4b69-4a3c-a4e8-a63b962e843f",
  // "d0e9c47b-4ce7-4f84-b182-eddcfa0b2658",
  // "6e8c5415-302c-492a-a5f9-f29c57ff18fb",
  // "38833785-fac5-48fd-944a-0f62a4c23ed1",
  // "b9fc3d70-5a72-4479-a046-c2cc1ab19efc",
  // "5e469121-c203-4775-962d-dcf2e5d6a472",
  // "2f75d249-1bec-459b-bf2b-b86221097ced",
  // "4f586cb6-972b-4ef7-a4ef-3c3800a3c004",
  // Test Rdev Collections:
  "423cdd2a-d096-4e87-a9a9-18e5168b45e7",
  "50b7e0d4-ecc8-4576-a815-6529e02c883c",
];

export const INCOMPLETE_GENE_MIGRATIONS_HEADER = "Important Note";
export const INCOMPLETE_GENE_MIGRATIONS_MSG = `
  The cellxgene Data Portal was recently migrated to an enhanced schema to 
  optimize reusability. The cellxgene team is still migrating a small subset 
  of genes with complex mappings to ENSEMBL 38, the human gene reference 
  required by the schema. As a result, more genes may be added in the coming 
  weeks. The current data can be used with confidence that the genes that are 
  present are correctly mapped.
`;
