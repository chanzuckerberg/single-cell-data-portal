import { H5 } from "@blueprintjs/core";
import { CollectionMigrationCallout as MigrationCallout } from "src/components/Collection/components/CollectionMigrationCallout/style";

const INCOMPLETE_GENE_MIGRATIONS_COLLECTIONS = [
  "2f75d249-1bec-459b-bf2b-b86221097ced",
];

interface Props {
  collectionId: string;
}

export default function CollectionMigrationCallout({
  collectionId,
}: Props): JSX.Element | null {
  const isIncompleteMigrationCollection =
    INCOMPLETE_GENE_MIGRATIONS_COLLECTIONS.includes(collectionId);

  return isIncompleteMigrationCollection ? (
    <MigrationCallout>
      <H5>Important Note</H5>
      <p>
        CELLxGENE Discover was migrated to an enhanced schema to optimize
        reusability. Some genes from this dataset were eliminated in this
        process. The data can be used with confidence that the genes that are
        present are correctly mapped.
      </p>
    </MigrationCallout>
  ) : null;
}
