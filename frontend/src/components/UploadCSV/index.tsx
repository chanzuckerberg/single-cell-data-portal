/**
 * (thuang): TEMP. Remove when we do https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/corpora-data-portal/917
 * This component is used temporarily for e2e testing and trying out different
 * CSV files
 */
import { ChangeEvent, FC, useState } from "react";
import { GeneSet } from "src/common/entities";
import parseGeneSetsCSV from "src/common/utils/parseGeneSetsCSV";

const MOCK_EXISTING_GENE_SETS: GeneSet[] = [
  {
    description: "foo description",
    id: "foo id",
    linked_datasets: [],
    name: "foo",
  },
];

const UploadCSV: FC = () => {
  const [csvResult, setCsvResult] = useState<unknown>(null);

  const handleCSV = async (event: ChangeEvent) => {
    const target = event.target as HTMLInputElement;
    const file = target?.files && target?.files[0];

    if (!file) return;

    try {
      const result = await parseGeneSetsCSV(file, MOCK_EXISTING_GENE_SETS);

      console.log("parseGeneSetsCSV result", result);

      setCsvResult(result);
    } catch (error) {
      console.log("------something went wrong", error);
    }
  };

  return (
    <>
      <input
        data-testid="upload-csv"
        type="file"
        accept=".csv"
        onChange={handleCSV}
      />
      {csvResult && (
        <div data-testid="csv-result">{JSON.stringify(csvResult)}</div>
      )}
    </>
  );
};

export default UploadCSV;
