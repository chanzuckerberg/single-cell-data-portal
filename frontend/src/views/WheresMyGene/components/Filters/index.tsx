import { ComplexFilter, DefaultMenuSelectOption } from "czifui";
import { memo, useMemo } from "react";
import { Filters as IFilters } from "../../common/types";

interface Props {
  filters: IFilters;
  onFiltersChange: (
    key: keyof IFilters,
    options: DefaultMenuSelectOption[] | null
  ) => void;
}

const MOCK_DATASETS = [{ name: "dataset 1" }, { name: "dataset 2" }];
const MOCK_DEVELOPMENT_STAGES = [
  { name: "development stage 1" },
  { name: "development stage 2" },
];
const MOCK_DISEASES = [{ name: "disease 1" }, { name: "disease 2" }];
const MOCK_ETHNICITIES = [{ name: "ethnicity 1" }, { name: "ethnicity 2" }];
const MOCK_SEXES = [{ name: "sex 1" }, { name: "sex 2" }];

export default memo(function Filters({
  filters,
  onFiltersChange,
}: Props): JSX.Element {
  const handleDatasetsChange = useMemo(
    () => handleFilterChange("datasets", onFiltersChange),
    [onFiltersChange]
  );

  const handleDevelopmentStagesChange = useMemo(
    () => handleFilterChange("developmentStages", onFiltersChange),
    [onFiltersChange]
  );

  const handleDiseasesChange = useMemo(
    () => handleFilterChange("diseases", onFiltersChange),
    [onFiltersChange]
  );

  const handleEthnicitiesChange = useMemo(
    () => handleFilterChange("ethnicities", onFiltersChange),
    [onFiltersChange]
  );

  const handleSexesChange = useMemo(
    () => handleFilterChange("sexes", onFiltersChange),
    [onFiltersChange]
  );

  return (
    <div>
      <ComplexFilter
        multiple
        label="Dataset"
        options={MOCK_DATASETS}
        onChange={handleDatasetsChange}
        value={filters.datasets}
      />
      <ComplexFilter
        multiple
        label="Development Stage"
        options={MOCK_DEVELOPMENT_STAGES}
        onChange={handleDevelopmentStagesChange}
        value={filters.developmentStages}
      />
      <ComplexFilter
        multiple
        label="Disease"
        options={MOCK_DISEASES}
        onChange={handleDiseasesChange}
        value={filters.diseases}
      />
      <ComplexFilter
        multiple
        label="Ethnicity"
        options={MOCK_ETHNICITIES}
        onChange={handleEthnicitiesChange}
        value={filters.ethnicities}
      />
      <ComplexFilter
        multiple
        label="Sex"
        options={MOCK_SEXES}
        onChange={handleSexesChange}
        value={filters.sexes}
      />
    </div>
  );
});

function handleFilterChange(
  key: keyof IFilters,
  onFiltersChange: Props["onFiltersChange"]
): (options: DefaultMenuSelectOption[] | null) => void {
  return (options: DefaultMenuSelectOption[] | null): void => {
    onFiltersChange(key, options);
  };
}
