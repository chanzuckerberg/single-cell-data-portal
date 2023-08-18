import {
  InputDropdownProps as IInputDropdownProps,
  Tooltip,
} from "@czi-sds/components";
import { useContext, useMemo } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectSortBy } from "src/views/WheresMyGene/common/store/actions";
import { SORT_BY } from "src/views/WheresMyGene/common/types";
import { Wrapper, FilterLabel, StyledDropdown } from "../common/style";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { ViewOptionsWrapper } from "./style";
import {
  StyledIconImage,
  StyledTooltip,
  TooltipButton,
} from "../../../CellInfoSideBar/style";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { ROUTES } from "src/common/constants/routes";
import {
  SELECT_TISSUE_GENE_TEXT,
  SORT_CELL_TYPES_TOOLTIP_TEXT,
  SORT_GENES_TOOLTIP_TEXT,
} from "src/views/WheresMyGene/common/constants";

const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

const CELL_TYPE_OPTIONS = [
  { id: SORT_BY.CELL_ONTOLOGY, name: "Cell Ontology" },
  { id: SORT_BY.H_CLUSTER, name: "Hierarchical" },
];

const GENE_OPTIONS = [
  { id: SORT_BY.USER_ENTERED, name: "As Entered" },
  { id: SORT_BY.H_CLUSTER, name: "Hierarchical" },
];

interface Props {
  areFiltersDisabled: boolean;
}

export default function Sort({ areFiltersDisabled }: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { sortBy } = useContext(StateContext);

  const InputDropdownProps = useMemo(
    () => ({
      ...DEFAULT_INPUT_DROPDOWN_PROPS,
      disabled: areFiltersDisabled,
    }),
    [areFiltersDisabled]
  );

  const cellTypeSelectedOption = useMemo(() => {
    return (
      CELL_TYPE_OPTIONS.find((option) => option.id === sortBy.cellTypes) ||
      CELL_TYPE_OPTIONS[0]
    );
  }, [sortBy]);

  const geneSelectedOption = useMemo(() => {
    return (
      GENE_OPTIONS.find((option) => option.id === sortBy.genes) ||
      GENE_OPTIONS[0]
    );
  }, [sortBy]);

  return (
    <ViewOptionsWrapper>
      <Wrapper>
        <FilterLabel>
          Sort Cell Types
          <Tooltip
            id="sort-cell-types-tooltip-icon"
            className="sort-cell-types-tooltip-icon"
            sdsStyle="dark"
            placement="right"
            width="default"
            arrow
            title={
              <StyledTooltip>
                {areFiltersDisabled && (
                  <p>
                    <em>{SELECT_TISSUE_GENE_TEXT}</em>
                  </p>
                )}
                <p>{SORT_CELL_TYPES_TOOLTIP_TEXT}</p>
                <a
                  href={ROUTES.WMG_DOCS_ORDERING}
                  rel="noopener"
                  target="_blank"
                >
                  Click to read more.
                </a>
              </StyledTooltip>
            }
          >
            <TooltipButton
              data-testid="sort-cell-types-tooltip-icon"
              sdsStyle="minimal"
              sdsType="secondary"
              isAllCaps={false}
            >
              <StyledIconImage src={questionMarkIcon} />
            </TooltipButton>
          </Tooltip>
        </FilterLabel>

        <StyledDropdown
          data-testid="cell-type-sort-dropdown"
          onChange={cellTypesOnChange}
          label={cellTypeSelectedOption.name}
          options={CELL_TYPE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
          value={cellTypeSelectedOption}
        />
      </Wrapper>
      <Wrapper>
        <FilterLabel>
          Sort Genes
          <Tooltip
            className="sort-genes-tooltip-icon"
            sdsStyle="dark"
            placement="right"
            width="default"
            arrow
            title={
              <StyledTooltip>
                {areFiltersDisabled && (
                  <p>
                    <em>{SELECT_TISSUE_GENE_TEXT}</em>
                  </p>
                )}
                <p>{SORT_GENES_TOOLTIP_TEXT}</p>
                <a
                  href={ROUTES.WMG_DOCS_ORDERING}
                  rel="noopener"
                  target="_blank"
                >
                  Click to read more.
                </a>
              </StyledTooltip>
            }
          >
            <TooltipButton
              data-testid="sort-genes-tooltip-icon"
              sdsStyle="minimal"
              sdsType="secondary"
              isAllCaps={false}
            >
              <StyledIconImage src={questionMarkIcon} />
            </TooltipButton>
          </Tooltip>
        </FilterLabel>
        <StyledDropdown
          data-testid="gene-sort-dropdown"
          onChange={genesOnChange}
          label={geneSelectedOption.name}
          options={GENE_OPTIONS}
          InputDropdownProps={InputDropdownProps}
          value={geneSelectedOption}
        />
      </Wrapper>
    </ViewOptionsWrapper>
  );

  function cellTypesOnChange(
    value: { id?: SORT_BY; name: string } | null
  ): void {
    if (!dispatch || !value || cellTypeSelectedOption.name === value.name)
      return;

    track(EVENTS.WMG_OPTION_SELECT_CELL_TYPES, {
      sort_cell_types_view_option: value.name,
    });

    dispatch(selectSortBy({ cellTypes: value.id as SORT_BY }));
  }

  function genesOnChange(value: { id?: SORT_BY; name: string } | null): void {
    if (!dispatch || !value || geneSelectedOption.name === value.name) return;

    track(EVENTS.WMG_OPTION_SELECT_SORT_GENES, {
      sort_genes_view_option: value.name,
    });

    dispatch(selectSortBy({ genes: value.id as SORT_BY }));
  }
}
