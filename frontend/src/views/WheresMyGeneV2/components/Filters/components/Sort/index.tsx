import { Tooltip } from "@czi-sds/components";
import { Wrapper, FilterLabel, StyledDropdown } from "../common/style";
import { ViewOptionsWrapper } from "./style";
import {
  StyledIconImage,
  StyledTooltip,
  TooltipButton,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import questionMarkIcon from "src/common/images/question-mark-icon.svg";
import { ROUTES } from "src/common/constants/routes";
import {
  SELECT_TISSUE_GENE_TEXT,
  SORT_CELL_TYPES_TOOLTIP_TEXT,
  SORT_GENES_TOOLTIP_TEXT,
} from "src/views/WheresMyGene/common/constants";
import { Props } from "./types";
import { useConnect } from "./connect";
import { CELL_TYPE_OPTIONS, GENE_OPTIONS } from "./constants";

export default function Sort({ areFiltersDisabled }: Props): JSX.Element {
  const {
    InputDropdownProps,
    cellTypeSelectedOption,
    geneSelectedOption,
    cellTypesOnChange,
    genesOnChange,
  } = useConnect({ areFiltersDisabled });

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
}
