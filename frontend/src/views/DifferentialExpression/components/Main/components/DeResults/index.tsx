import { useState, useEffect, useContext, useMemo } from "react";
import { useDifferentialExpression } from "src/common/queries/differentialExpression";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import {
  CellGroupTitle,
  CellGroupTitleWrapper,
  CellGroupWrapper,
  EffectSizeHeaderWrapper,
  EffectSizeIndicator,
  FilterTagsWrapper,
  InstructionsBody,
  InstructionsHeader,
  InstructionsWrapper,
  ResultsHeader,
  ResultsWrapper,
  StyledTag,
  StyledTextField,
  TableHeaderWrapper,
  TableWrapper,
} from "./style";
import { QueryGroup } from "src/views/DifferentialExpression/common/store/reducer";
import { clearSubmittedQueryGroups } from "src/views/DifferentialExpression/common/store/actions";
import { Pagination } from "@mui/material";
import Table from "src/views/CellGuide/components/CellGuideCard/components/common/Table";
import { ButtonIcon, Tooltip } from "@czi-sds/components";

interface DifferentialExpressionRow {
  name: string;
  pValue: string;
  effectSize: string;
}

interface Props {
  setIsLoading: (isLoading: boolean) => void;
}
export default function DeResults({ setIsLoading }: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { data, isLoading } = useDifferentialExpression();
  const { differentialExpressionResults: rawDifferentialExpressionResults } =
    data;

  const [differentialExpressionResults, setDifferentialExpressionResults] =
    useState<DifferentialExpressionRow[]>([]);

  const {
    organismId,
    submittedQueryGroups,
    queryGroups,
    submittedQueryGroupsWithNames: queryGroupsWithNames,
  } = useContext(StateContext);

  useEffect(() => {
    if (!organismId || isLoading) return;

    // map ids to name
    const formattedDeResults = rawDifferentialExpressionResults.map(
      (diffExpResult) => {
        return {
          name: diffExpResult.gene_symbol,
          pValue: diffExpResult.p_value.toExponential(3),
          effectSize: diffExpResult.effect_size.toFixed(3),
        };
      }
    );

    setDifferentialExpressionResults(formattedDeResults);
  }, [rawDifferentialExpressionResults, isLoading, organismId]);

  const namesToShow: string[][] = [];
  const { queryGroup1, queryGroup2 } = queryGroupsWithNames ?? {};
  for (const [index, queryGroupWithNames] of [
    queryGroup1,
    queryGroup2,
  ].entries()) {
    namesToShow.push([]);
    for (const key in queryGroupWithNames) {
      for (const value of queryGroupWithNames[key as keyof QueryGroup]) {
        namesToShow[index].push(value);
      }
    }
  }

  useEffect(() => {
    if (!submittedQueryGroups || !dispatch) return;
    if (JSON.stringify(queryGroups) !== JSON.stringify(submittedQueryGroups))
      dispatch(clearSubmittedQueryGroups());
  }, [dispatch, queryGroups, submittedQueryGroups]);

  const showEmpty = !submittedQueryGroups || isLoading;

  useEffect(() => {
    setIsLoading(isLoading);
  }, [isLoading, setIsLoading]);

  return (
    <div>
      {!showEmpty ? (
        <DifferentialExpressionResults
          results={differentialExpressionResults}
        />
      ) : (
        <InstructionsWrapper>
          <InstructionsHeader>Instructions</InstructionsHeader>
          <InstructionsBody>
            <ol>
              <li>
                Select a cell group of interest within the Cell Group 1 box by
                using the dropdown selectors.
                <br />
                <br />
                To copy the same selection over to Cell Group 2, click the copy
                button to the right of each dropdown in Cell Group 1.
              </li>
              <li>
                Within Cell Group 2, select a group that the cell group of
                interest will be compared to.
                <br />
                <br />
                To easily select the inverse of Cell Group 1, click the
                overlapping circle icon to the left of each dropdown in Cell
                Group 2.
              </li>
            </ol>
          </InstructionsBody>
        </InstructionsWrapper>
      )}
    </div>
  );
}

const ROWS_PER_PAGE = 15;
const QUERY_GROUP_KEY_TO_TAG_SUFFIX_MAP = {
  diseases: { single: "disease", plural: "diseases" },
  ethnicities: { single: "ethnicity", plural: "ethnicities" },
  sexes: { single: "sex", plural: "sexes" },
  tissues: { single: "tissue", plural: "tissues" },
  cellTypes: { single: "cell type", plural: "cell types" },
  publicationCitations: { single: "publication", plural: "publications" },
  developmentStages: {
    single: "development stage",
    plural: "development stages",
  },
};

interface DifferentialExpressionResultsProps {
  results: DifferentialExpressionRow[];
}
const DifferentialExpressionResults = ({
  results,
}: DifferentialExpressionResultsProps) => {
  const [page, setPage] = useState(1);
  const [searchQuery, setSearchQuery] = useState("");
  const [pvalueFilter, setPvalueFilter] = useState("");
  const [effectSizeFilter, setEffectSizeFilter] = useState("");
  const [sortDirection, setSortDirection] = useState<"asc" | "desc">("desc");

  const handlePageChange = (
    _event: React.ChangeEvent<unknown>,
    page: number
  ) => {
    setPage(page);
  };

  const sortedAndFilteredResults = useMemo(() => {
    return results
      .filter((result) =>
        result.name.toLowerCase().includes(searchQuery.toLowerCase())
      )
      .filter((result) => applyFilter(parseFloat(result.pValue), pvalueFilter))
      .filter((result) =>
        applyFilter(parseFloat(result.effectSize), effectSizeFilter)
      )
      .sort((a, b) => {
        if (sortDirection === "asc") {
          return parseFloat(a.effectSize) - parseFloat(b.effectSize);
        } else if (sortDirection === "desc") {
          return parseFloat(b.effectSize) - parseFloat(a.effectSize);
        }
        return 0;
      });
  }, [results, searchQuery, pvalueFilter, effectSizeFilter, sortDirection]);

  const pageCount = Math.ceil(sortedAndFilteredResults.length / ROWS_PER_PAGE);

  const columnIdToName: Record<
    keyof DifferentialExpressionRow,
    string | JSX.Element
  > = useMemo(() => {
    const handleSearch = (event: React.ChangeEvent<HTMLInputElement>) => {
      setSearchQuery(event.target.value);
      setPage(1);
    };
    const handlePvalueFilter = (event: React.ChangeEvent<HTMLInputElement>) => {
      setPvalueFilter(event.target.value);
      setPage(1);
    };
    const handleEffectSizeFilter = (
      event: React.ChangeEvent<HTMLInputElement>
    ) => {
      setEffectSizeFilter(event.target.value);
      setPage(1);
    };

    const handleSortDirectionChange = () => {
      setSortDirection((prevDirection) =>
        prevDirection === "asc" ? "desc" : "asc"
      );
    };
    return {
      name: (
        <TableHeaderWrapper>
          Gene
          <StyledTextField
            variant="outlined"
            onChange={handleSearch}
            placeholder="e.g. JCHAIN"
          />
        </TableHeaderWrapper>
      ),
      pValue: (
        <TableHeaderWrapper>
          P-value
          <StyledTextField
            placeholder="e.g <0.05"
            variant="outlined"
            onChange={handlePvalueFilter}
          />
        </TableHeaderWrapper>
      ),
      effectSize: (
        <TableHeaderWrapper>
          <EffectSizeHeaderWrapper>
            Effect Size{" "}
            <ButtonIcon
              onClick={handleSortDirectionChange}
              sdsIcon={sortDirection === "asc" ? "chevronUp" : "chevronDown"}
              sdsSize="small"
            />
          </EffectSizeHeaderWrapper>
          <StyledTextField
            placeholder="e.g >1.0"
            variant="outlined"
            onChange={handleEffectSizeFilter}
          />
        </TableHeaderWrapper>
      ),
    };
  }, [sortDirection]);
  return (
    <ResultsWrapper>
      <ResultsHeader>Results</ResultsHeader>
      <CellGroupWrapper>
        <CellGroupTitleWrapper>
          <CellGroupTitle>Cell Group 1</CellGroupTitle>
          <EffectSizeIndicator>{"(+) Effect Size"}</EffectSizeIndicator>
        </CellGroupTitleWrapper>
        <FilterTagsWrapper>
          <QueryGroupTags isQueryGroup1 />
        </FilterTagsWrapper>
      </CellGroupWrapper>
      <CellGroupWrapper>
        <CellGroupTitleWrapper>
          <CellGroupTitle>Cell Group 2</CellGroupTitle>
          <EffectSizeIndicator>{"(-) Effect Size"}</EffectSizeIndicator>
        </CellGroupTitleWrapper>
        <FilterTagsWrapper>
          <QueryGroupTags />
        </FilterTagsWrapper>
      </CellGroupWrapper>
      <TableWrapper>
        <Table<DifferentialExpressionRow>
          columns={["name", "pValue", "effectSize"]}
          rows={sortedAndFilteredResults.slice(
            (page - 1) * ROWS_PER_PAGE,
            page * ROWS_PER_PAGE
          )}
          columnIdToName={columnIdToName}
        />

        <Pagination count={pageCount} page={page} onChange={handlePageChange} />
      </TableWrapper>
    </ResultsWrapper>
  );
};

const QueryGroupTags = ({ isQueryGroup1 }: { isQueryGroup1?: boolean }) => {
  const { queryGroupsWithNames } = useContext(StateContext);
  const queryGroup = isQueryGroup1
    ? queryGroupsWithNames.queryGroup1
    : queryGroupsWithNames.queryGroup2;
  const nonEmptyQueryGroupKeys = useMemo(() => {
    return Object.keys(queryGroup).filter(
      (key) => queryGroup[key as keyof QueryGroup].length > 0
    );
  }, [queryGroup]);

  return (
    <>
      {nonEmptyQueryGroupKeys.map((key) => {
        const queryGroupKey = key as keyof QueryGroup;
        const selected = queryGroup[queryGroupKey];
        const suffix = QUERY_GROUP_KEY_TO_TAG_SUFFIX_MAP[queryGroupKey];
        const label = `${selected.length} ${
          selected.length > 1 ? suffix.plural : suffix.single
        }`;
        return (
          <Tooltip
            key={`${key}-tooltip`}
            sdsStyle="light"
            placement="top"
            width="wide"
            leaveDelay={0}
            title={selected.map((value, index) => (
              <div key={`value-${value}-${index}`}>{value}</div>
            ))}
          >
            <span>
              <StyledTag
                key={key}
                color="gray"
                sdsType="secondary"
                label={label}
              />
            </span>
          </Tooltip>
        );
      })}
    </>
  );
};

const parseExpressions = (expression: string) => {
  const expressions = expression
    .split(",")
    .map((expr) => {
      const match = expr.match(/([<>]=?)(\d*\.?\d+(?:e[+-]?\d+)?)/);
      if (!match) return null;

      const [, operator, valueStr] = match;
      const value = parseFloat(valueStr);
      return { operator, value };
    })
    .filter((expr) => expr !== null);

  return expressions.length > 0 ? expressions : null;
};

const applyFilter = (resultValue: number, filterExpression: string) => {
  const parsedExpressions = parseExpressions(filterExpression);
  if (!parsedExpressions) return true;
  const filteredParsedExpressions = parsedExpressions.filter(
    (expr) => expr !== null
  ) as {
    operator: string;
    value: number;
  }[];
  return filteredParsedExpressions.every(({ operator, value }) => {
    switch (operator) {
      case "<":
        return resultValue < value;
      case "<=":
        return resultValue <= value;
      case ">":
        return resultValue > value;
      case ">=":
        return resultValue >= value;
      default:
        return true;
    }
  });
};
