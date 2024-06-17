import {
  Wrapper,
  StepHeader,
  WordPop,
  StepSubHeader,
  CellGroupTitle,
  RunButton,
  RunButtonWrapper,
  FlexRow,
  QuerySelectorWrapper,
  TwoPanelLayout,
  CellCountTitle,
  Spinner,
} from "./style";
import QueryGroupFilters from "./components/Filters";
import Organism from "./components/Organism";
import DeResults from "./components/DeResults";
import Loader from "./components/Loader";
import Method from "./components/Method";
import OverlapBehavior from "./components/OverlapBehavior";
import {
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER,
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER,
  DIFFERENTIAL_EXPRESSION_FILTERS_LOADING_SPINNER,
  DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT,
  DIFFERENTIAL_EXPRESSION_FIND_GENES_BUTTON,
  DIFFERENTIAL_EXPRESSION_METHOD_INFO_TEXT,
} from "../../common/constants";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { useConnect } from "./connect";

export default function DifferentialExpression(): JSX.Element {
  const {
    isLoading,
    setIsLoading,
    queryGroup1,
    queryGroup2,
    canRunDifferentialExpression,
    handleRunDifferentialExpression,
    nCellsGroup1,
    isLoadingGroup1,
    nCellsGroup2,
    isLoadingGroup2,
  } = useConnect();

  return (
    <TwoPanelLayout>
      <div className="leftPanel">
        <Wrapper>
          {isLoading && <Loader />}

          <QuerySelectorWrapper>
            <StepHeader>
              <WordPop>Differential</WordPop> Expression
            </StepHeader>
            <StepSubHeader
              data-testid={DIFFERENTIAL_EXPRESSION_METHOD_INFO_TEXT}
            >
              Find differentially expressed genes between custom group of cells
              across the CELLxGENE data corpus. For additional help and
              information, read our documentation.
              <br />
              <br />
              This tool uses Welch&apos;s t-test to identify differentially
              expressed genes between groups of cell in the CELLxGENE Census.
              While the t-test performs reasonably well on individual datasets{" "}
              <a
                href="https://www.mdpi.com/2073-4425/12/12/1947"
                target="_blank"
                rel="noopener noreferrer"
              >
                [1]
              </a>
              <a
                href="https://www.nature.com/articles/s41467-019-12266-7"
                target="_blank"
                rel="noopener noreferrer"
              >
                [2]
              </a>
              <a
                href="https://www.nature.com/articles/nmeth.4612"
                target="_blank"
                rel="noopener noreferrer"
              >
                [3]
              </a>
              , its performance on concatenated, non-integrated datasets has not
              been extensively evaluated. We recommend using this tool for
              preliminary investigations and following up with a more robust
              method for formal analysis. Learn more about our data filtering
              and normalization{" "}
              <a
                href="/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_3__Gene%20Expression%20Data%20Processing"
                target="_blank"
                rel="noopener noreferrer"
                onClick={() => track(EVENTS.DE_DOCUMENTATION_CLICKED)}
              >
                here
              </a>
              .
            </StepSubHeader>
            <FlexRow>
              <Organism />
              <Method />
              <OverlapBehavior />
            </FlexRow>
            <FlexRow>
              <div data-testid={DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER}>
                <CellGroupTitle>
                  Cell Group 1
                  <CellCountTitle
                    data-testid={DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT}
                  >
                    {isLoadingGroup1 && (
                      <Spinner
                        data-testid={
                          DIFFERENTIAL_EXPRESSION_FILTERS_LOADING_SPINNER
                        }
                      />
                    )}
                    {nCellsGroup1.toLocaleString()} cells
                  </CellCountTitle>
                </CellGroupTitle>
                <QueryGroupFilters
                  key={`query-group-1`}
                  queryGroup={queryGroup1}
                  isQueryGroup1={true}
                />
              </div>

              <div data-testid={DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER}>
                <CellGroupTitle>
                  Cell Group 2
                  <CellCountTitle
                    data-testid={DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT}
                  >
                    {isLoadingGroup2 && (
                      <Spinner
                        data-testid={
                          DIFFERENTIAL_EXPRESSION_FILTERS_LOADING_SPINNER
                        }
                      />
                    )}
                    {nCellsGroup2.toLocaleString()} cells
                  </CellCountTitle>
                </CellGroupTitle>
                <QueryGroupFilters
                  key={`query-group-2`}
                  queryGroup={queryGroup2}
                  isQueryGroup1={false}
                />
              </div>
            </FlexRow>
            <RunButtonWrapper>
              <RunButton
                color="primary"
                size="large"
                variant="contained"
                onClick={handleRunDifferentialExpression}
                disabled={!canRunDifferentialExpression}
                data-testid={DIFFERENTIAL_EXPRESSION_FIND_GENES_BUTTON}
              >
                Find genes
              </RunButton>
            </RunButtonWrapper>
          </QuerySelectorWrapper>
        </Wrapper>
      </div>
      <div className="rightPanel">
        <DeResults setIsLoading={setIsLoading} />
      </div>
    </TwoPanelLayout>
  );
}

// interface CirclesOverlapBothProps {
//   isActive: boolean;
//   onClick: () => void;
// }
// const CirclesOverlapBoth = ({ isActive, onClick }: CirclesOverlapBothProps) => {
//   const fillColor = isActive ? "blue" : "gray"; // Change colors as needed
//   return (
//     <div onClick={onClick}>
//       <Icon sdsIcon="CirclesOverlap2" sdsSize="s" sdsType="button" />
//     </div>
//   );
// };

// interface CirclesOverlapLeftProps {
//   isActive: boolean;
//   onClick: () => void;
// }
// const CirclesOverlapLeft = ({ isActive, onClick }: CirclesOverlapLeftProps) => {
//   const fillColor = isActive ? "blue" : "gray"; // Change colors as needed
//   return (
//     <SvgIcon onClick={onClick}>
//       <svg
//         width="16"
//         height="16"
//         viewBox="0 0 16 16"
//         xmlns="http://www.w3.org/2000/svg"
//       >
//         <path
//           d="M15.25 8C15.25 10.3333 13.3166 12.25 10.9 12.25C8.48354 12.25 6.55005 10.3333 6.55005 8C6.55005 5.66672 8.48354 3.75 10.9 3.75C13.3166 3.75 15.25 5.66672 15.25 8Z"
//           fill="white"
//           stroke={fillColor}
//           stroke-width="1.5"
//         />
//         <path
//           d="M9.45 8C9.45 10.3333 7.51651 12.25 5.1 12.25C2.68349 12.25 0.75 10.3333 0.75 8C0.75 5.66672 2.68349 3.75 5.1 3.75C7.51651 3.75 9.45 5.66672 9.45 8Z"
//           fill={fillColor}
//           stroke={fillColor}
//           stroke-width="1.5"
//         />
//       </svg>
//     </SvgIcon>
//   );
// };

// interface CirclesOverlapRightProps {
//   isActive: boolean;
//   onClick: () => void;
// }
// const CirclesOverlapRight = ({
//   isActive,
//   onClick,
// }: CirclesOverlapRightProps) => {
//   const fillColor = isActive ? "blue" : "gray"; // Change colors as needed
//   return (
//     <SvgIcon onClick={onClick}>
//       <svg
//         width="16"
//         height="16"
//         viewBox="0 0 16 16"
//         xmlns="http://www.w3.org/2000/svg"
//       >
//         <path
//           d="M9.45 8C9.45 10.3333 7.51651 12.25 5.1 12.25C2.68349 12.25 0.75 10.3333 0.75 8C0.75 5.66672 2.68349 3.75 5.1 3.75C7.51651 3.75 9.45 5.66672 9.45 8Z"
//           fill="white"
//           stroke={fillColor}
//           stroke-width="1.5"
//         />
//         <path
//           d="M15.25 8C15.25 10.3333 13.3166 12.25 10.9 12.25C8.48354 12.25 6.55005 10.3333 6.55005 8C6.55005 5.66672 8.48354 3.75 10.9 3.75C13.3166 3.75 15.25 5.66672 15.25 8Z"
//           fill={fillColor}
//           stroke={fillColor}
//           stroke-width="1.5"
//         />
//       </svg>
//     </SvgIcon>
//   );
// };
