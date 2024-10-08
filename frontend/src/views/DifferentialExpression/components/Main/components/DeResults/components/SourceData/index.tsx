import { Tooltip } from "@czi-sds/components";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  CellGroupHeader,
  Content,
  InfoText,
  StyledDivTableCell,
  StyledDivTableRow,
  StyledLabel,
  StyledLink,
  StyledSpan,
  StyledValue,
  Wrapper,
} from "./style";
import { DivTable, DivTableHead } from "src/views/WheresMyGeneV2/common/styles";
import { ROUTES } from "src/common/constants/routes";
import { useConnect } from "./connect";
import { SOURCE_DATA_DISCLAIMER, SOURCE_DATA_INFO_TEXT } from "./constants";
import {
  StyledTooltip,
  TooltipContent,
  TooltipButton,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";
import { DIFFERENTIAL_EXPRESSION_SOURCE_DATA_SIDEBAR } from "src/views/DifferentialExpression/common/constants";

export default function SourceData(): JSX.Element {
  const { collections1, collections2, setHoverStartTime, handleBadgeHoverEnd } =
    useConnect();

  return (
    <Wrapper data-testid={DIFFERENTIAL_EXPRESSION_SOURCE_DATA_SIDEBAR}>
      <Content>
        <InfoText>
          {SOURCE_DATA_INFO_TEXT}{" "}
          <a
            href={ROUTES.WMG_DOCS_DATA_PROCESSING}
            target="_blank"
            rel="noreferrer noopener"
            data-testid="documentation-link"
          >
            exceptions and processing notes
          </a>
          .
        </InfoText>
        <InfoText>{SOURCE_DATA_DISCLAIMER}</InfoText>
        <DivTable data-testid="source-data-list">
          {[collections1, collections2].map((collections, index) => (
            <>
              <CellGroupHeader>
                {index === 0
                  ? "Cell Group 1 Collections"
                  : "Cell Group 2 Collections"}
              </CellGroupHeader>
              <DivTableHead>
                <StyledDivTableCell>Collection</StyledDivTableCell>
                <StyledDivTableCell align>Datasets</StyledDivTableCell>
              </DivTableHead>
              {Object.values(collections).map(({ name, url, datasets }) => (
                <StyledDivTableRow key={name}>
                  <StyledDivTableCell>
                    {name ? (
                      <StyledLink
                        href={url}
                        target="_blank"
                        rel="noopener"
                        onClick={() => {
                          track(EVENTS.VIEW_COLLECTION_PAGE_CLICKED, {
                            collection_name: name,
                          });
                        }}
                      >
                        {name}
                      </StyledLink>
                    ) : (
                      <StyledSpan>
                        Failed to load collection information
                      </StyledSpan>
                    )}
                  </StyledDivTableCell>
                  <StyledDivTableCell align>
                    <Tooltip
                      sdsStyle="light"
                      placement="left-end"
                      width="wide"
                      className="fmg-tooltip-icon"
                      arrow
                      onOpen={() => setHoverStartTime(Date.now())}
                      onClose={handleBadgeHoverEnd}
                      title={
                        <StyledTooltip>
                          <TooltipContent>
                            {datasets.map((dataset) => (
                              <StyledValue key={dataset.id}>
                                {dataset.label}
                              </StyledValue>
                            ))}
                          </TooltipContent>
                        </StyledTooltip>
                      }
                    >
                      <TooltipButton
                        sdsStyle="minimal"
                        sdsType="secondary"
                        isAllCaps={false}
                      >
                        <StyledLabel>{datasets.length}</StyledLabel>
                      </TooltipButton>
                    </Tooltip>
                  </StyledDivTableCell>
                </StyledDivTableRow>
              ))}
            </>
          ))}
        </DivTable>
      </Content>
    </Wrapper>
  );
}
