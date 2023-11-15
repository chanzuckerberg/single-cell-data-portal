import { Tooltip } from "@czi-sds/components";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  Content,
  InfoText,
  StyledDivTableCell,
  StyledDivTableRow,
  StyledLabel,
  StyledLink,
  StyledValue,
  Wrapper,
} from "./style";
import { DivTable, DivTableHead } from "src/views/WheresMyGeneV2/common/styles";
import { ROUTES } from "src/common/constants/routes";
import { useConnect } from "./connect";
import { SOURCE_DATA_INFO_TEXT } from "./constants";
import {
  StyledTooltip,
  TooltipContent,
  TooltipButton,
} from "src/views/WheresMyGeneV2/components/CellInfoSideBar/style";

export default function SourceData(): JSX.Element {
  const { collections, setHoverStartTime, handleBadgeHoverEnd } = useConnect();

  return (
    <Wrapper>
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
        <DivTable data-testid="source-data-list">
          <DivTableHead>
            <StyledDivTableCell>Collection</StyledDivTableCell>
            <StyledDivTableCell align>Datasets</StyledDivTableCell>
          </DivTableHead>
          {Object.values(collections).map(({ name, url, datasets }) => (
            <StyledDivTableRow key={name}>
              <StyledDivTableCell>
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
        </DivTable>
      </Content>
    </Wrapper>
  );
}
