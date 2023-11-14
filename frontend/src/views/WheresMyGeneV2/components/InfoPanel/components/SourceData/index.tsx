import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import {
  Content,
  InfoText,
  StyledDivTableCell,
  StyledDivTableRow,
  StyledLabel,
  StyledLink,
  Wrapper,
} from "./style";
import { DivTable, DivTableHead } from "../../../CellInfoSideBar/style";
import { ROUTES } from "src/common/constants/routes";
import { useConnect } from "./connect";
import { SOURCE_DATA_INFO_TEXT } from "./constants";

export default function SourceData(): JSX.Element {
  const { collections } = useConnect();

  return (
    <Wrapper>
      <Content>
        <InfoText>
          {SOURCE_DATA_INFO_TEXT}{" "}
          <a
            href={ROUTES.WMG_DOCS_DATA_PROCESSING}
            target="_blank"
            rel="noreferrer noopener"
          >
            exceptions and processing notes
          </a>
          .
        </InfoText>
        <DivTable>
          <DivTableHead>
            <StyledDivTableCell>Collection</StyledDivTableCell>
            <StyledDivTableCell align>Datasets</StyledDivTableCell>
          </DivTableHead>
          {Object.values(collections).map(
            ({ name, url, datasets, total_count }) => (
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
                  <StyledLabel>{`${datasets.length} of ${total_count}`}</StyledLabel>
                </StyledDivTableCell>
              </StyledDivTableRow>
            )
          )}
        </DivTable>
      </Content>
    </Wrapper>
  );
}
