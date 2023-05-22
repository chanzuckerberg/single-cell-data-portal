import CellCardSearchBar from "../CellCardSearchBar";
import { StyledHeader, Wrapper } from "./style";

export const LANDING_PAGE_HEADER = "landing-page-header";

export default function LandingPage(): JSX.Element {
  return (
    <Wrapper>
      <StyledHeader data-testid={LANDING_PAGE_HEADER}>
        CellGuide is a comprehensive resource for knowledge about cell types
      </StyledHeader>
      <CellCardSearchBar />
    </Wrapper>
  );
}
