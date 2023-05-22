import CellCardSearchBar from "../CellCardSearchBar";
import { StyledHeader, Wrapper } from "./style";

export default function LandingPage(): JSX.Element {
  return (
    <Wrapper>
      <StyledHeader data-testid="landing-page-header">
        CellGuide is a comprehensive resource for knowledge about cell types
      </StyledHeader>
      <CellCardSearchBar />
    </Wrapper>
  );
}
