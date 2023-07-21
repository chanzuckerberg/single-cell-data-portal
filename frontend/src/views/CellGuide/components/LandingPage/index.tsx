import Head from "next/head";
import CellGuideCardSearchBar from "../CellGuideCardSearchBar";
import { StyledHeader, Wrapper } from "./style";

export const LANDING_PAGE_HEADER = "landing-page-header";

const TITLE = "CellGuide Cell Types and Cell Tissues - CZ CELLxGENE";
const DESCRIPTION =
  "Explore single-cell transcriptomics data in CellGuide, a comprehensive resource that empowers researchers with deep insights into the intricacies of cell types";

export default function LandingPage(): JSX.Element {
  return (
    <>
      <Head>
        <title>{TITLE}</title>
        <meta property="title" key="title" content={TITLE} />
        <meta property="og:title" key="og:title" content={TITLE} />
        <meta property="twitter:title" key="twitter:title" content={TITLE} />

        <meta name="description" key="description" content={DESCRIPTION} />
        <meta
          property="og:description"
          key="og:description"
          content={DESCRIPTION}
        />
        <meta
          property="twitter:description"
          key="twitter:description"
          content={DESCRIPTION}
        />
      </Head>
      <Wrapper>
        <StyledHeader data-testid={LANDING_PAGE_HEADER}>
          CellGuide is a comprehensive resource for knowledge about cell types
        </StyledHeader>
        <CellGuideCardSearchBar />
      </Wrapper>
    </>
  );
}
