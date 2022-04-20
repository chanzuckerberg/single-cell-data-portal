import Link from "next/link";
import { EXTERNAL_LINKS, ROUTES } from "src/common/constants/routes";
import Step from "./components/Step";
import { Content, Details, Header, Step3Details } from "./style";

export default function GetStarted(): JSX.Element {
  return (
    <>
      <Header>Getting Started</Header>
      <Details>
        Visualize and compare the expression of genes across cell types using a
        dot plot. A dot plot is an effective way of summarizing expression data
        for genes across cells in categories. It gives an overview of the data
        by condensing many single cells into a mean expression (color scale) and
        percentage of cells (circle size) within each cell type that express
        each gene. Data are sourced from the{" "}
        <Link href={ROUTES.HOMEPAGE} passHref>
          <a href="passHref" rel="noopener" target="_blank">
            cellxgene Data Portal
          </a>
        </Link>{" "}
        data corpus.{" "}
      </Details>
      <Content>
        <Step
          step={1}
          header="Add Tissues"
          details="Add tissues you are interested in exploring. Cell types included in these tissues will automatically be added to the visualization. Cell types are defined by the annotations of the original authors."
        />
        <Step
          step={2}
          header="Add Genes"
          details="Add genes of interest to the visualization."
        />
        <Step
          step={3}
          header="Interpret Results"
          details={
            <>
              <Step3Details>
                Darker dots represent higher relative gene expression. The
                larger the dot, the more cells express that gene.
              </Step3Details>
              <a href={EXTERNAL_LINKS.WMG_DOC} rel="noopener" target="_blank">
                Learn More
              </a>
            </>
          }
        />
      </Content>
    </>
  );
}
