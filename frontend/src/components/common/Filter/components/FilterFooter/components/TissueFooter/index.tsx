import { TooltipLink, FooterContent } from "../../style";
import { LINKS } from "src/common/constants/links";
import { FilterFooterTooltip } from "../components/FilterFooterTooltip";
import { OpenFilterButton } from "../components/OpenFilterButton";

export function TissueFooter() {
  return (
    <FooterContent>
      This filter includes tissues from multiple organisms which may reference{" "}
      <FilterFooterTooltip
        content={
          <>
            Datasets for most organisms use the{" "}
            <TooltipLink
              href={LINKS.UBERON_ONTOLOGY_LINK}
              rel="noopener"
              target="_blank"
            >
              UBERON ontology
            </TooltipLink>{" "}
            to annotate tissues. Contributors may use organism-specific
            ontologies for:
            <ul>
              <li>Caenorhabditis elegans</li>
              <li>Danio rerio</li>
              <li>Drosophila melanogaster</li>
            </ul>
          </>
        }
      >
        organism-specific ontologies
      </FilterFooterTooltip>
      . Values may be filtered with the{" "}
      <OpenFilterButton>organism filter</OpenFilterButton>.
    </FooterContent>
  );
}
