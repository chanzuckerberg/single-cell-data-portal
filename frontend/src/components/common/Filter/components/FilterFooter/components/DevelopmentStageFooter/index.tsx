import { TooltipLink, FooterContent } from "../../style";
import { FilterFooterTooltip } from "../components/FilterFooterTooltip";
import { OpenFilterButton } from "../components/OpenFilterButton";

import { LINKS } from "src/common/constants/links";

export function DevelopmentStageFooter() {
  return (
    <FooterContent>
      By default, this filter shows <span className="italic">Homo sapiens</span>{" "}
      and <span className="italic">Mus musculus</span>.{" "}
      <OpenFilterButton>Filter by organism</OpenFilterButton> to show other{" "}
      <FilterFooterTooltip
        content={
          <>
            Datasets for most organisms use the{" "}
            <TooltipLink
              href={LINKS.CL_ONTOLOGY_LINK}
              rel="noopener"
              target="_blank"
            >
              UBERON ontology
            </TooltipLink>{" "}
            to annotate developmental stage. Contributors may use
            organism-specific ontologies for:
            <ul>
              <li>Homo sapiens</li>
              <li>Mus musculus</li>
              <li>Caenorhabditis elegans</li>
              <li>Danio rerio</li>
              <li>Drosophila melanogaster</li>
            </ul>
          </>
        }
      >
        organism-specific developmental stages.
      </FilterFooterTooltip>
    </FooterContent>
  );
}
