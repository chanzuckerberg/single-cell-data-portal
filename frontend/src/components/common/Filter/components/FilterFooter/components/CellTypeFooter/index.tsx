import { Tooltip } from "@czi-sds/components";
import {
  StyledTooltip,
  TooltipContent,
  TooltipLink,
  FooterContent,
  TooltipTrigger,
  LinkButton,
} from "../../style";
import { LINKS } from "src/common/constants/links";

export function CellTypeFooter() {
  const openOrganismFilter = () => {
    // TODO: Implement openOrganismFilter function
  };
  return (
    <FooterContent>
      <p>
        This filter includes cell types from multiple organisms which may
        reference{" "}
        <Tooltip
          sdsStyle="dark"
          placement="top"
          arrow
          slotProps={{
            tooltip: {
              style: {
                maxWidth: 395, // Override the max-width specification for dark sdsStyle.
              },
            },
          }}
          title={
            <StyledTooltip>
              <TooltipContent>
                Datasets for most organisms use the{" "}
                <TooltipLink
                  href={LINKS.CL_ONTOLOGY_LINK}
                  rel="noopener"
                  target="_blank"
                >
                  CL ontology
                </TooltipLink>{" "}
                to annotate cell types. Contributors may use organism-specific
                ontologies for:
                <ul>
                  <li>Caenorhabditis elegans</li>
                  <li>Danio rerio</li>
                  <li>Drosophila melanogaster</li>
                </ul>
              </TooltipContent>
            </StyledTooltip>
          }
        >
          <TooltipTrigger>organism-specific ontologies</TooltipTrigger>
        </Tooltip>
        . Values may be filtered with the{" "}
        <LinkButton onClick={openOrganismFilter}>organism filter</LinkButton>.
      </p>
    </FooterContent>
  );
}
