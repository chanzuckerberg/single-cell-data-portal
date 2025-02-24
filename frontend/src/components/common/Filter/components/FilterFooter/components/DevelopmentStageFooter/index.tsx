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

export function DevelopmentStageFooter() {
  const openOrganismFilter = () => {
    // TODO: Implement openOrganismFilter function
  };
  return (
    <FooterContent>
      <p>
        By default, this filter shows Homo sapiens and Mus musculus.{" "}
        <LinkButton onClick={openOrganismFilter}>Filter by organism</LinkButton>{" "}
        to show other{" "}
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
              </TooltipContent>
            </StyledTooltip>
          }
        >
          <TooltipTrigger>
            organism-specific developmental stages.
          </TooltipTrigger>
        </Tooltip>
      </p>
    </FooterContent>
  );
}
