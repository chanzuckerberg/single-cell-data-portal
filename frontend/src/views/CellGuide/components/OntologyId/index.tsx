import {
  Label,
  StyledOntologyId,
  Wrapper,
} from "src/views/CellGuide/components/OntologyId/style";
import { SourceLink } from "src/views/CellGuide/components/CellGuideCard/components/Description/style";
import Link from "src/views/CellGuide/components/CellGuideCard/components/common/Link";

interface Props {
  ontologyId: string;
  url: string;
}

export default function OntologyId({
  ontologyId = "",
  url = "",
  ...rest
}: Props) {
  return (
    <Wrapper {...rest}>
      <Label>Ontology ID</Label>
      <StyledOntologyId>
        <SourceLink>
          <Link url={url} label={ontologyId} />
        </SourceLink>
      </StyledOntologyId>
    </Wrapper>
  );
}
