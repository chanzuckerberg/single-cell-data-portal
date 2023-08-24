import {
  Label,
  StyledOntologyId,
  Wrapper,
} from "src/components/OntologyId/style";
import { SourceLink } from "src/views/CellGuide/components/CellGuideCard/components/Description/style";
import Link from "src/views/CellGuide/components/CellGuideCard/components/common/Link";

interface Props {
  ontologyId: string;
  cellTypeIdRaw: string | string[] | undefined;
}

export default function OntologyId({
  ontologyId = "",
  cellTypeIdRaw = "",
  ...rest
}: Props) {
  return (
    <Wrapper {...rest}>
      <Label>Ontology ID</Label>
      <StyledOntologyId>
        <SourceLink>
          <Link
            url={`https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F${cellTypeIdRaw}`}
            label={ontologyId}
          />
        </SourceLink>
      </StyledOntologyId>
    </Wrapper>
  );
}
