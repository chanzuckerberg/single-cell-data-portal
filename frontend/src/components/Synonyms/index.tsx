import { EMPTY_ARRAY } from "src/common/constants/utils";
import { Label, Wrapper, Synonym } from "src/components/Synonyms/style";

interface Props {
  synonyms?: string[];
}

export default function Synonyms({ synonyms = EMPTY_ARRAY, ...rest }: Props) {
  return (
    <Wrapper {...rest}>
      <Label>Synonyms</Label>
      <Synonym>{synonyms.join(", ") || "N/A"}</Synonym>
    </Wrapper>
  );
}
