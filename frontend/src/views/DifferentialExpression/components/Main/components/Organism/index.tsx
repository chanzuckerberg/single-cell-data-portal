import {
  DefaultMenuSelectOption,
  InputDropdownProps as RawInputDropdownProps,
} from "@czi-sds/components";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { StyledDropdown, Wrapper, Label } from "../common/style";
import { DIFFERENTIAL_EXPRESSION_ORGANISM_DROPDOWN } from "src/views/DifferentialExpression/common/constants";
import { useConnect } from "./connect";

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

export default function Organism(): JSX.Element {
  const { organism, organisms, handleOnChange } = useConnect();
  return (
    <Wrapper>
      <Label>Organism</Label>
      <StyledDropdown
        label={organism?.name || "Select"}
        options={organisms || EMPTY_ARRAY}
        onChange={handleOnChange as tempOnChange}
        InputDropdownProps={InputDropdownProps}
        data-testid={DIFFERENTIAL_EXPRESSION_ORGANISM_DROPDOWN}
        value={organism}
      />
    </Wrapper>
  );
}

// (HACK): Not sure why styled Dropdown changes `onChange` type
type tempOnChange = (
  options: DefaultMenuSelectOption | DefaultMenuSelectOption[] | null
) => void;
