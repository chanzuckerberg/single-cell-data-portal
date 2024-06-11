import { InputDropdownProps as RawInputDropdownProps } from "@czi-sds/components";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { StyledDropdown, Wrapper, Label } from "../common/style";
import { DIFFERENTIAL_EXPRESSION_ORGANISM_DROPDOWN } from "src/views/DifferentialExpression/common/constants";
import { useConnect } from "./connect";
import { Organism as IOrganism } from "src/views/DifferentialExpression/common/types";

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

export default function Organism(): JSX.Element {
  const { organism, organisms, handleOnChange } = useConnect();
  return (
    <Wrapper>
      <Label>Organism</Label>
      <StyledDropdown<IOrganism, false, false, false>
        label={organism?.name || "Select"}
        options={organisms || EMPTY_ARRAY}
        onChange={handleOnChange}
        InputDropdownProps={InputDropdownProps}
        data-testid={DIFFERENTIAL_EXPRESSION_ORGANISM_DROPDOWN}
        value={organism}
      />
    </Wrapper>
  );
}
