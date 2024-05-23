import { InputDropdownProps as RawInputDropdownProps } from "@czi-sds/components";

import { EMPTY_ARRAY, noop } from "src/common/constants/utils";

import { StyledDropdown, Label } from "../common/style";
import { Wrapper } from "./style";
import { DIFFERENTIAL_EXPRESSION_METHOD_DROPDOWN } from "src/views/DifferentialExpression/common/constants";

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

const METHODS = [{ name: "T-test" }];

export default function Method(): JSX.Element {
  const method = METHODS[0];
  return (
    <Wrapper>
      <Label>Method</Label>
      <StyledDropdown
        label={method.name || "Select"}
        options={METHODS || EMPTY_ARRAY}
        InputDropdownProps={InputDropdownProps}
        data-testid={DIFFERENTIAL_EXPRESSION_METHOD_DROPDOWN}
        value={method}
        onChange={noop}
        disabled={true}
      />
    </Wrapper>
  );
}
