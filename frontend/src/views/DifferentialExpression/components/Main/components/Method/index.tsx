import { InputDropdownProps as RawInputDropdownProps } from "@czi-sds/components";

import { EMPTY_ARRAY, noop } from "src/common/constants/utils";

import { StyledDropdown, Label } from "../common/style";
import { FlexRow, Wrapper } from "./style";
import HelpTooltip from "src/views/CellGuide/components/CellGuideCard/components/common/HelpTooltip";

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

const METHODS = [{ name: "T-test" }];

export default function Method(): JSX.Element {
  const method = METHODS[0];
  return (
    <Wrapper>
      <FlexRow>
        <Label>Method</Label>
        <HelpTooltip
          title="Method"
          dark
          buttonDataTestId={"dummy-tooltip-test-id"}
          text={<>dummy text</>}
        />
      </FlexRow>
      <StyledDropdown
        label={method.name || "Select"}
        options={METHODS || EMPTY_ARRAY}
        InputDropdownProps={InputDropdownProps}
        data-testid="method-de"
        value={method}
        onChange={noop}
        disabled={true}
      />
    </Wrapper>
  );
}
