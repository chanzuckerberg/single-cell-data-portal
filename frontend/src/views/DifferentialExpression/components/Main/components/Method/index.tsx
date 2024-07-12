import { InputDropdownProps as RawInputDropdownProps } from "@czi-sds/components";

import { EMPTY_ARRAY, noop } from "src/common/constants/utils";

import { StyledDropdown, Label } from "../common/style";
import { StyledP, StyledTooltipText, Wrapper } from "./style";
import { DIFFERENTIAL_EXPRESSION_METHOD_DROPDOWN } from "src/views/DifferentialExpression/common/constants";
import HelpTooltip from "src/views/CellGuide/components/CellGuideCard/components/common/HelpTooltip";

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

const METHODS = [{ name: "T-test" }];

export default function Method(): JSX.Element {
  const method = METHODS[0];
  return (
    <Wrapper>
      <Label>
        Method{" "}
        <HelpTooltip
          text={
            <StyledTooltipText>
              <StyledP>Currently, this tool only supports T-test.</StyledP>
              <p>
                While the t-test performs well on individual datasets{" "}
                <a
                  href="https://www.mdpi.com/2073-4425/12/12/1947"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  [1]
                </a>
                <a
                  href="https://www.nature.com/articles/s41467-019-12266-7"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  [2]
                </a>
                <a
                  href="https://www.nature.com/articles/nmeth.4612"
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  [3]
                </a>
                , its performance on concatenated, non-integrated datasets has
                not been extensively evaluated.
              </p>
              <StyledP>More methods will be added soon.</StyledP>
            </StyledTooltipText>
          }
          dark
          placement="top"
        />
      </Label>
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
