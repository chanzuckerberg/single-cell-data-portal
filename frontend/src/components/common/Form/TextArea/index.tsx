import {
  Intent,
  TextAreaProps,
  TextArea as RawTextArea,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import noop from "lodash/noop";
import { FC, useRef, useState } from "react";
import {
  FormLabelText,
  StyledFormLabel,
} from "src/components/common/Form/common/style";
import { Value } from "../common/constants";
import { StyledDangerIcon, Wrapper } from "./style";

const ErrorIcon = () => {
  return <StyledDangerIcon icon={IconNames.ISSUE} intent={Intent.DANGER} />;
};

interface Props extends TextAreaProps {
  isValid?: boolean;
  syncValidation?: Array<(value: string) => boolean>;
  handleChange?: ({ isValid, value, name }: Value) => void;
}

const TextArea: FC<Props> = (props) => {
  const [isValid, setIsValid] = useState(true);
  const inputRef = useRef<HTMLTextAreaElement>(null);

  const { handleChange = noop, ...restProps } = props;

  const { syncValidation = [], name = "default-text-area" } = props;

  return (
    <StyledFormLabel htmlFor={name}>
      <FormLabelText>Description</FormLabelText>
      <Wrapper>
        <RawTextArea
          growVertically
          inputRef={inputRef}
          intent={(!isValid && Intent.DANGER) || undefined}
          {...restProps}
          onChange={handleChange_}
        />
        {!isValid && <ErrorIcon />}
      </Wrapper>
    </StyledFormLabel>
  );

  function handleChange_() {
    if (!inputRef.current) return;

    const value = inputRef.current.value;

    const isSyncValid = syncValidation
      .map((validate) => validate(value))
      .every((result) => result);

    const result = Boolean(value?.length) && isSyncValid;

    handleChange({ isValid: result, name, value });

    if (result !== isValid) {
      setIsValid(result);
    }
  }
};

export default TextArea;
