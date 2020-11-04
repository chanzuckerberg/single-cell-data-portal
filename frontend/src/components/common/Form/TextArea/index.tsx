import {
  Intent,
  ITextAreaProps,
  TextArea as RawTextArea,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { noop } from "lodash-es";
import React, { FC, useRef, useState } from "react";
import { Value } from "../common/constants";
import { StyledIcon, Wrapper } from "./style";

const DangerIcon = () => {
  return <StyledIcon icon={IconNames.ISSUE} intent={Intent.DANGER} />;
};

interface Props extends ITextAreaProps {
  isValid?: boolean;
  syncValidation?: Array<(value: string) => boolean>;
  handleChange: ({ isValid, value, name }: Value) => void;
}

const TextArea: FC<Props> = (props) => {
  const [isValid, setIsValid] = useState(true);
  const inputRef = useRef<HTMLTextAreaElement>(null);

  const { handleChange = noop, ...restProps } = props;

  const { syncValidation = [], name = "default-text-area" } = props;

  return (
    <Wrapper>
      <RawTextArea
        inputRef={inputRef}
        intent={(!isValid && Intent.DANGER) || undefined}
        {...restProps}
        onChange={handleChange_}
      />
      {!isValid && <DangerIcon />}
    </Wrapper>
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
