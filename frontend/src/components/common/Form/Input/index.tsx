import { FormGroup, InputGroup, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import debounce from "lodash/debounce";
import React, { FC, useRef, useState } from "react";
import { DEBOUNCE_TIME_MS } from "../../../CreateCollectionModal/components/Content/common/constants";
import { Value } from "../common/constants";
import { LabelText, StyledIcon, StyledLabel } from "./style";

interface Props {
  name: string;
  text: string;
  percentage?: number;
  disabled?: boolean;
  handleChange: ({ isValid, value, name }: Value) => void;
  syncValidation?: Array<(value: string) => true | string>;
  noNameAttr?: boolean;
  className?: string;
  placeholder?: string;
}

const DangerIcon = () => {
  return <StyledIcon icon={IconNames.ISSUE} intent={Intent.DANGER} />;
};

const Input: FC<Props> = ({
  name,
  text,
  percentage = 100,
  disabled,
  handleChange,
  syncValidation = [],
  noNameAttr = false,
  className,
  placeholder = "",
}) => {
  const [isValid, setIsValid] = useState(true);
  const [errors, setErrors] = useState<string[]>([]);
  const inputRef = useRef<HTMLInputElement>(null);

  return (
    <StyledLabel percentage={percentage} htmlFor={name} className={className}>
      <FormGroup
        helperText={errors.join(", ")}
        intent={(!isValid && Intent.DANGER) || undefined}
      >
        <LabelText>{text}</LabelText>
        <InputGroup
          inputRef={inputRef}
          intent={(!isValid && Intent.DANGER) || undefined}
          id={name}
          name={noNameAttr ? undefined : name}
          rightElement={(!isValid && <DangerIcon />) || undefined}
          disabled={disabled}
          onChange={debounce(handleChange_, DEBOUNCE_TIME_MS)}
          placeholder={placeholder}
        />
      </FormGroup>
    </StyledLabel>
  );

  function handleChange_() {
    if (!inputRef.current) return;

    const value = inputRef.current.value;

    const validation = syncValidation.map((validate) => validate(value));

    const isSyncValid = validation.every((result) => result === true);

    const errors = validation.filter((error) => error !== true) as string[];

    const result = Boolean(value?.length) && isSyncValid;

    handleChange({ isValid: result, name, value });

    setErrors(errors);

    if (result !== isValid) {
      setIsValid(result);
    }
  }
};

export default Input;
