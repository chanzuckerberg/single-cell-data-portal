import { Icon, InputGroup, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { FC, useState } from "react";
import { LabelText, StyledLabel } from "./style";

interface Props {
  name: string;
  text: string;
  percentage?: number;
  disabled?: boolean;
  onChange: (event: React.ChangeEvent<HTMLInputElement>) => void;
  inputRef: any;
}

const DangerIcon = () => {
  return <Icon icon={IconNames.ISSUE} intent={Intent.DANGER} />;
};

const Input: FC<Props> = ({
  name,
  text,
  percentage = 100,
  disabled,
  onChange,
  inputRef,
}) => {
  const [isInvalid] = useState(false);

  return (
    <StyledLabel percentage={percentage} htmlFor={name}>
      <LabelText>{text}</LabelText>
      <InputGroup
        inputRef={inputRef}
        intent={(isInvalid && Intent.DANGER) || undefined}
        id={name}
        name={name}
        rightElement={(isInvalid && <DangerIcon />) || undefined}
        disabled={disabled}
        onChange={onChange}
      />
    </StyledLabel>
  );
};

export default Input;
