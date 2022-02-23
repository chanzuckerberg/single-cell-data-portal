import {
  Intent,
  ITextAreaProps,
  TextArea as RawTextArea,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import noop from "lodash/noop";
import { FC, useRef, useState } from "react";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import {
  FormLabelText as StyledLabelText,
  StyledFormLabel,
} from "src/components/common/Form/common/style";
import { LabelText, StyledLabel } from "src/components/common/Form/Input/style";
import { Value } from "../common/constants";
import { StyledDangerIcon, StyledIcon, Wrapper } from "./style";

/**
 * @deprecated - supersede by ErrorIcon once filter feature flag is removed (#1718).
 */
const DangerIcon = () => {
  return <StyledIcon icon={IconNames.ISSUE} intent={Intent.DANGER} />;
};

const ErrorIcon = () => {
  return <StyledDangerIcon icon={IconNames.ISSUE} intent={Intent.DANGER} />;
};

interface Props extends ITextAreaProps {
  isValid?: boolean;
  syncValidation?: Array<(value: string) => boolean>;
  handleChange?: ({ isValid, value, name }: Value) => void;
}

const TextArea: FC<Props> = (props) => {
  const [isValid, setIsValid] = useState(true);
  const inputRef = useRef<HTMLTextAreaElement>(null);

  /* Temporary structure to enable both filter feature and existing functionality. */
  /* Fragments and any deprecated styled components can either be removed or replaced once filter feature flag is removed (#1718). */
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const FormLabel = isFilterEnabled ? StyledFormLabel : StyledLabel;
  const FormLabelText = isFilterEnabled ? StyledLabelText : LabelText;
  const FormIcon = isFilterEnabled ? ErrorIcon : DangerIcon;

  const { handleChange = noop, ...restProps } = props;

  const { syncValidation = [], name = "default-text-area" } = props;

  return (
    <FormLabel htmlFor={name}>
      <FormLabelText>Description</FormLabelText>
      <Wrapper>
        <RawTextArea
          growVertically
          inputRef={inputRef}
          intent={(!isValid && Intent.DANGER) || undefined}
          {...restProps}
          onChange={handleChange_}
        />
        {!isValid && <FormIcon />}
      </Wrapper>
    </FormLabel>
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
