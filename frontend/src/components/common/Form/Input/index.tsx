import { FormGroup, Icon, InputGroup, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import debounce from "lodash/debounce";
import { FC, useEffect, useRef, useState } from "react";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import {
  FormLabelText as StyledLabelText,
  StyledFormLabel,
} from "src/components/common/Form/common/style";
import { DEBOUNCE_TIME_MS } from "../../../CreateCollectionModal/components/Content/common/constants";
import { Value } from "../common/constants";
import { LabelText, StyledIcon, StyledInputGroup, StyledLabel } from "./style";

interface Props {
  name: string;
  text: string;
  percentage?: number;
  disabled?: boolean;
  handleChange: ({ isValid, value, name }: Value) => void;
  isRevalidationRequired?: boolean;
  syncValidation?: Array<(value: string) => true | string>;
  noNameAttr?: boolean;
  className?: string;
  placeholder?: string;
  defaultValue?: string;
  optionalField?: boolean;
}

/**
 * @deprecated - supersede by ErrorIcon once filter feature flag is removed (#1718).
 */
const DangerIcon = (): JSX.Element => {
  return <StyledIcon icon={IconNames.ISSUE} intent={Intent.DANGER} />;
};

const ErrorIcon = (): JSX.Element => {
  return <Icon icon={IconNames.ISSUE} intent={Intent.DANGER} />;
};

const Input: FC<Props> = ({
  name,
  text,
  percentage = 100,
  disabled,
  handleChange,
  isRevalidationRequired = false,
  syncValidation = [],
  noNameAttr = false,
  className,
  placeholder = "",
  defaultValue,
  optionalField = false,
}) => {
  const [isValid, setIsValid] = useState(true);
  const [errors, setErrors] = useState<string[]>([]);
  const inputRef = useRef<HTMLInputElement>(null);

  /* Temporary structure to enable both filter feature and existing functionality. */
  /* Fragments and any deprecated styled components can either be removed or replaced once filter feature flag is removed (#1718). */
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const labelProps = isFilterEnabled
    ? undefined
    : {
        className: className,
        percentage: percentage,
      }; /* remove labelProps from <FormLabel/> once filter feature flag is removed (#1718). */
  const FormLabel = isFilterEnabled ? StyledFormLabel : StyledLabel;
  const FormLabelText = isFilterEnabled ? StyledLabelText : LabelText;
  const FormInputGroup = isFilterEnabled ? InputGroup : StyledInputGroup;
  const FormIcon = isFilterEnabled ? ErrorIcon : DangerIcon;

  // Revalidation is necessary if the link type has changed for this link to or from a DOI link type (as DOIs have
  // different validation to the other link types). Revalidation is not required when switching between link types
  // other than DOI as they share the same validation.
  useEffect(() => {
    if (isRevalidationRequired) {
      handleChange_();
    }
  }, [handleChange_, isRevalidationRequired]);

  return (
    <FormLabel htmlFor={name} {...labelProps}>
      <FormGroup
        helperText={errors.join(", ")}
        intent={(!isValid && Intent.DANGER) || undefined}
      >
        <FormLabelText>
          <span>{text}</span>
          {optionalField && <i>(optional)</i>}
        </FormLabelText>
        <FormInputGroup
          // (thuang): `autoComplete="off"` and `type="search"` are needed to stop autofill
          // https://stackoverflow.com/a/30873633
          autoComplete="off"
          type="search"
          inputRef={inputRef}
          intent={(!isValid && Intent.DANGER) || undefined}
          id={name}
          name={noNameAttr ? undefined : name}
          rightElement={(!isValid && <FormIcon />) || undefined}
          disabled={disabled}
          onChange={debounce(handleChange_, DEBOUNCE_TIME_MS)}
          placeholder={placeholder}
          defaultValue={defaultValue}
        />
      </FormGroup>
    </FormLabel>
  );

  function handleChange_() {
    if (!inputRef.current) return;

    const value = inputRef.current.value;

    const validation = syncValidation.map((validate) => validate(value));

    const isSyncValid = validation.every((result) => result === true);

    const errors = validation.filter((error) => error !== true) as string[];

    const result = isSyncValid;

    handleChange({ isValid: result, name, value });

    setErrors(errors);

    if (result !== isValid) {
      setIsValid(result);
    }
  }
};

export default Input;
