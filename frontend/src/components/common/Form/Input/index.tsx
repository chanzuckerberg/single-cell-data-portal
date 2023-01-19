import { FormGroup, Icon, InputGroup, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import debounce from "lodash/debounce";
import {
  FC,
  ReactElement,
  useCallback,
  useEffect,
  useRef,
  useState,
} from "react";
import {
  FormLabelText,
  Optional,
  StyledFormLabel,
} from "src/components/common/Form/common/style";
import { DEBOUNCE_TIME_MS } from "../../../CreateCollectionModal/components/Content/common/constants";
import { Value } from "../common/constants";

interface Props {
  markAsError?: boolean; // True if input has server side errors
  name: string;
  text: string;
  disabled?: boolean;
  handleChange: ({ isValid, value, name }: Value) => void;
  isRevalidationRequired?: boolean;
  leftElement?: ReactElement;
  syncValidation?: Array<(value: string) => true | string>;
  noNameAttr?: boolean;
  placeholder?: string;
  defaultValue?: string;
  optionalField?: boolean;
}

const ErrorIcon = (): JSX.Element => {
  return <Icon icon={IconNames.ISSUE} intent={Intent.DANGER} />;
};

const Input: FC<Props> = ({
  markAsError,
  name,
  text,
  disabled,
  handleChange,
  isRevalidationRequired = false,
  leftElement,
  syncValidation = [],
  noNameAttr = false,
  placeholder = "",
  defaultValue,
  optionalField = false,
}) => {
  const [isValid, setIsValid] = useState(true);
  const [errors, setErrors] = useState<string[]>([]);
  const inputRef = useRef<HTMLInputElement>(null);

  const handleChange_ = useCallback(() => {
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
  }, [handleChange, isValid, name, syncValidation]);

  // Revalidation is necessary if the link type has changed for this link to or from a DOI link type (as DOIs have
  // different validation to the other link types). Revalidation is not required when switching between link types
  // other than DOI as they share the same validation.
  useEffect(() => {
    if (isRevalidationRequired) {
      handleChange_();
    }
  }, [handleChange_, isRevalidationRequired]);

  // Check for server-side error messages and update error state of input if necessary.
  useEffect(() => {
    if (markAsError) {
      setIsValid(false);
    }
  }, [markAsError]);

  return (
    <StyledFormLabel htmlFor={name}>
      <FormGroup
        helperText={errors.join(", ")}
        intent={(!isValid && Intent.DANGER) || undefined}
      >
        <FormLabelText>
          <span>{text}</span>
          {optionalField && <Optional>(optional)</Optional>}
        </FormLabelText>
        <InputGroup
          // (thuang): `autoComplete="off"` and `type="search"` are needed to stop autofill
          // https://stackoverflow.com/a/30873633
          autoComplete="off"
          type="search"
          inputRef={inputRef}
          intent={(!isValid && Intent.DANGER) || undefined}
          id={name}
          leftElement={leftElement}
          name={noNameAttr ? undefined : name}
          rightElement={(!isValid && <ErrorIcon />) || undefined}
          disabled={disabled}
          onChange={debounce(handleChange_, DEBOUNCE_TIME_MS)}
          placeholder={placeholder}
          defaultValue={defaultValue}
        />
      </FormGroup>
    </StyledFormLabel>
  );
};

export default Input;
