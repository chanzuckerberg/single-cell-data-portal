import {
  DefaultAutocompleteOption,
  DropdownProps,
  Dropdown as SDSDropdown,
} from "@czi-sds/components";
import { AutocompleteValue } from "@mui/base";
import {
  FormLabelText,
  Optional,
} from "src/components/common/Form/common/style";
import {
  DropdownForm,
  DropdownPopper,
} from "src/components/common/Form/Dropdown/style";

export type Value = AutocompleteValue<
  DefaultAutocompleteOption,
  boolean | undefined,
  false,
  false
>;

interface Props
  extends DropdownProps<
    DefaultAutocompleteOption,
    boolean | undefined,
    false,
    false
  > {
  disablePortal?: boolean;
  label: string;
  optionalField?: boolean;
  text?: string;
}

const DropdownMenuProps = {
  PopperBaseProps: {
    modifiers: [
      {
        name: "offset",
        options: {
          offset: [0, 0],
        },
      },
    ],
  },
};

export default function Dropdown({
  disablePortal = false,
  label,
  multiple = false,
  onChange,
  optionalField = false,
  options,
  text,
  value,
}: Props): JSX.Element {
  return (
    <DropdownForm isSelected={isOptionSelected(value)}>
      {!!text && (
        <FormLabelText>
          <span>{text}</span>
          {optionalField && <Optional>(optional)</Optional>}
        </FormLabelText>
      )}
      <SDSDropdown<DefaultAutocompleteOption, boolean | undefined, false, false>
        closeOnBlur
        DropdownMenuProps={DropdownMenuProps}
        InputDropdownProps={{ sdsStyle: "square" }}
        isTriggerChangeOnOptionClick
        label={getDropdownLabel(label, value)}
        multiple={multiple}
        onChange={onChange} // TODO(SDSv20): Come back to this
        options={options}
        PopperComponent={(popperProps) => {
          const { anchorEl } = popperProps;

          return (
            <DropdownPopper
              disablePortal={disablePortal}
              {...popperProps}
              style={{ width: (anchorEl as HTMLButtonElement)?.offsetWidth }}
            />
          );
        }}
        value={value}
      />
    </DropdownForm>
  );
}

/**
 * Returns the label of the dropdown as either a placeholder if no value is selected, or
 * the selected value or (concatenated) values as a displayable string.
 * @param value - Dropdown value or values.
 * @param label - Default "placeholder" label.
 * @returns dropdown label.
 */
function getDropdownLabel(label: string, value?: Value): string {
  if (!value || (Array.isArray(value) && value.length === 0)) {
    return label;
  }
  return Array.isArray(value)
    ? value?.map(({ name }) => name).join(", ")
    : value?.name;
}

/**
 * Returns true if an option is selected.
 * @param value - Selected value.
 * @returns true is an option is selected.
 */
function isOptionSelected(value?: Value): boolean {
  if (!value) {
    return false;
  }
  if (Array.isArray(value)) {
    return value.length > 0;
  }
  return !!value;
}
