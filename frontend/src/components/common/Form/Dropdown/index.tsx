import {
  DefaultDropdownMenuOption,
  Dropdown as SDSDropdown,
} from "@czi-sds/components";
import {
  FormLabelText,
  Optional,
} from "src/components/common/Form/common/style";
import {
  DropdownForm,
  DropdownPopper,
} from "src/components/common/Form/Dropdown/style";

type onChangeFn = (options: Value) => void;
export type Value =
  | DefaultDropdownMenuOption
  | DefaultDropdownMenuOption[]
  | null;

interface Props {
  disablePortal?: boolean;
  label: string;
  multiple?: boolean;
  onChange: onChangeFn;
  optionalField?: boolean;
  options: DefaultDropdownMenuOption[];
  text?: string;
  value: Value;
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
      {text && (
        <FormLabelText>
          <span>{text}</span>
          {optionalField && <Optional>(optional)</Optional>}
        </FormLabelText>
      )}
      <SDSDropdown
        closeOnBlur
        DropdownMenuProps={DropdownMenuProps}
        InputDropdownProps={{ sdsStyle: "square" }}
        isTriggerChangeOnOptionClick
        label={getDropdownLabel(value, label)}
        multiple={multiple}
        onChange={onChange}
        options={options}
        PopperComponent={({ ...props }) => (
          <DropdownPopper
            disablePortal={disablePortal}
            open={props.open}
            placement="bottom-start"
            style={{
              width: props.anchorEl?.offsetWidth, // Set Popper width to equal the dropdown button.
            }}
            {...props}
          />
        )}
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
function getDropdownLabel(value: Value, label: string): string {
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
function isOptionSelected(value: Value): boolean {
  if (!value) {
    return false;
  }
  if (Array.isArray(value)) {
    return value.length > 0;
  }
  return !!value;
}
