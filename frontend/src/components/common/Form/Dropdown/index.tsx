import {
  DefaultMenuSelectOption,
  Dropdown as SDSDropdown,
  InputDropdownProps,
} from "czifui";
import {
  FormLabelText,
  Optional,
} from "src/components/common/Form/common/style";
import {
  DropdownForm,
  DropdownPopper,
  MAX_DISPLAYABLE_MENU_ITEMS,
} from "src/components/common/Form/Dropdown/style";

type onChangeFn = (value: Value) => void;
export type Value = DefaultMenuSelectOption | DefaultMenuSelectOption[] | null;

interface Props {
  label: string;
  multiple?: boolean;
  onChange: onChangeFn;
  optionalField?: boolean;
  options: DefaultMenuSelectOption[];
  text?: string;
  value: Value;
}

const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<InputDropdownProps> = {
  sdsStage: "default",
  sdsStyle: "square",
  sdsType: "multiSelect",
};

export default function Dropdown({
  label,
  multiple = false,
  onChange,
  optionalField = false,
  options,
  text,
  value,
}: Props): JSX.Element {
  return (
    <DropdownForm>
      {text && (
        <FormLabelText>
          <span>{text}</span>
          {optionalField && <Optional>(optional)</Optional>}
        </FormLabelText>
      )}
      <SDSDropdown
        closeOnBlur
        InputDropdownProps={{ ...DEFAULT_INPUT_DROPDOWN_PROPS }}
        label={getDropdownLabel(value, label)}
        MenuSelectProps={{
          getOptionSelected,
        }}
        multiple={multiple}
        onChange={onChange}
        options={options}
        PopperComponent={({ ...props }) => (
          <DropdownPopper
            children={props.children}
            // Disabling portal is required for successful rendering of the menu. Blueprint Dialog prop "enforceFocus"
            // default is "true" preventing focus from leaving itself and so disabling portal for the Popper component
            // renders the component under the DOM hierarchy of the Dropdown component and therefore the Dialog component.
            disablePortal
            isScrollable={options.length > MAX_DISPLAYABLE_MENU_ITEMS}
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
 * Used to determine if an option is selected, considering the current value.
 * @param option - The option to test.
 * @param value - The value to test against.
 * @returns true if the specified option is selected.
 */
function getOptionSelected(
  option: { name: string },
  value: { name: string }
): boolean {
  return option.name === value.name;
}
