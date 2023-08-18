import { SelectChangeEvent } from "@mui/material/Select";
import FormControl from "@mui/material/FormControl";
import { MenuItem } from "@czi-sds/components";
import { StyledSelect } from "./style";
interface Props {
  handleChange: (event: SelectChangeEvent<unknown>) => void;
  options: string[];
  selectedOption: string;
  testId?: string;
}

const DropdownSelect = ({
  handleChange,
  options,
  selectedOption,
  testId,
}: Props) => {
  return (
    <FormControl size="small">
      <StyledSelect
        labelId="dropdown-label"
        id="dropdown"
        data-testid={testId}
        value={selectedOption}
        onChange={handleChange}
      >
        {options.map((option) => (
          <MenuItem key={option} value={option}>
            {option}
          </MenuItem>
        ))}
      </StyledSelect>
    </FormControl>
  );
};
export default DropdownSelect;
