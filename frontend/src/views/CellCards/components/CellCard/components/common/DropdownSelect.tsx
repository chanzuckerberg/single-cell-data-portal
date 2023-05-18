import Select, { SelectChangeEvent } from "@mui/material/Select";
import FormControl from "@mui/material/FormControl";
import { MenuItem } from "czifui";
interface Props {
  handleChange: (event: SelectChangeEvent) => void;
  options: string[];
  selectedOption: string;
}

const DropdownSelect = ({ handleChange, options, selectedOption }: Props) => {
  return (
    <FormControl size="small">
      <Select
        labelId="dropdown-label"
        id="dropdown"
        value={selectedOption}
        onChange={handleChange}
      >
        {options.map((option) => (
          <MenuItem key={option} value={option}>
            {option}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  );
};
export default DropdownSelect;
