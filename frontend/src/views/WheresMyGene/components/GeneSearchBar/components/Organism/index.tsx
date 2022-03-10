import {
  DefaultMenuSelectOption,
  InputDropdownProps as RawInputDropdownProps,
} from "czifui";
import { useReducer } from "react";
import { INITIAL_STATE, reducer } from "src/views/WheresMyGene/common/store";
import { selectOrganism } from "src/views/WheresMyGene/common/store/actions";
import { Organism as IOrganism } from "src/views/WheresMyGene/common/types";
import { Label } from "../../style";
import { StyledDropdown, Wrapper } from "./style";

const ORGANISMS: IOrganism[] = [
  { name: "Homo sapiens" },
  { name: "mus musculus" },
];

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

export default function Organism(): JSX.Element {
  const [state, dispatch] = useReducer(reducer, INITIAL_STATE);

  const { selectedOrganism } = state;

  return (
    <Wrapper>
      <Label>Organism</Label>
      <StyledDropdown
        label={selectedOrganism?.name || ""}
        options={ORGANISMS}
        onChange={handleOnChange as tempOnChange}
        InputDropdownProps={InputDropdownProps}
      />
    </Wrapper>
  );

  function handleOnChange(organism: DefaultMenuSelectOption | null): void {
    dispatch(selectOrganism(organism));
  }
}

// (HACK): Not sure why styled Dropdown changes `onChange` type
type tempOnChange = (
  options: DefaultMenuSelectOption | DefaultMenuSelectOption[] | null
) => void;
