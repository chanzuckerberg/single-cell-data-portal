import {
  DefaultMenuSelectOption,
  InputDropdownProps as RawInputDropdownProps,
} from "czifui";
import { useContext } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectOrganism } from "src/views/WheresMyGene/common/store/actions";
import { Organism as IOrganism } from "src/views/WheresMyGene/common/types";
import { Label } from "../../style";
import { StyledDropdown, Wrapper } from "./style";

const ORGANISMS: { name: IOrganism }[] = [
  { name: "Homo sapiens" },
  { name: "mus musculus" },
];

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

export default function Organism(): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { selectedOrganism } = useContext(StateContext);

  return (
    <Wrapper>
      <Label>Organism</Label>
      <StyledDropdown
        label={selectedOrganism || "Select"}
        options={ORGANISMS}
        onChange={handleOnChange as tempOnChange}
        InputDropdownProps={InputDropdownProps}
        data-test-id="add-organism"
      />
    </Wrapper>
  );

  function handleOnChange(organism: DefaultMenuSelectOption | null): void {
    if (!dispatch) return;

    dispatch(selectOrganism(organism?.name || null));
  }
}

// (HACK): Not sure why styled Dropdown changes `onChange` type
type tempOnChange = (
  options: DefaultMenuSelectOption | DefaultMenuSelectOption[] | null
) => void;
