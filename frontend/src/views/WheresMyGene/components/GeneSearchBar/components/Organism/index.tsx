import { DefaultMenuSelectOption } from "czifui";
import { useReducer } from "react";
import { INITIAL_STATE, reducer } from "src/views/WheresMyGene/common/store";
import { selectOrganism } from "src/views/WheresMyGene/common/store/actions";
import { Organism as IOrganism } from "src/views/WheresMyGene/common/types";
import { StyledDropdown } from "./style";

const ORGANISMS: IOrganism[] = [
  { name: "Homo sapiens" },
  { name: "mus musculus" },
];

export default function Organism(): JSX.Element {
  const [state, dispatch] = useReducer(reducer, INITIAL_STATE);

  const { selectedOrganism } = state;

  return (
    <StyledDropdown
      label={selectedOrganism?.name || ""}
      options={ORGANISMS}
      multiple={false}
      onChange={handleOnChange}
    />
  );

  function handleOnChange(organism: DefaultMenuSelectOption | null): void {
    dispatch(selectOrganism(organism));
  }
}
