import {
  DefaultMenuSelectOption,
  InputDropdownProps as RawInputDropdownProps,
} from "@czi-sds/components";
import { useContext, useEffect, useMemo } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { useAvailableOrganisms } from "src/common/queries/differentialExpression";
import {
  DispatchContext,
  StateContext,
} from "src/views/DifferentialExpression/common/store";
import { selectOrganism } from "src/views/DifferentialExpression/common/store/actions";
import { StyledDropdown, Wrapper, Label } from "../common/style";
import { Organism as IOrganism } from "src/views/DifferentialExpression/common/types";

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

export default function Organism(): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { organismId } = useContext(StateContext);
  const { data: organisms } = useAvailableOrganisms();

  // (thuang): Default to "Homo sapiens" on first load
  useEffect(() => {
    if (!organisms || !dispatch || organismId) return;

    const organism = organisms.find(
      (organism: IOrganism) => organism.name === "Homo sapiens"
    );

    if (!organism) return;

    dispatch(selectOrganism(organism.id));
  }, [organisms, dispatch, organismId]);

  const organismsById = useMemo(() => {
    const result: { [id: string]: IOrganism } = {};

    if (!organisms) return result;

    for (const organism of organisms) {
      result[organism.id] = organism;
    }

    return result;
  }, [organisms]);

  const organism = organismsById[organismId || ""];
  return (
    <Wrapper>
      <Label>Organism</Label>
      <StyledDropdown
        label={organism?.name || "Select"}
        options={organisms || EMPTY_ARRAY}
        onChange={handleOnChange as tempOnChange}
        InputDropdownProps={InputDropdownProps}
        data-testid="add-organism-de"
        value={organism}
      />
    </Wrapper>
  );

  function handleOnChange(organism: IOrganism | null): void {
    if (!dispatch || !organism || organismId === organism.id) return;

    track(EVENTS.WMG_SELECT_ORGANISM, { payload: organism?.name });

    dispatch(selectOrganism(organism?.id || null));
  }
}

// (HACK): Not sure why styled Dropdown changes `onChange` type
type tempOnChange = (
  options: DefaultMenuSelectOption | DefaultMenuSelectOption[] | null
) => void;
