import {
  DefaultMenuSelectOption,
  InputDropdownProps as RawInputDropdownProps,
} from "czifui";
import { useContext, useEffect, useMemo } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { usePrimaryFilterDimensions } from "src/common/queries/wheresMyGene";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectOrganism } from "src/views/WheresMyGene/common/store/actions";
import { Organism as IOrganism } from "src/views/WheresMyGene/common/types";
import { Label } from "../../style";
import { StyledDropdown, Wrapper } from "./style";

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

interface Props {
  isLoading: boolean;
}

export default function Organism({ isLoading }: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { selectedOrganismId } = useContext(StateContext);
  const { data } = usePrimaryFilterDimensions();
  const { organisms } = data || {};

  useEffect(() => {
    if (!organisms || !dispatch) return;

    const organism = organisms.find(
      (organism: IOrganism) => organism.name === "Homo sapiens"
    );

    if (!organism) return;

    dispatch(selectOrganism(organism.id));
  }, [organisms, dispatch]);

  const organismsById = useMemo(() => {
    const result: { [id: string]: IOrganism } = {};

    if (!organisms) return result;

    for (const organism of organisms) {
      result[organism.id] = organism;
    }

    return result;
  }, [organisms]);

  const organism = organismsById[selectedOrganismId || ""];

  return (
    <Wrapper>
      <Label>Organism</Label>
      <StyledDropdown
        label={organism?.name || "Select"}
        options={organisms || EMPTY_ARRAY}
        onChange={handleOnChange as tempOnChange}
        InputDropdownProps={{ ...InputDropdownProps, disabled: isLoading }}
        data-test-id="add-organism"
      />
    </Wrapper>
  );

  function handleOnChange(organism: IOrganism | null): void {
    if (!dispatch || !organism) return;

    track(EVENTS.WMG_SELECT_ORGANISM, { payload: organism?.name });

    dispatch(selectOrganism(organism?.id || null));
  }
}

// (HACK): Not sure why styled Dropdown changes `onChange` type
type tempOnChange = (
  options: DefaultMenuSelectOption | DefaultMenuSelectOption[] | null
) => void;
