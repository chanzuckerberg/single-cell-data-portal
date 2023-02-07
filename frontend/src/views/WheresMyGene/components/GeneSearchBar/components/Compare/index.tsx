import {
  DefaultDropdownMenuOption,
  DefaultMenuSelectOption,
  InputDropdownProps as RawInputDropdownProps,
} from "czifui";
import { useContext } from "react";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { FILTER_LABELS } from "src/views/WheresMyGene/common/constants";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectCompare } from "src/views/WheresMyGene/common/store/actions";
import { CompareDimensionOption } from "../../../../common/types";
import { CompareLabel, StyledDropdown, Wrapper } from "./style";

const InputDropdownProps: Partial<RawInputDropdownProps> = {
  sdsStyle: "square",
};

interface Props {
  isLoading: boolean;
}

const FILTERS: CompareDimensionOption[] = [
  {
    name: "None",
    id: "",
  },
  // {
  //   name: FILTER_LABELS.DATASET,
  //   id: "dataset"
  // },
  {
    name: FILTER_LABELS.DISEASE,
    id: "disease",
  },
  {
    name: FILTER_LABELS.ETHNICITY,
    id: "ethnicity",
  },
  {
    name: FILTER_LABELS.SEX,
    id: "gender",
  },
];

export default function Compare({ isLoading }: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { selectedCompare } = useContext(StateContext);

  const dimension = FILTERS.find((filter) => filter.id === selectedCompare);
  return (
    <Wrapper>
      <CompareLabel>Compare</CompareLabel>
      <StyledDropdown
        label={dimension?.name || "None"}
        options={
          (FILTERS as unknown as DefaultDropdownMenuOption[]) || EMPTY_ARRAY
        }
        onChange={handleOnChange as tempOnChange}
        InputDropdownProps={{ ...InputDropdownProps, disabled: isLoading }}
        data-test-id="compare-dropdown"
      />
    </Wrapper>
  );

  function handleOnChange(dimension: CompareDimensionOption | null): void {
    if (!dispatch || !dimension) return;

    console.log("Selected", dimension.id);

    // track(EVENTS.WMG_SELECT_ORGANISM, { payload: organism?.name });

    dispatch(selectCompare(dimension.id || null));
  }
}

// (HACK): Not sure why styled Dropdown changes `onChange` type
type tempOnChange = (
  options: DefaultMenuSelectOption | DefaultMenuSelectOption[] | null
) => void;
