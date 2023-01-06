import {
  Icon,
  InputDropdownProps as IInputDropdownProps,
  Tooltip,
} from "czifui";
import { useContext, useMemo } from "react";
import {
  DispatchContext,
  StateContext,
} from "src/views/WheresMyGene/common/store";
import { selectSortBy } from "src/views/WheresMyGene/common/store/actions";
import { SORT_BY } from "src/views/WheresMyGene/common/types";
import { FlexDiv, Label, LabelWrapper, StyledDropdown, Wrapper } from "./style";

interface Props {
  handleIsScaledChange: () => void;
}

const COLOR_SCALE_OPTIONS = [
  { id: SORT_BY.COLOR_SCALED, name: "Scaled" },
  { id: SORT_BY.COLOR_UNSCALED, name: "Unscaled" },
];

const DEFAULT_INPUT_DROPDOWN_PROPS: Partial<IInputDropdownProps> = {
  sdsStyle: "square",
};

export default function ColorScale({
  handleIsScaledChange,
}: Props): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { sortBy } = useContext(StateContext);

  const colorScaledOptionLabel = useMemo(() => {
    return (
      COLOR_SCALE_OPTIONS.find((option) => option.id === sortBy.scaled)?.name ||
      COLOR_SCALE_OPTIONS[0].name
    );
  }, [sortBy]);

  return (
    <Wrapper>
      <LabelWrapper>
        <Label>Color Scale</Label>

        <Tooltip title="Expression is scaled to the range [0,1]. Scaling is done by assigning the minimum value in the current view to 0 and the max is assigned to 1.">
          <FlexDiv>
            <Icon sdsIcon="infoCircle" sdsSize="s" sdsType="static" />
          </FlexDiv>
        </Tooltip>
      </LabelWrapper>

      <StyledDropdown
        data-test-id="color-scale-dropdown"
        onChange={colorScaleOnChange}
        label={colorScaledOptionLabel}
        options={COLOR_SCALE_OPTIONS}
        InputDropdownProps={DEFAULT_INPUT_DROPDOWN_PROPS}
      />
    </Wrapper>
  );

  function colorScaleOnChange(
    value: { id?: SORT_BY; name: string } | null
  ): void {
    if (!dispatch || !value) return;

    handleIsScaledChange();

    dispatch(selectSortBy({ scaled: value.id as SORT_BY }));
  }
}
