import CirclesOverlapLeft from "src/common/icons/circlesOverlapLeft.svg";
import CirclesOverlapRight from "src/common/icons/circlesOverlapRight.svg";
import CirclesOverlapBoth from "src/common/icons/circlesOverlapBoth.svg";

import { SegmentedControl } from "@czi-sds/components";
import { Label } from "../common/style";
import { Wrapper } from "./style";
import { useConnect } from "./connect";
import { DIFFERENTIAL_EXPRESSION_OVERLAP_BEHAVIOR } from "src/views/DifferentialExpression/common/constants";

export default function Method(): JSX.Element {
  const { activeValue, setActiveValue } = useConnect();

  return (
    <Wrapper>
      <Label>Overlap Behavior</Label>
      <SegmentedControl
        data-testid={DIFFERENTIAL_EXPRESSION_OVERLAP_BEHAVIOR}
        value={activeValue}
        onChange={(_, value) => value && setActiveValue(value)}
        defaultValue="excludeTwo"
        buttonDefinition={BUTTON_DEFINITION}
      />
    </Wrapper>
  );
}

const BUTTON_DEFINITION = [
  {
    icon: <CirclesOverlapBoth />,
    value: "retainBoth",
    tooltipText: "Retain overlapping cells in both groups",
  },
  {
    icon: <CirclesOverlapLeft />,
    value: "excludeTwo",
    tooltipText: "Remove overlapping cells from group 2",
  },
  {
    icon: <CirclesOverlapRight />,
    value: "excludeOne",
    tooltipText: "Remove overlapping cells from group 1",
  },
];
