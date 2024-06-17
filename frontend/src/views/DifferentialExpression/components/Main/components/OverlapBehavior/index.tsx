import { ReactComponent as CirclesOverlapLeft } from "src/common/icons/circlesOverlapLeft.svg";
import { ReactComponent as CirclesOverlapRight } from "src/common/icons/circlesOverlapRight.svg";
import { ReactComponent as CirclesOverlapBoth } from "src/common/icons/circlesOverlapBoth.svg";

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
        onChange={(_, value) => {
          setActiveValue(value);
        }}
        defaultValue={"excludeTwo"}
        buttonDefinition={[
          {
            icon: <CirclesOverlapLeft />,
            value: "excludeOne",
            tooltipText: "Exclude overlapping cells from group 1",
          },
          {
            icon: <CirclesOverlapRight />,
            value: "excludeTwo",
            tooltipText: "Exclude overlapping cells from group 2",
          },
          {
            icon: <CirclesOverlapBoth />,
            value: "retainBoth",
            tooltipText: "Retain overlapping cells in both groups",
          },
        ]}
      />
    </Wrapper>
  );
}
