import { ReactComponent as CirclesOverlapLeft } from "src/common/icons/circlesOverlapLeft.svg";
import { ReactComponent as CirclesOverlapRight } from "src/common/icons/circlesOverlapRight.svg";

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
            icon: <CirclesOverlapLeft data-element="excludeOne" />,
            value: "excludeOne",
            tooltipText: "Exclude overlapping cells from group 1",
          },
          {
            icon: <CirclesOverlapRight data-element="excludeTwo" />,
            value: "excludeTwo",
            tooltipText: "Exclude overlapping cells from group 2",
          },
          {
            icon: "CirclesOverlap2",
            value: "excludeBoth",
            tooltipText: "Exclude overlapping cells from both groups",
          },
        ]}
      />
    </Wrapper>
  );
}
