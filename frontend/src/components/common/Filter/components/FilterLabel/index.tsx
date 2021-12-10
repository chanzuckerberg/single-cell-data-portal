import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { CATEGORY_LABEL } from "src/components/common/Filter/common/entities";
import { CategoryButton } from "./style";

interface Props {
  isDisabled: boolean;
  label: CATEGORY_LABEL;
}

export default function FilterLabel({ isDisabled, label }: Props): JSX.Element {
  return (
    <CategoryButton>
      <Button
        disabled={isDisabled}
        minimal
        rightIcon={IconNames.CHEVRON_DOWN}
        text={label}
      />
    </CategoryButton>
  );
}
