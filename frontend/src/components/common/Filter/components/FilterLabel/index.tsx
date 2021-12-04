import { Button } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { CATEGORY_KEY } from "src/components/common/Filter/common/entities";
import { CategoryButton } from "./style";

interface Props {
  isDisabled: boolean;
  label: CATEGORY_KEY;
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
