import { Button, Checkbox, Collapse, Intent } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import React, { FC, useState } from "react";
import { GRAY } from "src/components/common/theme";
import { BulletWrapper, StyledIcon, Text, Wrapper } from "./style";

export const POLICY_BULLETS = {
  items: [
    "I give CZI permission to use, display and create derivative works (e.g. visualizations) of this Data for purposes of offering the cellxgene Public Portal, and I have the authority to give this permission.",
    "It is my responsibility to ensure that this Data is not identifiable. In particular, I commit that I will remove any direct personal identifiers in the metadata portions of the Data, and that CZI may further contact me if it believes more work is needed to de-identify it.",
    "If I choose to publish this data publicly, I understand that anyone with access to cellxgene will be able to access it subject to a CC-BY license, meaning they can download, share, and adapt the data without restriction beyond providing attribution to me.",
  ],
  version: "1.0",
};

interface Props {
  handleChange: (value: string) => void;
}

const Bullet = ({ text }: { text: string }) => {
  return (
    <BulletWrapper>
      <StyledIcon color={GRAY.A} iconSize={10} icon={IconNames.SYMBOL_CIRCLE} />
      <Text>{text}</Text>
    </BulletWrapper>
  );
};

const Policy: FC<Props> = ({ handleChange }) => {
  const [isChecked, setIsChecked] = useState(false);
  const [isOpen, setIsOpen] = useState(false);

  return (
    <Wrapper>
      <Checkbox checked={isChecked} onChange={handleChange_}>
        I agree to cellxgene&#39;s data submission policies.
        <Button minimal intent={Intent.PRIMARY} onClick={handleShowButtonClick}>
          {isOpen ? "Hide" : "Show"} Details
        </Button>
      </Checkbox>
      <Collapse isOpen={isOpen}>
        {POLICY_BULLETS.items.map((item) => (
          <Bullet key={item} text={item} />
        ))}
      </Collapse>
    </Wrapper>
  );

  function handleChange_() {
    const newIsChecked = !isChecked;

    handleChange(newIsChecked ? POLICY_BULLETS.version : "");
    setIsChecked(newIsChecked);
  }

  function handleShowButtonClick() {
    setIsOpen(!isOpen);
  }
};

export default Policy;
